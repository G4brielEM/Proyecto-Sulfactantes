import csv
import json
import time
from collections import Counter
from pathlib import Path

import requests
from rdkit import Chem

# CONFIGURACIÓN
CARBONOS_MIN = 10
CARBONOS_MAX = 30
BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

CARPETA = Path(__file__).resolve().parent
CACHE = CARPETA / "cache_pubchem"
CACHE.mkdir(exist_ok=True)

# Cada patrón contiene un éster dentro de un anillo.
# Los demás átomos del anillo son carbonos.
# "~" permite distintos tipos de enlace entre esos carbonos.
PATRONES = {
    "beta": "[C;X3]1(=[O;X1])-[O;X2]-[#6]~[#6]~1",
    "gamma": "[C;X3]1(=[O;X1])-[O;X2]-[#6]~[#6]~[#6]~1",
    "delta": "[C;X3]1(=[O;X1])-[O;X2]-[#6]~[#6]~[#6]~[#6]~1",
    "epsilon": "[C;X3]1(=[O;X1])-[O;X2]-[#6]~[#6]~[#6]~[#6]~[#6]~1",
}

CONSULTAS = {
    familia: Chem.MolFromSmarts(patron)
    for familia, patron in PATRONES.items()
}

SESION = requests.Session()


def consultar(ruta, parametros=None):
    """Consultar PubChem con pausas y reintentos."""
    ultimo_error = None

    for intento in range(5):
        time.sleep(0.4)

        try:
            respuesta = SESION.get(
                BASE + ruta,
                params=parametros,
                timeout=(20, 120),
            )

            if respuesta.status_code in (429, 500, 502, 503, 504):
                raise requests.RequestException(
                    f"PubChem respondió HTTP {respuesta.status_code}"
                )

            respuesta.raise_for_status()
            return respuesta.json()

        except (requests.RequestException, ValueError) as error:
            ultimo_error = error
            print(
                f"  Intento {intento + 1}/5 falló: {error}",
                flush=True,
            )
            time.sleep(2 ** intento)

    raise RuntimeError(f"No se pudo consultar PubChem: {ultimo_error}")


def buscar_familia(familia, patron):
    """Guardar los identificadores para poder reutilizarlos."""
    archivo = CACHE / f"cids_{familia}_v1.json"

    if archivo.exists():
        print(f"Usando búsqueda guardada de {familia}.", flush=True)
        return json.loads(archivo.read_text(encoding="utf-8"))

    print(f"Buscando lactonas {familia} en PubChem...", flush=True)

    datos = consultar(
        "/compound/fastsubstructure/smarts/cids/JSON",
        {"smarts": patron},
    )

    cids = datos["IdentifierList"]["CID"]

    # El límite documentado por defecto es de dos millones.
    # Si se alcanza, no tratamos el resultado como completo.
    if len(cids) >= 2_000_000:
        raise RuntimeError(
            f"La búsqueda de {familia} alcanzó el límite de resultados."
        )

    archivo.write_text(json.dumps(cids), encoding="utf-8")
    return cids


def descargar_lote(cids):
    """Descargar propiedades de hasta 100 registros."""
    identificadores = ",".join(str(cid) for cid in cids)
    archivo = CACHE / f"propiedades_{cids[0]}_{cids[-1]}.json"

    if archivo.exists():
        guardado = json.loads(archivo.read_text(encoding="utf-8"))
        if guardado["cids"] == cids:
            return guardado["propiedades"]

    datos = consultar(
        f"/compound/cid/{identificadores}"
        "/property/Title,SMILES,MolecularFormula,InChIKey/JSON"
    )

    propiedades = datos["PropertyTable"]["Properties"]

    recibidos = {dato["CID"] for dato in propiedades}
    if recibidos != set(cids):
        raise RuntimeError("PubChem devolvió un lote incompleto.")

    archivo.write_text(
        json.dumps({"cids": cids, "propiedades": propiedades}),
        encoding="utf-8",
    )
    return propiedades


def main():
    todos_cids = set()

    for familia, patron in PATRONES.items():
        cids = buscar_familia(familia, patron)
        print(f"  {familia}: {len(cids)} registros encontrados.", flush=True)
        todos_cids.update(cids)

    cids = sorted(todos_cids)
    print(f"\nRegistros únicos por revisar: {len(cids)}", flush=True)

    columnas = [
        "cid",
        "nombre",
        "familias",
        "carbonos_totales",
        "formula",
        "smiles_pubchem",
        "smiles_rdkit",
        "inchikey",
        "numero_fragmentos",
        "centros_quirales_sin_especificar",
        "fuente",
    ]

    salida = CARPETA / "lactonas_10a30.csv"
    errores_archivo = CARPETA / "errores_descarga.csv"

    vistos = set()
    conteos = Counter()
    guardadas = 0
    errores = 0

    with (
        salida.open("w", newline="", encoding="utf-8-sig") as archivo,
        errores_archivo.open(
            "w", newline="", encoding="utf-8-sig"
        ) as archivo_errores,
    ):
        escritor = csv.DictWriter(archivo, fieldnames=columnas)
        escritor.writeheader()

        escritor_errores = csv.writer(archivo_errores)
        escritor_errores.writerow(["cids", "error"])

        for inicio in range(0, len(cids), 100):
            lote = cids[inicio:inicio + 100]

            try:
                propiedades = descargar_lote(lote)
            except Exception as error:
                errores += len(lote)
                escritor_errores.writerow([
                    ";".join(map(str, lote)), str(error)
                ])
                archivo_errores.flush()
                print(f"  Lote pendiente por error: {error}", flush=True)
                continue

            for dato in propiedades:
                try:
                    smiles = dato["SMILES"]
                    mol = Chem.MolFromSmiles(smiles)

                    if mol is None:
                        raise ValueError("SMILES no interpretable por RDKit.")

                    carbonos = sum(
                        atomo.GetAtomicNum() == 6
                        for atomo in mol.GetAtoms()
                    )

                    if not CARBONOS_MIN <= carbonos <= CARBONOS_MAX:
                        continue

                    familias = [
                        familia
                        for familia, consulta in CONSULTAS.items()
                        if mol.HasSubstructMatch(consulta)
                    ]

                    if not familias:
                        continue

                    # Conserva la estereoquímica especificada.
                    smiles_normalizado = Chem.MolToSmiles(
                        mol, canonical=True, isomericSmiles=True
                    )

                    if smiles_normalizado in vistos:
                        continue

                    vistos.add(smiles_normalizado)

                    centros = Chem.FindMolChiralCenters(
                        mol, includeUnassigned=True
                    )

                    escritor.writerow({
                        "cid": dato["CID"],
                        "nombre": dato.get("Title", ""),
                        "familias": ";".join(familias),
                        "carbonos_totales": carbonos,
                        "formula": dato.get("MolecularFormula", ""),
                        "smiles_pubchem": smiles,
                        "smiles_rdkit": smiles_normalizado,
                        "inchikey": dato.get("InChIKey", ""),
                        "numero_fragmentos": len(Chem.GetMolFrags(mol)),
                        "centros_quirales_sin_especificar": sum(
                            etiqueta == "?" for _, etiqueta in centros
                        ),
                        "fuente": (
                            "https://pubchem.ncbi.nlm.nih.gov/compound/"
                            f"{dato['CID']}"
                        ),
                    })

                    guardadas += 1
                    conteos.update(familias)

                except Exception as error:
                    errores += 1
                    escritor_errores.writerow([
                        dato.get("CID", ""), str(error)
                    ])

            archivo.flush()
            archivo_errores.flush()

            revisados = min(inicio + 100, len(cids))
            print(
                f"Revisados {revisados}/{len(cids)}"
                f" | Moléculas guardadas: {guardadas}",
                flush=True,
            )

    print(f"\nMoléculas únicas guardadas: {guardadas}")
    for familia in PATRONES:
        print(f"  {familia}: {conteos[familia]}")

    print("Una molécula puede pertenecer a varias familias.")
    print(f"CSV: {salida}")

    if errores:
        print(f"ATENCIÓN: quedaron {errores} registros pendientes.")
        print(f"Revisa: {errores_archivo}")
        print("El dataset es parcial. Puedes volver a ejecutar el programa.")
    else:
        print("Se procesaron todos los registros recuperados.")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"\nEl programa se detuvo: {error}")
        print("Los datos descargados en cache_pubchem se conservan.")