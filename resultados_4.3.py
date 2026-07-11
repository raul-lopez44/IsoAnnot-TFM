#!/usr/bin/env python3
"""
Calcula las Tablas 4 y 5 del TFM (transcriptomas de referencia Ensembl/RefSeq).
Incluye anotacion a nivel transcrito (UTRscan, RepeatMasker, NMD, miRNA binding sites)
y a nivel proteina (InterProScan, UniProt, GO, Reactome, NLS).
Ejecutar: python resultados_4.3.py  (ajustar OUTPUT_ROOT al inicio).
"""

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

BASE = Path(__file__).resolve().parent
OUT_DIR = BASE / "referencias_tfm/stats_4_3"

# Directorio raiz con data/{prefix}/config/{db}/ y output/{db}/layers/
OUTPUT_ROOT = Path("/home/rlopez/tfm_resultados")

NOMBRES_ESPECIE = {
    "Hsapiens": "Homo sapiens",
    "Mmusculus": "Mus musculus",
    "Drerio": "Danio rerio",
    "Dmelanogaster": "Drosophila melanogaster",
    "Athaliana": "Arabidopsis thaliana",
}

# 5 especies x ensembl/refseq 
RUNS = [
    {"prefix": p, "reference": db}
    for p in NOMBRES_ESPECIE
    for db in ("ensembl", "refseq")
]

CAPAS_TRANSCRITO = {
    "layer_utrscan", "layer_repeatmasker", "layer_nmd", "layer_mirna_bs",
}
CAPAS_PROTEINA = {
    "layer_interproscan", "layer_uniprot", "layer_go",
    "layer_reactome", "layer_nls",
}
# Se excluyen layer_exons y layer_junctions (estructura genomica, no funcional)
CAPAS_INCLUIDAS = CAPAS_TRANSCRITO | CAPAS_PROTEINA


def rutas(prefix, reference):
    return (
        OUTPUT_ROOT / "data" / prefix / "config" / reference / "sqanti_classification.txt",
        OUTPUT_ROOT / "data" / prefix / "output" / reference / "layers",
    )


def leer_classification(path):
    todos, coding = set(), set()
    with open(path, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            iso = row["isoform"].strip()
            if not iso:
                continue
            todos.add(iso)
            if row.get("coding", "").strip().lower() == "coding":
                coding.add(iso)
    return todos, coding


def leer_capas(layers_dir, isoformas_validos):
    transcrito, proteina = set(), set()
    por_capa = defaultdict(set)

    gtfs = sorted(layers_dir.glob("layer_*.gtf")) + sorted(layers_dir.glob("layer_*.gtf.gz"))
    for gtf in gtfs:
        capa = gtf.name.replace(".gtf.gz", "").replace(".gtf", "")
        if capa not in CAPAS_INCLUIDAS:
            continue
        nivel = "transcrito" if capa in CAPAS_TRANSCRITO else "proteina"
        opener = open
        if str(gtf).endswith(".gz"):
            import gzip
            opener = lambda p: gzip.open(p, "rt", encoding="utf-8")
        with opener(gtf) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                iso = line.split("\t", 1)[0].strip()
                if iso in isoformas_validos:
                    por_capa[capa].add(iso)
                    (transcrito if nivel == "transcrito" else proteina).add(iso)
    return transcrito, proteina, por_capa


def pct(n, total):
    return round(100.0 * n / total, 2) if total else 0.0


def analizar_run(prefix, reference):
    cpath, ldir = rutas(prefix, reference)
    todos, coding = leer_classification(cpath)
    transcrito, proteina, por_capa = leer_capas(ldir, todos)
    total = transcrito | proteina
    return {
        "species": NOMBRES_ESPECIE[prefix],
        "reference": reference,
        "n_transcripts": len(todos),
        "n_coding": len(coding),
        "pct_coding": pct(len(coding), len(todos)),
        "pct_transcript_level": pct(len(transcrito), len(todos)),
        "pct_protein_level": pct(len(proteina & coding), len(coding)),
        "pct_total": pct(len(total), len(todos)),
        "per_layer_counts": {k: len(v) for k, v in sorted(por_capa.items())},
    }


def escribir_tsv(path, filas, columnas):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=columnas, extrasaction="ignore")
        w.writeheader()
        w.writerows(filas)


def main():
    resultados = []
    for spec in RUNS:
        cpath, ldir = rutas(spec["prefix"], spec["reference"])
        if not cpath.is_file():
            print(f"[AVISO] Sin classification: {cpath}", file=sys.stderr)
            continue
        if not ldir.is_dir():
            print(f"[AVISO] Sin capas: {ldir}", file=sys.stderr)
            continue
        resultados.append(analizar_run(spec["prefix"], spec["reference"]))

    if not resultados:
        sys.exit("No se calculo ningun run. Ajusta OUTPUT_ROOT en el script.")

    print("\n=== Tabla 4: Transcritos y % codificantes ===")
    for r in resultados:
        print(f"{r['species']:<22} {r['reference']:<8} {r['n_transcripts']:>10,} {r['pct_coding']:>8.2f}%"
              .replace(",", "."))

    print("\n=== Tabla 5: Cobertura funcional ===")
    for r in resultados:
        print(f"{r['species']:<22} {r['reference']:<8} "
              f"transcrito {r['pct_transcript_level']:>6.2f}%  "
              f"proteina {r['pct_protein_level']:>6.2f}%  "
              f"total {r['pct_total']:>6.2f}%")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    escribir_tsv(OUT_DIR / "tabla_4_transcripts_coding.tsv", resultados,
                 ["species", "reference", "n_transcripts", "n_coding", "pct_coding"])
    escribir_tsv(OUT_DIR / "tabla_5_annotation_coverage.tsv", resultados,
                 ["species", "reference", "pct_transcript_level",
                  "pct_protein_level", "pct_total"])
    with open(OUT_DIR / "reference_stats.json", "w", encoding="utf-8") as f:
        json.dump(resultados, f, indent=2, ensure_ascii=False)

    print(f"\nTablas en: {OUT_DIR.resolve()}")


if __name__ == "__main__":
    main()
