#!/usr/bin/env python3
"""
Estadisticas de novedad SQANTI (Tabla 6, Fig. 21-24). Ejecutar:
python novelty_category_stats.py  (ajustar OUTPUT_ROOT y RUNS al inicio).
"""

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

BASE = Path(__file__).resolve().parent
OUT_DIR = BASE / "tfm/stats_4_4_1_novelty"
GENERAR_FIGURAS = True
INCLUIR_MIRNA_BS = True

OUTPUT_ROOT = Path("/home/rlopez/tfm")
DB = "mytranscripts"

RUNS = [
    {"prefix": "Dmelanogaster", "species_name": "fly",       "species": "Drosophila melanogaster"},
    {"prefix": "Athaliana",     "species_name": "athaliana", "species": "Arabidopsis thaliana"},
    {"prefix": "Mmusculus",     "species_name": "mouse",     "species": "Mus musculus"},
    {"prefix": "Hsapiens",      "species_name": "human",     "species": "Homo sapiens"},
]

CATEGORIAS = ("Novel_gene", "Novel_known_gene", "Known")
ETIQUETAS = {
    "Known": "Known isoforms",
    "Novel_known_gene": "Novel isoforms from known genes",
    "Novel_gene": "Novel gene isoforms",
}

PROTEIN_FEATURES = (
    "ACT_SITE", "BINDING", "COILED", "COMPBIAS", "DISORDER", "DOMAIN",
    "INTRAMEM", "MOTIF", "NLS", "PTM", "SIGNAL", "TRANSMEM",
)
TRANSCRIPT_FEATURES_BASE = ("3UTRmotif", "5UTRmotif", "NMD", "PAS", "repeat", "uORF")
MIRNA_BS = "miRNA_bs"

# Tipografia figuras TFM
FONT_RC = {
    "font.size": 14, "axes.titlesize": 16, "axes.labelsize": 14,
    "xtick.labelsize": 12, "ytick.labelsize": 12, "legend.fontsize": 14,
    "figure.titlesize": 18,
}
LEGEND_Y = 0.90
SUBPLOTS_TOP = 0.78


def rutas(prefix, species_name):
    return (
        OUTPUT_ROOT / "data" / prefix / "config" / DB / "sqanti_classification.txt",
        OUTPUT_ROOT / "data" / prefix / f"{species_name}_tappas_{DB}_annotation_file.gff3_mod",
    )


def clasificar_novedad(gene, transcript):
    gene = (gene or "").strip()
    tx = (transcript or "").strip()
    if tx and tx.lower() != "novel" and tx.upper() != "NA":
        return "Known"
    if not gene or gene.lower() in ("novel", "na") or gene.lower().startswith("novelgene"):
        return "Novel_gene"
    return "Novel_known_gene"


def leer_sqanti(path):
    cat, coding = {}, {}
    with open(path, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            iso = row["isoform"].strip()
            if not iso:
                continue
            cat[iso] = clasificar_novedad(
                row.get("associated_gene", ""), row.get("associated_transcript", "")
            )
            coding[iso] = row.get("coding", "").strip().lower() == "coding"
    return cat, coding


def mapear_feature(source, feature):
    feat = feature.strip()
    if feat in ("3UTRmotif", "3'UTRmotif", "5UTRmotif", "NMD", "PAS", "uORF", "repeat"):
        return ("transcript", feat.replace("3'UTRmotif", "3UTRmotif"))
    if feat in PROTEIN_FEATURES:
        return ("protein", feat)
    if INCLUIR_MIRNA_BS and feat in ("miRNA_binding_site", "mirna_binding_site"):
        return ("transcript", MIRNA_BS)
    return None


def leer_gff3(path, validos):
    hits = defaultdict(set)
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            campos = line.rstrip("\n").split("\t")
            if len(campos) < 3:
                continue
            iso = campos[0].strip()
            if iso not in validos:
                continue
            bucket = mapear_feature(campos[1], campos[2])
            if bucket:
                hits[bucket].add(iso)
    return hits


def pct(n, total):
    return round(100.0 * n / total, 2) if total else 0.0


def features_transcript():
    if INCLUIR_MIRNA_BS:
        return TRANSCRIPT_FEATURES_BASE + (MIRNA_BS,)
    return TRANSCRIPT_FEATURES_BASE


def tabla_categorias(cat, coding):
    filas = []
    for c in CATEGORIAS:
        isos = [i for i, k in cat.items() if k == c]
        n_cod = sum(1 for i in isos if coding.get(i))
        filas.append({
            "category": c, "category_label": ETIQUETAS[c],
            "n_transcripts": len(isos), "n_coding": n_cod,
            "pct_coding": pct(n_cod, len(isos)),
        })
    return filas


def tabla_cobertura(cat, hits):
    filas = []
    feats_tx = features_transcript()
    for nivel, features in (("protein", PROTEIN_FEATURES), ("transcript", feats_tx)):
        for feat in features:
            anotados = hits.get((nivel, feat), set())
            for c in CATEGORIAS:
                isos_cat = {i for i, k in cat.items() if k == c}
                n_ann = len(anotados & isos_cat)
                filas.append({
                    "level": nivel, "feature": feat, "category": c,
                    "category_label": ETIQUETAS[c],
                    "n_annotated": n_ann, "n_total": len(isos_cat),
                    "pct_annotated": pct(n_ann, len(isos_cat)),
                })
    return filas


def escribir_tsv(path, filas, columnas):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=columnas, extrasaction="ignore")
        w.writeheader()
        w.writerows(filas)


def italic_label(text):
    text_math = text.replace(" ", r"\ ")
    return f"$\\it{{{text_math}}}$"


def graficar(species, cobertura, path):
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        return

    plt.rcParams.update(FONT_RC)
    colores = {"Known": "#4DBBD5", "Novel_known_gene": "#00A087", "Novel_gene": "#E64B35"}
    feats_tx = features_transcript()

    fig, axes = plt.subplots(1, 2, figsize=(16, 6.5), sharey=True)
    paneles = [("protein", PROTEIN_FEATURES, axes[0]), ("transcript", feats_tx, axes[1])]
    bar_w = 0.22

    for panel_idx, (nivel, features, ax) in enumerate(paneles):
        x = np.arange(len(features))
        for j, cat in enumerate(CATEGORIAS):
            vals = []
            for feat in features:
                match = [r for r in cobertura
                         if r["level"] == nivel and r["feature"] == feat and r["category"] == cat]
                vals.append(match[0]["pct_annotated"] if match else 0.0)
            offset = (j - 1) * bar_w
            ax.bar(x + offset, vals, bar_w,
                   label=italic_label(ETIQUETAS[cat]) if panel_idx == 0 else "",
                   color=colores[cat])
        ax.set_title("Proteina" if nivel == "protein" else "Transcrito")
        ax.set_ylabel("% Isoformas anotadas")
        ax.set_xticks(x)
        ax.set_xticklabels(features, rotation=45, ha="right")
        ax.set_ylim(0, 100)

    fig.tight_layout(rect=[0, 0, 1, SUBPLOTS_TOP])
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center",
               bbox_to_anchor=(0.5, LEGEND_Y), ncol=3, frameon=True)
    sp = species.replace(" ", r"\ ")
    fig.suptitle(
        f"Proporcion de isoformas anotadas a nivel de transcrito y proteina "
        f"$\\it{{({sp})}}$",
        y=0.98,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Figura guardada: {path}")


def main():
    filas_cat, filas_cov, todos = [], [], []

    for spec in RUNS:
        cpath, gpath = rutas(spec["prefix"], spec["species_name"])
        species = spec["species"]
        if not cpath.is_file() or not gpath.is_file():
            print(f"[AVISO] Omitido {species}: faltan archivos", file=sys.stderr)
            continue

        cat, coding = leer_sqanti(cpath)
        hits = leer_gff3(gpath, set(cat.keys()))
        cat_table = tabla_categorias(cat, coding)
        cobertura = tabla_cobertura(cat, hits)

        for r in cat_table:
            filas_cat.append({"species": species, **r})
        for r in cobertura:
            filas_cov.append({"species": species, **r})
        todos.append({"species": species, "category_table": cat_table, "coverage": cobertura})

        print(f"\n=== {species} ===")
        for r in cat_table:
            print(f"  {r['category_label']:<40} {r['n_transcripts']:>8,}  {r['pct_coding']:>5.1f}%"
                  .replace(",", "."))

        if GENERAR_FIGURAS:
            safe = species.replace(" ", "_")
            graficar(species, cobertura, OUT_DIR / f"fig_novelty_coverage_{safe}.png")

    if not todos:
        sys.exit("No se calculo ningun run. Ajusta OUTPUT_ROOT en el script.")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    escribir_tsv(OUT_DIR / "tabla_6_novelty_transcripts_by_category.tsv", filas_cat,
                 ["species", "category", "category_label",
                  "n_transcripts", "n_coding", "pct_coding"])
    escribir_tsv(OUT_DIR / "tabla_6_novelty_annotation_coverage_by_category.tsv", filas_cov,
                 ["species", "level", "feature", "category",
                  "n_annotated", "n_total", "pct_annotated"])
    with open(OUT_DIR / "novelty_category_stats.json", "w", encoding="utf-8") as f:
        json.dump(todos, f, indent=2, ensure_ascii=False)

    print(f"\nResultados en: {OUT_DIR.resolve()}")


if __name__ == "__main__":
    main()
