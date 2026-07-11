#!/usr/bin/env python3
"""
Genera las Fig. 19 y 20 del TFM a partir de benchmarks Snakemake de IsoAnnot.
Ejecutar: python extract_times.py 
"""

import csv
import re
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parent

# Cada entrada produce una figura PNG + TSV resumen
FIGURAS = [
    {
        "bench_root": BASE / "referencias_tfm/_test_figures/benchmarks_antes",
        "out_png": BASE / "referencias_tfm/_test_figures/figura_19_benchmark_antes.png",
        "title": (
            "Rendimiento computacional de IsoAnnot antes de las optimizaciones "
            "de NLS y sitios de union a miRNA"
        ),
    },
    {
        "bench_root": BASE / "referencias_tfm/_test_figures/benchmarks",
        "out_png": BASE / "referencias_tfm/_test_figures/figura_20_benchmark_despues.png",
        "title": (
            "Rendimiento computacional de IsoAnnot tras la paralelizacion de NucImport "
            "y el mapeo acelerado de sitios de union a miRNA"
        ),
    },
]

TOP_N = 20
NUCIMPORT_RE = re.compile(r"^run_nucimport_(\d+)$")

FONT_AXIS = 16
FONT_TITLE = 18
FONT_TICK = 13
FONT_BAR = 12


def leer_benchmark(path):
    with open(path, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row:
                return {
                    "s": float(row.get("s") or 0),
                    "max_rss": float(row.get("max_rss") or 0),
                    "cpu_time": float(row.get("cpu_time") or 0),
                }
    return None


def cargar_reglas(bench_root):
    if not bench_root.is_dir():
        sys.exit(f"No existe el directorio: {bench_root}")

    filas = []
    for tsv in sorted(bench_root.glob("**/*.tsv")):
        datos = leer_benchmark(tsv)
        if datos:
            filas.append({"rule": tsv.stem, **datos})

    # Colapsar chunks NucImport: max tiempo/RAM, suma CPU
    chunks = [r for r in filas if NUCIMPORT_RE.match(r["rule"])]
    otras = [r for r in filas if not NUCIMPORT_RE.match(r["rule"])]
    if chunks:
        otras.append({
            "rule": f"run_nucimport (max. de {len(chunks)} fragmentos)",
            "s": max(r["s"] for r in chunks),
            "max_rss": max(r["max_rss"] for r in chunks),
            "cpu_time": sum(r["cpu_time"] for r in chunks),
        })

    otras.sort(key=lambda r: r["s"], reverse=True)
    return otras


def resumir(filas):
    out = []
    for r in filas:
        out.append({
            "rule": r["rule"],
            "time_min": round(r["s"] / 60, 2),
            "cpu_time_min": round(r["cpu_time"] / 60, 2),
            "max_memory_gb": round(r["max_rss"] / 1024, 2),
        })
    return out


def guardar_tsv(filas, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["regla", "tiempo_min", "tiempo_cpu_min", "ram_max_gb"])
        for r in filas:
            w.writerow([r["rule"], r["time_min"], r["cpu_time_min"], r["max_memory_gb"]])


def graficar(filas, out_path, titulo, top_n):
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError:
        print("matplotlib/seaborn no instalados; se omite figura.")
        return

    top = [r for r in filas if r["time_min"] is not None][:top_n]
    if not top:
        print(f"[AVISO] Sin datos de tiempo en {out_path}")
        return

    reglas = [r["rule"] for r in top]
    tiempos = [r["time_min"] for r in top]
    ram = [r["max_memory_gb"] or 0.0 for r in top]

    sns.set_theme(style="whitegrid", font_scale=1.25)
    fig, axes = plt.subplots(1, 2, figsize=(16, max(6, 0.38 * len(top))), sharey=True)

    bars_t = sns.barplot(x=tiempos, y=reglas, ax=axes[0], color="steelblue")
    axes[0].set_xlabel("Tiempo de ejecucion (min)", fontsize=FONT_AXIS, fontweight="bold")
    axes[0].set_ylabel(
        f"Regla del pipeline (top {len(top)} mas lentas)",
        fontsize=FONT_AXIS, fontweight="bold",
    )
    axes[0].set_title("A. Cuello de botella de tiempo de ejecucion", fontsize=FONT_TITLE)
    axes[0].tick_params(labelsize=FONT_TICK)

    max_t = max(tiempos) or 1.0
    for patch in bars_t.patches:
        w = patch.get_width()
        if w > 0:
            axes[0].text(
                w + max_t * 0.02, patch.get_y() + patch.get_height() / 2,
                f"{int(w)}" if w >= 10 else f"{w:.1f}",
                ha="left", va="center", fontsize=FONT_BAR,
            )
    axes[0].set_xlim(0, max_t * 1.15)

    bars_m = sns.barplot(x=ram, y=reglas, ax=axes[1], color="indianred")
    axes[1].set_xlabel("Memoria RAM maxima (GB)", fontsize=FONT_AXIS, fontweight="bold")
    axes[1].set_title("B. Cuello de botella de uso de memoria", fontsize=FONT_TITLE)
    axes[1].tick_params(labelsize=FONT_TICK)

    max_m = max(ram) or 1.0
    for patch in bars_m.patches:
        w = patch.get_width()
        if w > 1.0:
            axes[1].text(
                w + max_m * 0.02, patch.get_y() + patch.get_height() / 2,
                f"{w:.1f}", ha="left", va="center", fontsize=FONT_BAR,
            )
    axes[1].set_xlim(0, max(max_m * 1.15, 1.0))

    fig.suptitle(titulo, fontsize=FONT_TITLE, y=1.02)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Figura guardada: {out_path}")


def main():
    for cfg in FIGURAS:
        print(f"\n=== {cfg['out_png'].name} ===")
        filas = resumir(cargar_reglas(cfg["bench_root"]))
        for r in filas[:TOP_N]:
            print(f"{r['rule']:<45} {r['time_min']:>8.1f} min  {r['max_memory_gb']:>6.1f} GB")
        guardar_tsv(filas, cfg["out_png"].with_suffix(".tsv"))
        graficar(filas, cfg["out_png"], cfg["title"], TOP_N)


if __name__ == "__main__":
    main()
