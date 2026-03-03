"""
================================================================
Day 18 — SNP & Variant Analysis (REAL DATA)
Author  : Subhadip Jana
Organism: Escherichia coli K-12 MG1655 (NC_000913.3)
Dataset : 3,800 SNPs | VCFv4.2 format

Analyses:
  1.  VCF parsing — full INFO/FORMAT field extraction
  2.  Transition / Transversion classification & Ts/Tv ratio
  3.  SNP density along genome (sliding window)
  4.  Allele frequency (AF) spectrum — site frequency spectrum
  5.  Read depth (DP) distribution
  6.  Genotype calls (hom-ref / het / hom-alt) breakdown
  7.  Base substitution spectrum (6 mutation classes)
  8.  Quality score (QUAL) distribution
  9.  Genomic context — SNP clustering detection
  10. Functional impact prediction (synonymous vs non-synonymous
      proxy from codon position)
================================================================
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from collections import Counter, defaultdict
import warnings
warnings.filterwarnings("ignore")

np.random.seed(42)

# ─────────────────────────────────────────────────────────────
# SECTION 1: PARSE VCF
# ─────────────────────────────────────────────────────────────

print("🔬 Parsing VCF — E. coli K-12 MG1655...")

records = []
meta    = {}

with open("data/ecoli_k12_variants.vcf") as f:
    for line in f:
        line = line.rstrip()
        if line.startswith("##"):
            if "=" in line:
                k = line[2:line.index("=")]
                meta[k] = meta.get(k, []) + [line]
            continue
        if line.startswith("#CHROM"):
            cols = line[1:].split("\t")
            continue
        parts = line.split("\t")
        if len(parts) < 8:
            continue

        chrom, pos, snp_id, ref, alt, qual, flt, info_str = parts[:8]
        fmt_str  = parts[8] if len(parts) > 8 else ""
        samp_str = parts[9] if len(parts) > 9 else ""

        # Parse INFO
        info = {}
        for item in info_str.split(";"):
            if "=" in item:
                k, v = item.split("=", 1)
                info[k] = v
            else:
                info[item] = True

        # Parse FORMAT/SAMPLE
        fmt_keys = fmt_str.split(":")
        samp_vals= samp_str.split(":")
        sample   = dict(zip(fmt_keys, samp_vals))

        records.append({
            "CHROM"  : chrom,
            "POS"    : int(pos),
            "ID"     : snp_id,
            "REF"    : ref.upper(),
            "ALT"    : alt.upper(),
            "QUAL"   : float(qual) if qual not in (".", "") else np.nan,
            "FILTER" : flt,
            "DP"     : int(info.get("DP", 0)),
            "AF"     : float(info.get("AF", 0)),
            "TYPE"   : info.get("TYPE", "SNP"),
            "GT"     : sample.get("GT", "./."),
            "SAMP_DP": int(sample.get("DP", 0)),
            "GQ"     : int(sample.get("GQ", 0)),
        })

df = pd.DataFrame(records)
L_GENOME = 4_641_652

print(f"✅ Parsed {len(df):,} variants")
print(f"   PASS filter: {(df['FILTER']=='PASS').sum():,}")
print(f"   AF range   : {df['AF'].min():.3f} – {df['AF'].max():.3f}")
print(f"   DP range   : {df['DP'].min()} – {df['DP'].max()}")

# ─────────────────────────────────────────────────────────────
# SECTION 2: TRANSITION / TRANSVERSION CLASSIFICATION
# ─────────────────────────────────────────────────────────────

print("\n🔄 Classifying Ts/Tv...")

PURINES     = {"A", "G"}
PYRIMIDINES = {"C", "T"}

# 6 canonical mutation classes (both strands collapsed)
MUT_CLASSES = {
    ("C","A"):"C>A", ("G","T"):"C>A",
    ("C","G"):"C>G", ("G","C"):"C>G",
    ("C","T"):"C>T", ("G","A"):"C>T",
    ("T","A"):"T>A", ("A","T"):"T>A",
    ("T","C"):"T>C", ("A","G"):"T>C",
    ("T","G"):"T>G", ("A","C"):"T>G",
}

def classify_snp(ref, alt):
    if ref in PURINES and alt in PURINES:   return "Transition"
    if ref in PYRIMIDINES and alt in PYRIMIDINES: return "Transition"
    return "Transversion"

def mut_class(ref, alt):
    return MUT_CLASSES.get((ref, alt), f"{ref}>{alt}")

def genotype_call(gt):
    if gt in ("0/0", "0|0"):   return "Hom-Ref"
    if gt in ("1/1", "1|1"):   return "Hom-Alt"
    if gt in ("0/1", "1/0", "0|1", "1|0"): return "Het"
    return "Unknown"

df["VAR_TYPE"]  = df.apply(lambda r: classify_snp(r["REF"], r["ALT"]), axis=1)
df["MUT_CLASS"] = df.apply(lambda r: mut_class(r["REF"], r["ALT"]), axis=1)
df["GT_CALL"]   = df["GT"].apply(genotype_call)

ts = (df["VAR_TYPE"] == "Transition").sum()
tv = (df["VAR_TYPE"] == "Transversion").sum()
tstv = ts / tv if tv > 0 else np.nan

print(f"   Transitions  : {ts:,} ({ts/len(df)*100:.1f}%)")
print(f"   Transversions: {tv:,} ({tv/len(df)*100:.1f}%)")
print(f"   Ts/Tv ratio  : {tstv:.3f}  (expected ~2.0–2.5 for genomic SNPs)")

# ─────────────────────────────────────────────────────────────
# SECTION 3: SNP DENSITY SLIDING WINDOW
# ─────────────────────────────────────────────────────────────

print("\n📈 Computing SNP density (50 kb windows)...")

WIN  = 50_000
STEP = 25_000
dens_vals, dens_pos = [], []

for start in range(0, L_GENOME - WIN, STEP):
    n_snps = ((df["POS"] >= start) & (df["POS"] < start+WIN)).sum()
    dens_vals.append(n_snps / WIN * 1000)   # SNPs per kb
    dens_pos.append((start + WIN//2) / 1e6)

dens_vals = np.array(dens_vals)
dens_pos  = np.array(dens_pos)
print(f"   Mean density : {dens_vals.mean():.3f} SNPs/kb")
print(f"   Peak density : {dens_vals.max():.3f} SNPs/kb "
      f"at {dens_pos[dens_vals.argmax()]:.2f} Mb")

# ─────────────────────────────────────────────────────────────
# SECTION 4: ALLELE FREQUENCY SPECTRUM
# ─────────────────────────────────────────────────────────────

print("\n📊 Computing allele frequency spectrum (SFS)...")

# Classify AF into population genetics categories
df["AF_CAT"] = pd.cut(df["AF"],
    bins=[0, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0],
    labels=["Rare\n(<10%)", "Low\n(10–20%)", "Intermediate\n(20–30%)",
            "Common\n(30–50%)", "High\n(50–80%)", "Fixed\n(>80%)"])

af_counts = df["AF_CAT"].value_counts().sort_index()
print("   Allele frequency categories:")
for cat, n in af_counts.items():
    print(f"   {str(cat):20s}: {n:,} ({n/len(df)*100:.1f}%)")

# ─────────────────────────────────────────────────────────────
# SECTION 5: CLUSTERING ANALYSIS
# ─────────────────────────────────────────────────────────────

print("\n🔍 Detecting SNP clusters (inter-SNP distance)...")

positions  = np.sort(df["POS"].values)
inter_dist = np.diff(positions)
clustered  = (inter_dist < 100).sum()
print(f"   Median inter-SNP distance : {np.median(inter_dist):.0f} bp")
print(f"   Mean inter-SNP distance   : {np.mean(inter_dist):.0f} bp")
print(f"   Expected (uniform)        : {L_GENOME//len(df)} bp")
print(f"   Clustered SNPs (<100 bp apart): {clustered} pairs")

# ─────────────────────────────────────────────────────────────
# SECTION 6: IMPACT PREDICTION (proxy)
# ─────────────────────────────────────────────────────────────

print("\n🧬 Predicting functional impact (codon position proxy)...")

# Codon position: pos % 3 == 0 → 3rd (often synonymous)
#                 pos % 3 != 0 → 1st or 2nd (often non-synonymous)
df["CODON_POS"] = (df["POS"] % 3).map({0:"3rd (syn-likely)",
                                         1:"1st (nonsyn-likely)",
                                         2:"2nd (nonsyn-likely)"})
# Transitions at 3rd codon position → more likely synonymous
df["IMPACT"] = np.where(
    (df["VAR_TYPE"] == "Transition") & (df["POS"] % 3 == 0),
    "Likely Synonymous",
    np.where(df["VAR_TYPE"] == "Transversion", "Likely Non-synonymous",
             "Intermediate")
)
impact_counts = df["IMPACT"].value_counts()
print("   Predicted impact:")
for k, v in impact_counts.items():
    print(f"   {k:25s}: {v:,} ({v/len(df)*100:.1f}%)")

# ─────────────────────────────────────────────────────────────
# SECTION 7: SAVE RESULTS
# ─────────────────────────────────────────────────────────────

df.to_csv("outputs/snp_annotated.csv", index=False)

summary = {
    "Total SNPs"           : len(df),
    "Transitions"          : int(ts),
    "Transversions"        : int(tv),
    "Ts/Tv ratio"          : round(tstv, 4),
    "Mean AF"              : round(df["AF"].mean(), 4),
    "Mean DP"              : round(df["DP"].mean(), 1),
    "Het calls"            : int((df["GT_CALL"]=="Het").sum()),
    "Hom-Alt calls"        : int((df["GT_CALL"]=="Hom-Alt").sum()),
    "Mean QUAL"            : round(df["QUAL"].mean(), 2),
    "SNP density (SNPs/kb)": round(dens_vals.mean(), 4),
    "Median inter-SNP (bp)": int(np.median(inter_dist)),
    "Clustered pairs"      : int(clustered),
}
pd.DataFrame(list(summary.items()),
             columns=["Metric","Value"]).to_csv("outputs/snp_summary.csv", index=False)

print("\n✅ SNP data saved → outputs/snp_annotated.csv")

# ─────────────────────────────────────────────────────────────
# SECTION 8: DASHBOARD (9 panels)
# ─────────────────────────────────────────────────────────────

print("\n🎨 Generating dashboard...")

TSTV_COLOR = {"Transition":"#E74C3C", "Transversion":"#3498DB"}
MUT_COLORS = {
    "C>A":"#3498DB","C>G":"#2C3E50","C>T":"#E74C3C",
    "T>A":"#95A5A6","T>C":"#2ECC71","T>G":"#F39C12",
}
GT_COLORS  = {"Hom-Ref":"#BDC3C7","Het":"#F39C12","Hom-Alt":"#E74C3C","Unknown":"#ECF0F1"}

fig = plt.figure(figsize=(24, 20))
fig.suptitle(
    "SNP & Variant Analysis — Escherichia coli K-12 MG1655\n"
    "3,800 SNPs | VCFv4.2 | NC_000913.3 | REAL DATA\n"
    "Ts/Tv · SFS · Density · Base spectrum · Impact prediction",
    fontsize=15, fontweight="bold", y=0.99
)

# ── Plot 1: SNP density along genome ──
ax1 = fig.add_subplot(3, 3, 1)
ax1.fill_between(dens_pos, dens_vals, alpha=0.4, color="#9B59B6")
ax1.plot(dens_pos, dens_vals, lw=1, color="#9B59B6", alpha=0.8)
ax1.axhline(dens_vals.mean(), color="#E74C3C", lw=1.5,
            linestyle="--", label=f"Mean: {dens_vals.mean():.3f} SNPs/kb")
# Mark peak
pk = dens_vals.argmax()
ax1.annotate(f"Peak\n{dens_pos[pk]:.2f} Mb",
             xy=(dens_pos[pk], dens_vals[pk]),
             xytext=(dens_pos[pk]+0.3, dens_vals[pk]-0.05),
             fontsize=8, arrowprops=dict(arrowstyle="->", color="black"))
ax1.set_xlabel("Genome position (Mb)")
ax1.set_ylabel("SNP density (SNPs/kb)")
ax1.set_title("SNP Density Along Genome\n(50 kb windows, 25 kb step)",
              fontweight="bold", fontsize=10)
ax1.legend(fontsize=8)
ax1.set_xlim(0, L_GENOME/1e6)

# ── Plot 2: Ts/Tv pie ──
ax2 = fig.add_subplot(3, 3, 2)
ts_tv_counts = df["VAR_TYPE"].value_counts()
wedge_colors = [TSTV_COLOR[k] for k in ts_tv_counts.index]
wedges, texts, autotexts = ax2.pie(
    ts_tv_counts.values, labels=ts_tv_counts.index,
    colors=wedge_colors, autopct="%1.1f%%",
    startangle=90, pctdistance=0.75,
    wedgeprops={"edgecolor":"white","linewidth":2})
for at in autotexts: at.set_fontsize(11); at.set_fontweight("bold")
ax2.set_title(f"Transition / Transversion\nTs/Tv ratio = {tstv:.3f}",
              fontweight="bold", fontsize=10)

# ── Plot 3: 6-class mutation spectrum ──
ax3 = fig.add_subplot(3, 3, 3)
mut_counts = df["MUT_CLASS"].value_counts().reindex(
    ["C>A","C>G","C>T","T>A","T>C","T>G"], fill_value=0)
bars3 = ax3.bar(mut_counts.index, mut_counts.values,
                color=[MUT_COLORS[k] for k in mut_counts.index],
                edgecolor="black", linewidth=0.5, alpha=0.9, width=0.6)
for bar, val in zip(bars3, mut_counts.values):
    ax3.text(bar.get_x()+bar.get_width()/2, bar.get_height()+10,
             f"{val:,}", ha="center", fontsize=9, fontweight="bold")
ax3.set_ylabel("Count")
ax3.set_title("6-Class Base Substitution Spectrum\n(SBS — strand-collapsed)",
              fontweight="bold", fontsize=10)
ax3.set_ylim(0, mut_counts.max()*1.18)

# ── Plot 4: Allele frequency spectrum ──
ax4 = fig.add_subplot(3, 3, 4)
ax4.hist(df["AF"], bins=40, color="#2ECC71", edgecolor="black",
         linewidth=0.4, alpha=0.85)
ax4.axvline(df["AF"].mean(), color="#E74C3C", lw=2,
            linestyle="--", label=f"Mean AF={df['AF'].mean():.3f}")
ax4.axvline(0.5, color="#F39C12", lw=1.5,
            linestyle=":", label="AF=0.5")
ax4.set_xlabel("Allele Frequency (AF)")
ax4.set_ylabel("Number of SNPs")
ax4.set_title("Allele Frequency Spectrum (SFS)\n(Site Frequency Spectrum)",
              fontweight="bold", fontsize=10)
ax4.legend(fontsize=9)

# ── Plot 5: Read depth distribution ──
ax5 = fig.add_subplot(3, 3, 5)
dp_vals = df["DP"].clip(0, 200)
ax5.hist(dp_vals, bins=40, color="#3498DB", edgecolor="black",
         linewidth=0.4, alpha=0.85)
ax5.axvline(df["DP"].mean(), color="#E74C3C", lw=2, linestyle="--",
            label=f"Mean DP={df['DP'].mean():.1f}×")
ax5.axvline(df["DP"].median(), color="#2ECC71", lw=2, linestyle=":",
            label=f"Median DP={df['DP'].median():.0f}×")
ax5.set_xlabel("Read Depth (DP)")
ax5.set_ylabel("Number of SNPs")
ax5.set_title("Read Depth Distribution\n(coverage at variant sites)",
              fontweight="bold", fontsize=10)
ax5.legend(fontsize=9)

# ── Plot 6: Genotype call breakdown ──
ax6 = fig.add_subplot(3, 3, 6)
gt_counts = df["GT_CALL"].value_counts()
gt_labels = gt_counts.index.tolist()
gt_colors_plot = [GT_COLORS.get(g,"#BDC3C7") for g in gt_labels]
wedges6, _, autotexts6 = ax6.pie(
    gt_counts.values,
    labels=[f"{g}\n(n={v:,})" for g, v in zip(gt_labels, gt_counts.values)],
    colors=gt_colors_plot,
    autopct="%1.1f%%", startangle=90,
    pctdistance=0.72,
    wedgeprops={"edgecolor":"white","linewidth":2})
for at in autotexts6: at.set_fontsize(10); at.set_fontweight("bold")
ax6.set_title("Genotype Call Distribution\n(Hom-Alt / Het / Hom-Ref)",
              fontweight="bold", fontsize=10)

# ── Plot 7: Quality score distribution ──
ax7 = fig.add_subplot(3, 3, 7)
ax7.hist(df["QUAL"].dropna(), bins=35, color="#F39C12",
         edgecolor="black", linewidth=0.4, alpha=0.85)
ax7.axvline(df["QUAL"].mean(), color="#E74C3C", lw=2, linestyle="--",
            label=f"Mean QUAL={df['QUAL'].mean():.1f}")
ax7.axvline(30, color="#9B59B6", lw=1.5, linestyle=":",
            label="QUAL=30 threshold")
ax7.set_xlabel("Variant Quality Score (QUAL)")
ax7.set_ylabel("Number of SNPs")
ax7.set_title("Variant Quality Distribution\n(higher = more confident call)",
              fontweight="bold", fontsize=10)
ax7.legend(fontsize=9)

# ── Plot 8: Inter-SNP distance distribution ──
ax8 = fig.add_subplot(3, 3, 8)
clipped_dist = inter_dist[inter_dist < np.percentile(inter_dist, 98)]
ax8.hist(clipped_dist, bins=50, color="#1ABC9C",
         edgecolor="black", linewidth=0.4, alpha=0.85)
ax8.axvline(np.median(inter_dist), color="#E74C3C", lw=2,
            linestyle="--",
            label=f"Median={np.median(inter_dist):.0f} bp")
ax8.axvline(L_GENOME//len(df), color="#F39C12", lw=1.5,
            linestyle=":", label=f"Expected={L_GENOME//len(df)} bp\n(uniform)")
ax8.set_xlabel("Inter-SNP Distance (bp)")
ax8.set_ylabel("Count")
ax8.set_title("Inter-SNP Distance Distribution\n(detects clusters & cold spots)",
              fontweight="bold", fontsize=10)
ax8.legend(fontsize=8)

# ── Plot 9: Summary table ──
ax9 = fig.add_subplot(3, 3, 9)
ax9.axis("off")
rows9 = list(summary.items())
tbl9 = ax9.table(cellText=rows9,
                  colLabels=["Metric","Value"],
                  cellLoc="left", loc="center")
tbl9.auto_set_font_size(False); tbl9.set_fontsize(9); tbl9.scale(1.8, 1.95)
for j in range(2):
    tbl9[(0,j)].set_facecolor("#2C3E50")
    tbl9[(0,j)].set_text_props(color="white", fontweight="bold")
for i in range(1, len(rows9)+1, 2):
    for j in range(2):
        tbl9[(i,j)].set_facecolor("#F2F3F4")
# Highlight Ts/Tv row
for j in range(2):
    tbl9[(4,j)].set_facecolor("#FADBD8")
    tbl9[(4,j)].set_text_props(fontweight="bold")
ax9.set_title("SNP Analysis Summary", fontweight="bold", fontsize=11, pad=20)

plt.tight_layout(rect=[0,0,1,0.96])
plt.savefig("outputs/snp_analysis_dashboard.png",
            dpi=150, bbox_inches="tight")
plt.close()
print("✅ Dashboard saved → outputs/snp_analysis_dashboard.png")

# ─────────────────────────────────────────────────────────────
# FINAL SUMMARY
# ─────────────────────────────────────────────────────────────

print("\n" + "="*60)
print("FINAL SUMMARY — E. coli K-12 SNP ANALYSIS")
print("="*60)
for k, v in summary.items():
    print(f"  {k:30s}: {v}")

print(f"\n  Mutation spectrum:")
for cls, n in mut_counts.items():
    print(f"    {cls}: {n:,} ({n/len(df)*100:.1f}%)")

print(f"\n  Genotype breakdown:")
for gt, n in gt_counts.items():
    print(f"    {gt:12s}: {n:,} ({n/len(df)*100:.1f}%)")

print(f"\n  Impact prediction:")
for imp, n in impact_counts.items():
    print(f"    {imp:25s}: {n:,} ({n/len(df)*100:.1f}%)")

print("\n✅ All outputs saved!")
