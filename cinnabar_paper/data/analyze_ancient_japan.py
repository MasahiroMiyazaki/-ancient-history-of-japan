#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analyze_ancient_japan.py
────────────────────────────────────────────────────────────────────────────
物部系神社(n=17)・百名山神社（n=18）・一之宮（n=51）と古代鉱山（n=27）の
空間統計分析スクリプト

【分析内容】
  1. Haversine式による大圏距離計算
  2. 記述統計（平均・中央値・標準偏差・30km圏内率など）
  3. Mann-Whitney U検定（ノンパラメトリック）
  4. Cohen's d 効果量
  5. Bonferroni補正による多重検定評価
  6. 検出力分析（Post-hoc power analysis）
  7. 距離分布ヒストグラム、箱ひげ図、CDF、ランキング図、閾値分析図の生成

【対照群の選定方針】
  百名山神社群（n=18）: 深田久弥『日本百名山』（1964年）に選定された山岳のうち
  福島以南の本州に所在し、山岳信仰を示す神社が現存するものを機械的全件抽出。
  箱根神社（箱根山）・戸隠神社（戸隠山）はこの基準を満たさないため含めない。

【入力ファイル】
  data/mononobe_shrines_honshu.csv   : 物部系神社 n=17
  data/hyakumeizan_shrines.csv       : 百名山神社  n=18
  data/ichinomiya_shrines.csv        : 一之宮      n=51
  data/mines_ancient_honshu.csv      : 古代鉱山    n=27

【出力】
  images/0_histograms_distribution.png  : 距離分布ヒストグラム（三群）
  images/1_threshold_analysis.png       : 距離閾値別カバレッジ
  images/2_boxplot_detailed.png         : 詳細箱ひげ図（三群）
  images/3_mononobe_ranking.png         : 物部系神社17社 距離ランキング
  images/4_cdf_analysis.png             : 累積分布関数（CDF）

────────────────────────────────────────────────────────────────────────────
"""

import csv
import math
import os
from typing import List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import matplotlib.patches as mpatches
import numpy as np
from scipy import stats

# ── フォント設定（日本語対応） ────────────────────────────────────────────
import platform
import os
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

def set_universal_font():
    """OSごとの日本語フォントを自動探索して適用する"""
    system = platform.system()
    
    # OSごとの代表的な日本語フォントパスリスト
    font_paths = []
    if system == "Windows":
        font_paths = [
            "C:\\Windows\\Fonts\\msgothic.ttc",  # MS Gothic
            "C:\\Windows\\Fonts\\meiryo.ttc"     # Meiryo
        ]
    elif system == "Darwin":  # macOS
        font_paths = [
            "/System/Library/Fonts/ヒラギノ角ゴシック W3.ttc",
            "/Library/Fonts/Arial Unicode.ttf"
        ]
    elif system == "Linux":   # Linux (Ubuntu, Debian, etc.)
        font_paths = [
            "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc", # Debian/Ubuntu
            "/usr/share/fonts/truetype/noto/NotoSansCJK-Regular.ttc",
            "/usr/share/fonts/opentype/ipafont-gothic/ipag.ttf"       # IPAフォント
        ]

    # フォントの適用
    found = False
    for path in font_paths:
        if os.path.exists(path):
            # 外部ファイルとして読み込む
            fm.fontManager.addfont(path)
            prop = fm.FontProperties(fname=path)
            plt.rcParams["font.family"] = prop.get_name()
            print(f"[Font Info] Set font to: {prop.get_name()} ({path})")
            found = True
            break
    
    if not found:
        print("[Warning] No Japanese font found. Falling back to default.")
        plt.rcParams["font.family"] = "sans-serif"

    # マイナス記号の文字化け対策
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["figure.dpi"] = 150

set_universal_font()

# ── 定数 ─────────────────────────────────────────────────────────────────
EARTH_R = 6371.0          # km
THRESHOLD_KM = 30         # 実効支配距離の閾値
COLOR_M = "#c0392b"       # 物部系神社（赤）
COLOR_H = "#2980b9"       # 百名山神社（青）
COLOR_I = "#27ae60"       # 一之宮（緑）
ALPHA = 0.01              # 有意水準
ALPHA_FWER = 0.05         # 族全体有意水準（Bonferroni補正の出発点）
BONFERRONI_N = 2          # 検定回数
ALPHA_BONF = ALPHA_FWER / BONFERRONI_N  # Bonferroni補正後 α' = 0.025


# ════════════════════════════════════════════════════════════════════════════
#  1. ユーティリティ関数
# ════════════════════════════════════════════════════════════════════════════

def haversine(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """Haversine式による大圏距離（km）"""
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = math.sin(dlat / 2) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
    return 2 * math.asin(math.sqrt(a)) * EARTH_R


def nearest_mine_distance(shrine_lat: float, shrine_lon: float,
                           mines: List[Tuple]) -> Tuple[float, str]:
    """最近接鉱山までの距離と鉱山名を返す"""
    min_d = float("inf")
    min_name = ""
    for m in mines:
        d = haversine(shrine_lat, shrine_lon, m["lat"], m["lon"])
        if d < min_d:
            min_d = d
            min_name = m["name"]
    return min_d, min_name


import pandas as pd # ファイルの先頭（import群の中）にこれが必要です

def load_csv(path: str) -> List[dict]:
    # sep=None と engine='python' にすると、タブでもカンマでも自動判定してくれます
    df = pd.read_csv(path, sep=None, engine='python')
    return df.to_dict(orient='records')


def cohen_d(x: np.ndarray, y: np.ndarray) -> float:
    """Pooled SD による Cohen's d"""
    nx, ny = len(x), len(y)
    pooled_std = math.sqrt(((nx - 1) * np.std(x, ddof=1) ** 2 +
                             (ny - 1) * np.std(y, ddof=1) ** 2) / (nx + ny - 2))
    return abs(np.mean(x) - np.mean(y)) / pooled_std if pooled_std > 0 else 0.0


def effect_label(d: float) -> str:
    if d >= 0.8:
        return "Large"
    elif d >= 0.5:
        return "Medium"
    elif d >= 0.2:
        return "Small"
    return "Negligible"


def post_hoc_power(d: float, n1: int, n2: int, alpha: float = 0.01) -> float:
    """
    Two-sample t-test の事後検出力（Mann-Whitney の近似として利用）
    Cohen's d と両サンプルサイズから計算。
    """
    se = math.sqrt(1 / n1 + 1 / n2)
    ncp = d / se  # non-centrality parameter
    z_alpha = stats.norm.ppf(1 - alpha / 2)
    return float(stats.norm.cdf(ncp - z_alpha) + stats.norm.cdf(-ncp - z_alpha))


# ════════════════════════════════════════════════════════════════════════════
#  2. データ読み込みと距離計算
# ════════════════════════════════════════════════════════════════════════════

BASE = os.path.dirname(os.path.abspath(__file__))
DATA = BASE  # スクリプトと同じ場所にあると指定する
IMG  = os.path.join(BASE, "images")
os.makedirs(IMG, exist_ok=True)

print("=" * 60)
print("  古代日本 空間統計分析")
print("  物部神社・百名山神社・一之宮 vs 古代鉱山")
print("=" * 60)

# --- 読み込み
raw_mononobe   = load_csv(os.path.join(DATA, "mononobe_shrines_honshu.csv"))
raw_hyakumeizan = load_csv(os.path.join(DATA, "hyakumeizan_shrines.csv"))
raw_ichinomiya = load_csv(os.path.join(DATA, "ichinomiya_shrines_full.csv"))
raw_mines      = load_csv(os.path.join(DATA, "mines_ancient_honshu.csv"))


mines = []
for m in raw_mines:
    # 辞書のキー名を、先ほどCSVの1行目に付けた名前に合わせます
    mines.append({
        "name": m["name"], 
        "lat": float(m["lat"]), 
        "lon": float(m["lon"])
    })

# --- 距離計算
def calc_distances(shrines, mines, label):
    results = []
    for s in shrines:
        d, mine_name = nearest_mine_distance(float(s["lat"]), float(s["lon"]), mines)
        results.append({
            "shrine": s["name"],
            "lat": float(s["lat"]),
            "lon": float(s["lon"]),
            "dist_km": round(d, 2),
            "nearest_mine": mine_name,
        })
    dists = np.array([r["dist_km"] for r in results])
    return results, dists

mononobe_data,    d_mono  = calc_distances(raw_mononobe,    mines, "物部系神社")
hyakumeizan_data, d_hyaku = calc_distances(raw_hyakumeizan, mines, "百名山神社")
ichinomiya_data,  d_ichi  = calc_distances(raw_ichinomiya,  mines, "一之宮")

print(f"\nデータ確認:")
print(f"  物部系神社  : {len(d_mono)}社")
print(f"  百名山神社  : {len(d_hyaku)}社")
print(f"  一之宮      : {len(d_ichi)}社")
print(f"  古代鉱山    : {len(mines)}鉱山")


# ════════════════════════════════════════════════════════════════════════════
#  3. 記述統計
# ════════════════════════════════════════════════════════════════════════════

def describe(arr: np.ndarray, label: str) -> dict:
    n = len(arr)
    w30 = (arr <= THRESHOLD_KM).sum()
    w10 = (arr <= 10).sum()
    return {
        "label": label,
        "n": n,
        "mean": round(arr.mean(), 2),
        "median": round(float(np.median(arr)), 2),
        "std": round(arr.std(ddof=1), 2),
        "min": round(arr.min(), 2),
        "max": round(arr.max(), 2),
        "within10_n": int(w10),
        "within10_pct": round(w10 / n * 100, 1),
        "within30_n": int(w30),
        "within30_pct": round(w30 / n * 100, 1),
    }

stats_mono  = describe(d_mono,  "物部系神社 (n=17)")
stats_hyaku = describe(d_hyaku, "百名山神社 (n=18)")
stats_ichi  = describe(d_ichi,  "一之宮     (n=51)")

print("\n" + "─" * 60)
print("【記述統計】")
print("─" * 60)
header = f"{'統計量':<16} {'物部系':>10} {'百名山':>10} {'一之宮':>10}"
print(header)
print("─" * 50)
rows = [
    ("n",         stats_mono["n"],        stats_hyaku["n"],        stats_ichi["n"]),
    ("平均 (km)",  stats_mono["mean"],     stats_hyaku["mean"],     stats_ichi["mean"]),
    ("中央値 (km)",stats_mono["median"],   stats_hyaku["median"],   stats_ichi["median"]),
    ("標準偏差",   stats_mono["std"],      stats_hyaku["std"],      stats_ichi["std"]),
    ("最小値",     stats_mono["min"],      stats_hyaku["min"],      stats_ichi["min"]),
    ("最大値",     stats_mono["max"],      stats_hyaku["max"],      stats_ichi["max"]),
    ("10km以内数", stats_mono["within10_n"], stats_hyaku["within10_n"], stats_ichi["within10_n"]),
    ("10km以内%",  stats_mono["within10_pct"], stats_hyaku["within10_pct"], stats_ichi["within10_pct"]),
    ("30km以内数", stats_mono["within30_n"], stats_hyaku["within30_n"], stats_ichi["within30_n"]),
    ("30km以内%",  stats_mono["within30_pct"], stats_hyaku["within30_pct"], stats_ichi["within30_pct"]),
]
for name, v1, v2, v3 in rows:
    print(f"  {name:<14} {str(v1):>10} {str(v2):>10} {str(v3):>10}")


# ════════════════════════════════════════════════════════════════════════════
#  4. 推測統計 — Mann-Whitney U 検定 + Cohen's d
# ════════════════════════════════════════════════════════════════════════════

print("\n" + "─" * 60)
print("【推測統計】")
print("─" * 60)

# 分析1: 物部 vs 百名山神社
u1, p1 = stats.mannwhitneyu(d_mono, d_hyaku, alternative="two-sided")
d1 = cohen_d(d_mono, d_hyaku)
pw1 = post_hoc_power(d1, len(d_mono), len(d_hyaku), alpha=ALPHA)

print(f"\n  ■ 物部系神社 (n={len(d_mono)}) vs 百名山神社 (n={len(d_hyaku)})")
print(f"    Mann-Whitney U統計量 : {u1:.1f}")
print(f"    p値 (両側検定): {p1:.6f}")
print(f"    Cohen's d            : {d1:.4f}  [{effect_label(d1)}]")
print(f"    有意水準 α={ALPHA}    : {'有意 ✓' if p1 < ALPHA else '非有意'}")
print(f"    Bonferroni α'={ALPHA_BONF}: {'有意 ✓' if p1 < ALPHA_BONF else '非有意'}")
print(f"    事後検出力            : {pw1:.3f} ({pw1*100:.1f}%)")

# 分析2: 物部 vs 一之宮
u2, p2 = stats.mannwhitneyu(d_mono, d_ichi, alternative="two-sided")
d2 = cohen_d(d_mono, d_ichi)
pw2 = post_hoc_power(d2, len(d_mono), len(d_ichi), alpha=ALPHA)

print(f"\n  ■ 物部系神社 (n={len(d_mono)}) vs 一之宮 (n={len(d_ichi)})")
print(f"    Mann-Whitney U統計量 : {u2:.1f}")
print(f"    p値 (両側検定) : {p2:.6f}")
print(f"    Cohen's d            : {d2:.4f}  [{effect_label(d2)}]")
print(f"    有意水準 α={ALPHA}    : {'有意 ✓' if p2 < ALPHA else '非有意'}")
print(f"    Bonferroni α'={ALPHA_BONF}: {'有意 ✓' if p2 < ALPHA_BONF else '非有意'}")
print(f"    事後検出力            : {pw2:.3f} ({pw2*100:.1f}%)")

# 感度分析: n=14 (物部社名のみ) vs 百名山神社 / 一之宮
exclude_extended = ["大津神社（飛騨市）", "彌彦神社（弥彦村）", "石上神宮"]
d_mono14 = np.array([
    mononobe_data[i]["dist_km"] for i, s in enumerate(raw_mononobe)
    if s["name"] not in exclude_extended
])

u1b, p1b = stats.mannwhitneyu(d_mono14, d_hyaku, alternative="two-sided")
d1b = cohen_d(d_mono14, d_hyaku)
u2b, p2b = stats.mannwhitneyu(d_mono14, d_ichi, alternative="two-sided")
d2b = cohen_d(d_mono14, d_ichi)

print(f"\n  ■ 感度分析: 物部社名のみ (n={len(d_mono14)}) vs 百名山神社 (n={len(d_hyaku)})")
print(f"    p値: {p1b:.6f}  Cohen's d: {d1b:.4f}  [{effect_label(d1b)}]")
print(f"    有意 α={ALPHA}: {'✓' if p1b < ALPHA else '✗ (α=0.05: ' + ('✓' if p1b < 0.05 else '✗') + ')'}")

print(f"\n  ■ 感度分析: 物部社名のみ (n={len(d_mono14)}) vs 一之宮 (n={len(d_ichi)})")
print(f"    p値: {p2b:.6f}  Cohen's d: {d2b:.4f}  [{effect_label(d2b)}]")
print(f"    有意 α={ALPHA}: {'✓' if p2b < ALPHA else '✗'}")

# 正規性検定
_, p_shapiro_mono  = stats.shapiro(d_mono)
_, p_shapiro_hyaku = stats.shapiro(d_hyaku)
_, p_shapiro_ichi  = stats.shapiro(d_ichi)
print(f"\n  ■ Shapiro-Wilk正規性検定")
print(f"    物部系: p={p_shapiro_mono:.4f}  {'正規' if p_shapiro_mono > 0.05 else '非正規'}")
print(f"    百名山: p={p_shapiro_hyaku:.4f}  {'正規' if p_shapiro_hyaku > 0.05 else '非正規'}")
print(f"    一之宮: p={p_shapiro_ichi:.4f}  {'正規' if p_shapiro_ichi > 0.05 else '非正規'}")

# 等分散性検定
lev_stat, lev_p = stats.levene(d_mono, d_hyaku)
print(f"\n  ■ Levene等分散性検定 (物部 vs 百名山)")
print(f"    F={lev_stat:.3f}, p={lev_p:.4f}  {'等分散' if lev_p > 0.05 else '非等分散'}")


# ════════════════════════════════════════════════════════════════════════════
#  5. グラフ描画
# ════════════════════════════════════════════════════════════════════════════

LABEL_M = f"物部系神社 (n={len(d_mono)})"
LABEL_H = f"百名山神社 (n={len(d_hyaku)})"
LABEL_I = f"一之宮 (n={len(d_ichi)})"


# ─── 図0：距離分布ヒストグラム（三群） ─────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=False)
fig.suptitle("三群の最近接鉱山距離分布", fontsize=14, fontweight="bold", y=1.01)

bins = np.arange(0, 145, 10)
configs = [
    (d_mono,  COLOR_M, LABEL_M),
    (d_hyaku, COLOR_H, LABEL_H),
    (d_ichi,  COLOR_I, LABEL_I),
]
for ax, (data, color, label) in zip(axes, configs):
    ax.hist(data, bins=bins, color=color, alpha=0.75, edgecolor="white", linewidth=0.8)
    med = np.median(data)
    ax.axvline(med, color="#333", linestyle="--", linewidth=1.5,
               label=f"中央値: {med:.1f} km")
    ax.axvline(THRESHOLD_KM, color="gray", linestyle=":", linewidth=1.2, alpha=0.8)
    ax.set_title(label, fontsize=11, pad=8)
    ax.set_xlabel("最近接鉱山距離 (km)", fontsize=10)
    ax.set_ylabel("神社数", fontsize=10)
    ax.legend(fontsize=9)
    ax.set_xlim(0, 140)
    ax.grid(axis="y", alpha=0.3)

plt.tight_layout()
path0 = os.path.join(IMG, "0_histograms_distribution.png")
plt.savefig(path0, dpi=150, bbox_inches="tight")
plt.close()
print(f"\n  → 保存: {path0}")


# ─── 図1：距離閾値別カバレッジ ───────────────────────────────────────────
thresholds = [10, 20, 30, 50, 100]
cov_mono  = [(d_mono  <= t).mean() * 100 for t in thresholds]
cov_hyaku = [(d_hyaku <= t).mean() * 100 for t in thresholds]
cov_ichi  = [(d_ichi  <= t).mean() * 100 for t in thresholds]

fig, ax = plt.subplots(figsize=(9, 6))
x = np.arange(len(thresholds))
w = 0.26
rects1 = ax.bar(x - w, cov_mono,  w, label=LABEL_M, color=COLOR_M, alpha=0.85)
rects2 = ax.bar(x,     cov_hyaku, w, label=LABEL_H, color=COLOR_H, alpha=0.85)
rects3 = ax.bar(x + w, cov_ichi,  w, label=LABEL_I, color=COLOR_I, alpha=0.85)

for rects in [rects1, rects2, rects3]:
    for rect in rects:
        h = rect.get_height()
        ax.annotate(f"{h:.0f}%",
                    xy=(rect.get_x() + rect.get_width() / 2, h),
                    xytext=(0, 3), textcoords="offset points",
                    ha="center", va="bottom", fontsize=8)

ax.set_xticks(x)
ax.set_xticklabels([f"{t} km" for t in thresholds])
ax.set_xlabel("距離閾値", fontsize=11)
ax.set_ylabel("閾値以内の神社の割合 (%)", fontsize=11)
ax.set_title("距離閾値別・鉱山圏内神社の割合（三群比較）", fontsize=13, fontweight="bold")
ax.legend(fontsize=10)
ax.set_ylim(0, 110)
ax.grid(axis="y", alpha=0.3)
fig.tight_layout()
path1 = os.path.join(IMG, "1_threshold_analysis.png")
plt.savefig(path1, dpi=150, bbox_inches="tight")
plt.close()
print(f"  → 保存: {path1}")


# ─── 図2：詳細箱ひげ図 ───────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 7))

bp = ax.boxplot(
    [d_mono, d_hyaku, d_ichi],
    patch_artist=True,
    notch=False,
    widths=0.5,
    medianprops=dict(color="#1a3a6b", linewidth=2.5),
    whiskerprops=dict(linewidth=1.5),
    capprops=dict(linewidth=1.5),
    flierprops=dict(marker="o", markersize=5, alpha=0.6),
)
colors = [COLOR_M, COLOR_H, COLOR_I]
for patch, color in zip(bp["boxes"], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.65)

# 平均値プロット
means = [d_mono.mean(), d_hyaku.mean(), d_ichi.mean()]
ax.scatter([1, 2, 3], means, color="red", zorder=5, s=60, label="平均値", marker="D")

ax.axhline(THRESHOLD_KM, color="gray", linestyle="--", linewidth=1.3, alpha=0.7,
           label=f"{THRESHOLD_KM} km 閾値")
ax.set_xticks([1, 2, 3])
ax.set_xticklabels([LABEL_M, LABEL_H, LABEL_I], fontsize=10)
ax.set_ylabel("最近接鉱山距離 (km)", fontsize=11)
ax.set_title("三群の鉱山距離分布（箱ひげ図）", fontsize=13, fontweight="bold")
ax.legend(fontsize=10)
ax.grid(axis="y", alpha=0.3)

# p値注記
y_max = max(d_mono.max(), d_hyaku.max(), d_ichi.max()) * 1.05
ax.annotate(
    f"p={p1:.4f} (Mann-Whitney)\nd={d1:.3f} [{effect_label(d1)}]",
    xy=(1.5, y_max * 0.88), ha="center", fontsize=9,
    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor=COLOR_M, alpha=0.8),
)

fig.tight_layout()
path2 = os.path.join(IMG, "2_boxplot_detailed.png")
plt.savefig(path2, dpi=150, bbox_inches="tight")
plt.close()
print(f"  → 保存: {path2}")


# ─── 図3：物部系神社17社 ランキング ─────────────────────────────────────
sorted_mono = sorted(mononobe_data, key=lambda x: x["dist_km"])
names  = [s["shrine"] for s in sorted_mono]
dists  = [s["dist_km"] for s in sorted_mono]
mines_ = [s["nearest_mine"] for s in sorted_mono]
colors_bar = [COLOR_M if d <= THRESHOLD_KM else "#e8a0a0" for d in dists]

fig, ax = plt.subplots(figsize=(13, 7))
bars = ax.barh(range(len(names)), dists, color=colors_bar, alpha=0.85, edgecolor="white")
ax.axvline(THRESHOLD_KM, color="#333", linestyle="--", linewidth=1.8,
           label=f"{THRESHOLD_KM} km 閾値")

for i, (d, mine) in enumerate(zip(dists, mines_)):
    ax.text(d + 0.8, i, f"{d:.1f} km  [{mine}]", va="center", fontsize=8.5, color="#222")

ax.set_yticks(range(len(names)))
ax.set_yticklabels(names, fontsize=9)
ax.set_xlabel("最近接鉱山距離 (km)", fontsize=11)
ax.set_title("物部系神社17社 — 最近接鉱山距離（昇順）", fontsize=13, fontweight="bold")

within30 = sum(1 for d in dists if d <= THRESHOLD_KM)
pct30 = within30 / len(dists) * 100
patch1 = mpatches.Patch(color=COLOR_M, alpha=0.85, label=f"30km以内 ({within30}社, {pct30:.1f}%)")
patch2 = mpatches.Patch(color="#e8a0a0", alpha=0.85, label=f"30km超 ({len(dists)-within30}社)")
ax.legend(handles=[patch1, patch2], fontsize=10, loc="upper right")
ax.set_xlim(0, max(dists) * 1.35)
ax.grid(axis="x", alpha=0.3)
ax.invert_yaxis()
fig.tight_layout()
path3 = os.path.join(IMG, "3_mononobe_ranking.png")
plt.savefig(path3, dpi=150, bbox_inches="tight")
plt.close()
print(f"  → 保存: {path3}")


# ─── 図4：累積分布関数（CDF） ────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 7))

for data, color, label in [(d_mono, COLOR_M, LABEL_M),
                            (d_hyaku, COLOR_H, LABEL_H),
                            (d_ichi, COLOR_I, LABEL_I)]:
    xs = np.sort(data)
    ys = np.arange(1, len(xs) + 1) / len(xs) * 100
    ax.step(xs, ys, where="post", color=color, linewidth=2.2, label=label)

ax.axvline(THRESHOLD_KM, color="#555", linestyle="--", linewidth=1.5, alpha=0.8,
           label=f"{THRESHOLD_KM} km 閾値")

# 30km地点の読み取り値注記
for data, color in [(d_mono, COLOR_M), (d_hyaku, COLOR_H), (d_ichi, COLOR_I)]:
    y30 = (data <= THRESHOLD_KM).mean() * 100
    ax.scatter([THRESHOLD_KM], [y30], color=color, zorder=5, s=55)
    ax.annotate(f"{y30:.1f}%", xy=(THRESHOLD_KM, y30),
                xytext=(2.5, 2.5), textcoords="offset points",
                color=color, fontsize=9, fontweight="bold")

ax.set_xlabel("最近接鉱山距離 (km)", fontsize=11)
ax.set_ylabel("神社の累積割合 (%)", fontsize=11)
ax.set_title("三群の鉱山距離 累積分布関数（CDF）", fontsize=13, fontweight="bold")
ax.legend(fontsize=10)
ax.set_xlim(0, 140)
ax.set_ylim(0, 105)
ax.grid(alpha=0.3)
fig.tight_layout()
path4 = os.path.join(IMG, "4_cdf_analysis.png")
plt.savefig(path4, dpi=150, bbox_inches="tight")
plt.close()
print(f"  → 保存: {path4}")


# ════════════════════════════════════════════════════════════════════════════
#  6. 最終サマリー
# ════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 60)
print("【最終統計サマリー】")
print("=" * 60)
print(f"\n  物部系神社 (n=17): 平均 {d_mono.mean():.2f}km / 中央値 {np.median(d_mono):.2f}km")
print(f"\n  対 百名山神社 (n={len(d_hyaku)}):")
print(f"    Mann-Whitney U検定 p値: {p1:.6f}")
print(f"    効果量 (Cohen's d)    : {d1:.4f}  [{effect_label(d1)}]")
print(f"    U統計量               : {u1:.1f}")
print(f"\n  対 一之宮 (n={len(d_ichi)}):")
print(f"    Mann-Whitney U検定 p値: {p2:.6f}")
print(f"    効果量 (Cohen's d)    : {d2:.4f}  [{effect_label(d2)}]")
print(f"    U統計量               : {u2:.1f}")
print(f"\n  30km圏内到達率:")
print(f"    物部系神社: {stats_mono['within30_pct']}% ({stats_mono['within30_n']}/{stats_mono['n']})")
print(f"    百名山神社: {stats_hyaku['within30_pct']}% ({stats_hyaku['within30_n']}/{stats_hyaku['n']})")
print(f"    一之宮    : {stats_ichi['within30_pct']}% ({stats_ichi['within30_n']}/{stats_ichi['n']})")
print()
print("  生成画像:")
for p in [path0, path1, path2, path3, path4]:
    print(f"    {os.path.basename(p)}")
print("\n  完了。")