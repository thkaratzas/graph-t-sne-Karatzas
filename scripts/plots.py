# Figures for the report, from results/times-<machine>.csv
#   python3 scripts/plots.py results/times-MacBook-Air-Thanasis.csv
# Writes figures/*.pdf and figures/*.png, and prints the numbers behind them.
import sys
import math
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV = sys.argv[1] if len(sys.argv) > 1 else "results/times-MacBook-Air-Thanasis.csv"
OUT = "figures/"
PCORES = 4            # Apple M1: 4 performance + 4 efficiency cores

# fixed colors per version (same everywhere), text stays black/gray
COLOR = {"omp": "#2a78d6", "cilk": "#eb6834", "pthreads": "#1baf7a"}
NAME = {"omp": "OpenMP", "cilk": "OpenCilk", "pthreads": "Pthreads"}
GRAPHS = ["k49_norm_10NN", "com-Amazon", "coPapersDBLP"]
GCOLOR = {"k49_norm_10NN": "#2a78d6", "com-Amazon": "#eb6834", "coPapersDBLP": "#1baf7a"}
SIZE = {"k49_norm_10NN": (38547, 618158),
        "com-Amazon": (334863, 1851744),
        "coPapersDBLP": (540486, 30491458)}

plt.rcParams.update({
    "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#888888", "axes.labelcolor": "#333333",
    "xtick.color": "#555555", "ytick.color": "#555555",
    "axes.grid": True, "grid.color": "#e6e6e6", "grid.linewidth": 0.6,
    "lines.linewidth": 2, "lines.markersize": 5, "legend.frameon": False,
})

d = pd.read_csv(CSV)
d["bh"] = d["build"] + d["com"] + d["walk"]
phases = ["P", "attr", "build", "com", "walk", "grad", "bh", "direct"]
# minimum of the repetitions for every configuration
t = d.groupby(["graph", "version", "schedule", "threads"])[phases].min()


def seq(graph, phase):
    return t.loc[(graph, "seq", "start", 1), phase]


def curve(graph, version, schedule, phase):
    c = t.loc[(graph, version, schedule)][phase]
    return c.index.values, c.values


def ecore_band(ax, pmax):
    ax.axvspan(PCORES + 0.5, pmax + 0.5, color="#f2f2f2", zorder=0)
    y0, y1 = ax.get_ylim()
    ax.text(PCORES + 0.7, y0 + 0.03 * (y1 - y0), "E-cores", va="bottom",
            color="#777777", fontsize=8)


def save(fig, name):
    fig.tight_layout()
    fig.savefig(OUT + name + ".pdf")
    fig.savefig(OUT + name + ".png", dpi=200)
    plt.close(fig)
    print("wrote", OUT + name + ".pdf")


# ---------- 1. speedup per phase, dynamic scheduling ----------
panels = [("coPapersDBLP", "P", "P (DBLP)"),
          ("coPapersDBLP", "attr", "attraction (DBLP)"),
          ("coPapersDBLP", "walk", "Barnes-Hut walk (DBLP)"),
          ("k49_norm_10NN", "direct", "direct O(n^2) (k49)")]
fig, axes = plt.subplots(1, 4, figsize=(10, 2.8), sharey=True)
for ax, (g, ph, title) in zip(axes, panels):
    pmax = 8
    ax.plot([1, pmax], [1, pmax], color="#bbbbbb", linewidth=1, linestyle="--")
    for v in ["omp", "cilk", "pthreads"]:
        sch = "ws" if v == "cilk" else "dynamic"
        p, y = curve(g, v, sch, ph)
        ax.plot(p, seq(g, ph) / y, marker="o", color=COLOR[v], label=NAME[v])
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("threads")
    ax.set_xticks(range(1, pmax + 1))
    ax.set_ylim(0, pmax + 0.5)
    ecore_band(ax, pmax)
axes[0].set_ylabel("speedup over sequential")
axes[0].legend(loc="upper left")
save(fig, "speedup")

# ---------- 2. static vs dynamic, phase P, all graphs (OpenMP) ----------
fig, axes = plt.subplots(1, 3, figsize=(10, 2.8))
for ax, g in zip(axes, GRAPHS):
    for sch, style in [("dynamic", "-"), ("static", "--")]:
        p, y = curve(g, "omp", sch, "P")
        ax.plot(p, y * 1000, marker="o", linestyle=style, color=COLOR["omp"],
                label=sch + (",256" if sch == "dynamic" else ""))
    ax.set_title("P, " + g, fontsize=9)
    ax.set_xlabel("threads")
    ax.set_xticks(range(1, 9))
    ax.set_ylim(bottom=0)
    ecore_band(ax, 8)
axes[0].set_ylabel("time (ms)")
axes[0].legend(loc="upper right")
save(fig, "schedule")

# ---------- 3. time per phase vs graph, sequential and 8 threads ----------
order = ["P", "attr", "build", "com", "walk", "grad"]
fig, axes = plt.subplots(1, 2, figsize=(10, 2.8), sharey=True)
for ax, (label, ver, sch, p) in zip(axes, [("sequential", "seq", "start", 1),
                                           ("OpenMP, 8 threads", "omp", "dynamic", 8)]):
    for k, g in enumerate(GRAPHS):
        y = [t.loc[(g, ver, sch, p), ph] * 1000 for ph in order]
        x = [i + (k - 1) * 0.22 for i in range(len(order))]
        ax.scatter(x, y, color=GCOLOR[g], s=30, zorder=3, label=g)
    ax.set_yscale("log")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order)
    ax.set_title(label, fontsize=9)
axes[0].set_ylabel("time (ms, log scale)")
axes[0].legend(loc="upper center")
save(fig, "phases")

# ---------- 4. Barnes-Hut vs direct on k49 ----------
g = "k49_norm_10NN"
fig, ax = plt.subplots(figsize=(4.2, 2.8))
for col, label, style in [("direct", "direct O(n^2)", "--"),
                          ("bh", "Barnes-Hut (build+com+walk)", "-")]:
    p, y = curve(g, "omp", "dynamic", col)
    ax.plot(p, y * 1000, marker="o", linestyle=style, color=COLOR["omp"], label=label)
ax.set_yscale("log")
ax.set_xlabel("threads (OpenMP)")
ax.set_ylabel("time (ms, log scale)")
ax.set_xticks(range(1, 9))
ax.set_title("repulsion on k49_norm_10NN", fontsize=9)
ax.legend(loc="upper right")
save(fig, "bh_vs_direct")

# ---------- numbers behind the figures ----------
print("\nspeedup at 4 / 8 threads (dynamic, Cilk: work stealing)")
for g, ph, title in panels:
    row = []
    for v in ["omp", "cilk", "pthreads"]:
        sch = "ws" if v == "cilk" else "dynamic"
        s4 = seq(g, ph) / t.loc[(g, v, sch, 4), ph]
        s8 = seq(g, ph) / t.loc[(g, v, sch, 8), ph]
        row.append("%s %.2f / %.2f" % (NAME[v], s4, s8))
    print("  %-24s %s" % (title, "   ".join(row)))

print("\nstatic / dynamic time ratio, phase P (OpenMP)")
for g in GRAPHS:
    r = [t.loc[(g, "omp", "static", p), "P"] / t.loc[(g, "omp", "dynamic", p), "P"]
         for p in [2, 4, 8]]
    print("  %-15s p=2: %.2f  p=4: %.2f  p=8: %.2f" % (g, r[0], r[1], r[2]))

# crossover: direct ~ a n^2, Barnes-Hut ~ b n log2 n, constants from k49
g = "k49_norm_10NN"
n = SIZE[g][0]
a = seq(g, "direct") / n ** 2
b = seq(g, "bh") / (n * math.log2(n))
lo, hi = 2.0, float(n)
for _ in range(100):
    mid = (lo + hi) / 2
    if a * mid ** 2 < b * mid * math.log2(mid):
        lo = mid
    else:
        hi = mid
print("\nBarnes-Hut vs direct on k49 (sequential): %.1fx faster" %
      (seq(g, "direct") / seq(g, "bh")))
print("model crossover (direct faster below): n ~ %.0f" % lo)

print("\nthrottling check: sequential end / start - 1")
for g in GRAPHS:
    e = t.loc[(g, "seq", "end", 1)]
    s = t.loc[(g, "seq", "start", 1)]
    print("  %-15s P %+.1f%%  walk %+.1f%%" %
          (g, 100 * (e["P"] / s["P"] - 1), 100 * (e["walk"] / s["walk"] - 1)))
