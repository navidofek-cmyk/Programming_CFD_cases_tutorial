"""Generate CFD domain figure for the laminar channel SIMPLE case."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
import numpy as np

fig, axes = plt.subplots(1, 2, figsize=(13, 4.5),
                         gridspec_kw={"width_ratios": [3, 1.6]})

# ── LEFT: domain overview ─────────────────────────────────────────────────────
ax = axes[0]
ax.set_xlim(-0.4, 7.0)
ax.set_ylim(-0.55, 1.65)
ax.set_aspect("equal")
ax.axis("off")

Lx, Ly = 6.0, 1.0

# Channel walls
wall_kw = dict(color="k", linewidth=2)
ax.plot([0, Lx], [0,  0 ], **wall_kw)   # bottom wall
ax.plot([0, Lx], [Ly, Ly], **wall_kw)   # top wall
ax.plot([0,  0], [0,  Ly], **wall_kw)   # inlet
ax.plot([Lx, Lx], [0, Ly], linewidth=1.5, color="gray", linestyle="--")  # outlet

# Dimension annotations
ax.annotate("", xy=(Lx, -0.35), xytext=(0, -0.35),
            arrowprops=dict(arrowstyle="<->", color="k", lw=1.2))
ax.text(Lx/2, -0.47, r"$L_x = 6$", ha="center", va="top", fontsize=11)

ax.annotate("", xy=(-0.28, Ly), xytext=(-0.28, 0),
            arrowprops=dict(arrowstyle="<->", color="k", lw=1.2))
ax.text(-0.38, Ly/2, r"$L_y = 1$", ha="right", va="center", fontsize=11, rotation=90)

# Inlet arrows + profile
y_vals = np.linspace(0.05, 0.95, 9)
U_in = 0.25
for y in y_vals:
    ax.annotate("", xy=(U_in, y), xytext=(0, y),
                arrowprops=dict(arrowstyle="->", color="steelblue", lw=1.4))
ax.text(-0.05, Ly/2, r"$u = U_{in}$" + "\n" + r"$v = 0$",
        ha="right", va="center", fontsize=10, color="steelblue")

# Outlet dashed arrows + profile
y_fine = np.linspace(0, 1, 60)
u_pois = 1.5 * y_fine * (1 - y_fine) * 4
scale = 0.5
for y, u in zip(y_fine[5:-5:6], u_pois[5:-5:6]):
    ax.annotate("", xy=(Lx + u*scale, y), xytext=(Lx, y),
                arrowprops=dict(arrowstyle="->", color="tomato", lw=1.2,
                                linestyle="dashed"))
ax.plot(Lx + u_pois*scale, y_fine, color="tomato", lw=2)
ax.text(Lx + 0.05, Ly/2, r"$\partial u/\partial x = 0$" + "\n" + r"$p = 0$",
        ha="left", va="center", fontsize=10, color="tomato")

# Wall labels
ax.text(Lx/2, Ly + 0.1, r"No-slip wall  ($u = 0$, ghost cell)",
        ha="center", va="bottom", fontsize=10, color="k")
ax.text(Lx/2, -0.12, r"No-slip wall  ($u = 0$, ghost cell)",
        ha="center", va="top", fontsize=10, color="k")

# Coordinate
ax.annotate("", xy=(0.4, 0.12), xytext=(0.0, 0.12),
            arrowprops=dict(arrowstyle="->", color="k", lw=1))
ax.annotate("", xy=(0.0, 0.52), xytext=(0.0, 0.12),
            arrowprops=dict(arrowstyle="->", color="k", lw=1))
ax.text(0.42, 0.10, r"$x$", fontsize=10)
ax.text(0.02, 0.54, r"$y$", fontsize=10)

# Shading for domain interior
ax.fill_between([0, Lx], [0, 0], [Ly, Ly], alpha=0.06, color="steelblue")

ax.set_title("Channel domain  ($\\mathrm{Re}=100$, $N_x\\times N_y = 120\\times 40$)",
             fontsize=12, pad=8)

# ── RIGHT: staggered grid cell ────────────────────────────────────────────────
ax2 = axes[1]
ax2.set_xlim(-0.1, 3.1)
ax2.set_ylim(-0.1, 2.1)
ax2.set_aspect("equal")
ax2.axis("off")
ax2.set_title("Staggered (MAC) cell", fontsize=12, pad=8)

dx, dy = 1.0, 1.0
cells_x, cells_y = 3, 2

# Grid lines
for i in range(cells_x + 1):
    ax2.plot([i*dx, i*dx], [0, cells_y*dy], color="lightgray", lw=0.8)
for j in range(cells_y + 1):
    ax2.plot([0, cells_x*dx], [j*dy, j*dy], color="lightgray", lw=0.8)

# Central cell highlight
ci, cj = 1, 1
rect = mpatches.Rectangle((ci*dx, cj*dy), dx, dy,
                           linewidth=1.5, edgecolor="k", facecolor="lightyellow", zorder=2)
ax2.add_patch(rect)

# p at centre
px, py = (ci + 0.5)*dx, (cj + 0.5)*dy
ax2.plot(px, py, "s", color="navy", ms=9, zorder=3)
ax2.text(px + 0.08, py + 0.07, r"$p_{i,j}$", fontsize=11, color="navy")

# u faces (east/west of highlighted cell)
for (ix, label) in [(ci, r"$u_{i,j}$"), (ci+1, r"$u_{i+1,j}$")]:
    ax2.plot(ix*dx, (cj+0.5)*dy, ">", color="steelblue", ms=10, zorder=3)
    ax2.text(ix*dx + 0.05, (cj+0.5)*dy + 0.1, label, fontsize=10, color="steelblue")

# v faces (north/south)
for (jv, label, off) in [(cj, r"$v_{i,j}$", -0.18), (cj+1, r"$v_{i,j+1}$", 0.06)]:
    ax2.plot((ci+0.5)*dx, jv*dy, "^", color="tomato", ms=10, zorder=3)
    ax2.text((ci+0.5)*dx + 0.05, jv*dy + off, label, fontsize=10, color="tomato")

# Legend
handles = [
    mpatches.Patch(color="lightyellow", ec="k", label="pressure cell"),
    plt.Line2D([0],[0], marker="s", color="navy", ms=8, ls="", label=r"$p$ (cell centre)"),
    plt.Line2D([0],[0], marker=">", color="steelblue", ms=8, ls="", label=r"$u$ (x-face)"),
    plt.Line2D([0],[0], marker="^", color="tomato", ms=8, ls="", label=r"$v$ (y-face)"),
]
ax2.legend(handles=handles, loc="lower left", fontsize=9, framealpha=0.9)

plt.tight_layout()
plt.savefig("domain.png", dpi=150, bbox_inches="tight")
print("Saved domain.png")
