#!/usr/bin/env python3
"""
Study-design schematic for 2026WR044716 (revised).

Five columns: fracture intensity, cation supply, mineralogy, injection strategy,
and aperture evolution. Revision R2: results from the corrected coupled runs\n(gravity omitted, uniform outflow pressure).

Usage
-----
    python3 src/fig_study_design.py --out figures
"""
from __future__ import annotations
import argparse, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon, PathPatch
from matplotlib.path import Path
from matplotlib.colors import LinearSegmentedColormap, to_rgb
import numpy as np

FIG_W = 5.40      # inches; single-column width

# ===================================================================
# NUMBERS -- keep in step with finalised.csv and aperture analysis
# ===================================================================
TITLE = ["Fracture intensity", "Cation supply",
         "Mineralogy", "Injection strategy",
         "Aperture evolution"]

CARD = [
    "$P_{32}$ $\\times$ 0.75$-$2.00\nfive levels, ten\nnetworks each",
    "where the cations\nthat form carbonate\nare released",
    "four carbonates\ncomposition, kinetics\nrate sweeping",
    "continuous injection\nvs shut-in after\n10 years",
    "fixed porosity\nvs coupled\nporosity-permeability",
]
NRUN = ["50 networks", "50 networks", "16 variants",
        "50 matched pairs", "50 matched pairs"]
FIND = [
    "Carbonate follows\npore volume;\nCV 43$-$69 %",
    "Cations released\nacross the network;\n1.25 % captured",
    "Calcite 84$-$87 %;\nno dawsonite\n(sodium-limited)",
    "Shut-in raises\nthe carbonate\n(median 6.2$\\times$)",
    "Evolving-aperture\nmodel primary; fixed\nporosity overfills\ncells in 46 of\n50 networks",
]

BAND_TITLE = ("Carbonate forms where Ca/Mg complexes, carbonate species and a higher local pH coincide:\n"
              "3$-$6 % of the fracture volume contains >90 % of the new carbonate")
BULLETS = [
    "Fracture intensity increases mineralization amount but does not increase "
    "mineralization efficiency",
    "Forsterite dissolves first (95$-$96 % in 2 years), "
    "then diopside and anorthite (>95 % by 50 years)",
    "Magnesite forms first, calcite dominates by 50 years (84$-$87 %). "
    "Surface area has the largest effect",
    "Continuous injection fixes 0.004$-$0.006 % of the carbon. "
    "After shut-in, 97$-$100 % of the dissolved carbon becomes carbonate",
    "",
]

TIER_LABELS = ["Control\nvariable", "Parameter\nswept",
               "Principal\nresult", "Conclusion"]

# ---- palette ------------------------------------------------------
ORANGE = ("#c9781f", "#fbf0e2", "#c9781f")
GREEN  = ("#2e8b57", "#e8f5ee", "#2e8b57")
TEAL   = ("#1a7f8e", "#e4f2f4", "#1a7f8e")
PURPLE = ("#7b4fa3", "#f3ecf9", "#7b4fa3")
RED    = ("#c44e52", "#fce8e9", "#c44e52")
GREY_F, GREY_E, INK = "#eef1f4", "#8794a2", "#333333"
COLS = [ORANGE, GREEN, TEAL, PURPLE, RED]

# ---- typography ---------------------------------------------------
# One knob. Base sizes below are the original 4.6-6.0 pt set; FONT_SCALE
# multiplies all of them, and the box widths and row heights follow.
# Single-column (FIG_W 5.40) tops out near 1.12x because the boxes are then
# at the full usable width; above that, widen FIG_W to double column.
FONT_SCALE = 1.12

FS_LBL, FS_BAN, FS_CARD, FS_RUN, FS_RES = [
    s * FONT_SCALE for s in (5.6, 5.7, 5.0, 4.8, 5.0)]
FS_BAND_TITLE, FS_BULLET = [s * FONT_SCALE for s in (6.0, 5.7)]

# ---- geometry, in points on the final canvas ----------------------
X_LEFT, X_RIGHT = 0.095, 0.980
LBL_X = X_LEFT / 2.0     # centred in the label column, not a loose guess
NCOL = len(TITLE)

BOX_PAD_X = 4.7 * FONT_SCALE     # padding each side of the longest text line
GUTTER_PT = 9.0                  # white space between neighbouring boxes

M_TOP, M_BOT = 4.0, 4.0
# row heights scale with the text they hold; gaps and margins do not
H_BAN = 10.0 * FONT_SCALE
H_CARD = 33.0 * FONT_SCALE
H_RES = 38.0 * FONT_SCALE
G_BC, G_CR, G_RB = [g * FONT_SCALE for g in (14.0, 14.0, 17.0)]
SHOW_ARROWS = True
# One spreading funnel sits in each gap between rows: concave sides sweep in
# from the outer box centres to a shaft, then a block head points at the row
# below. Fill fades from white at the top edge to full colour at the tip.
ARROW_MS = 7.0 * FONT_SCALE   # head length is 0.4 x this, in points
ARROW_DASH = (0, (2.0, 1.5))  # dash pattern for the row connectors
ARROW_INSET = 0.6        # points clear of the box edge at each end
ARROW_SPAN = 0.62        # top width, fraction of the outer-centre distance
ARROW_SHAFT = 11.0 * FONT_SCALE     # shaft width, points
ARROW_HEAD = 24.0 * FONT_SCALE      # head width, points
ARROW_HEAD_F = 0.38                 # head as a fraction of total height
ARROW_COLOR = GREY_E                # colour at the tip
ARROW_TOP = 0.04                    # colour strength at the top edge
# card internals: text hangs from the top edge, run count sits on the bottom
CARD_PAD_T, CARD_PAD_B = 4.0, 3.5
BAND_PAD, BAND_HANG = 0.012, 0.010 * FONT_SCALE


def _wrap(text, max_pt, fontsize, fig):
    """Greedy wrap measured on a real canvas of the final width."""
    r = fig.canvas.get_renderer()
    probe = fig.text(0, 0, "", fontsize=fontsize)
    scale = 72.0 / fig.dpi
    lines, cur = [], ""
    for w in text.split():
        trial = (cur + " " + w).strip()
        probe.set_text(trial)
        if cur and probe.get_window_extent(renderer=r).width * scale > max_pt:
            lines.append(cur); cur = w
        else:
            cur = trial
    if cur:
        lines.append(cur)
    probe.remove()
    return lines


def build(outdir):
    RC = {"font.family": "sans-serif",
          "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans",
                              "TeX Gyre Heros", "Nimbus Sans", "DejaVu Sans"],
          "mathtext.fontset": "stixsans",
          "savefig.dpi": 300, "savefig.bbox": "standard"}
    with matplotlib.rc_context(RC):
        # ---- box width follows the longest line it has to hold ----
        tmp2 = plt.figure(figsize=(FIG_W, 4.0)); tmp2.canvas.draw()
        _r2 = tmp2.canvas.get_renderer(); _sc2 = 72.0 / tmp2.dpi

        def _w(s, fs, **kw):
            t = tmp2.text(0, 0, s, fontsize=fs, **kw)
            v = t.get_window_extent(renderer=_r2).width * _sc2
            t.remove()
            return v

        longest = 0.0
        for s in TITLE:
            longest = max(longest, _w(s, FS_BAN, weight="bold"))
        for blk, fs, kw in ((CARD, FS_CARD, {}),
                            (FIND, FS_RES, {"weight": "bold"})):
            for cell in blk:
                for line in cell.split("\n"):
                    longest = max(longest, _w(line, fs, **kw))
        for s in NRUN:
            longest = max(longest, _w(s, FS_RUN, weight="bold", style="italic"))
        plt.close(tmp2)

        PT = FIG_W * 72.0
        BOX_W_PT = longest + 2 * BOX_PAD_X
        room = (X_RIGHT - X_LEFT) * PT
        max_box = (room + GUTTER_PT) / NCOL - GUTTER_PT
        squeezed = BOX_W_PT > max_box
        BOX_W_PT = min(BOX_W_PT, max_box)
        pitch = BOX_W_PT + GUTTER_PT
        block = NCOL * pitch - GUTTER_PT
        BOX_W = BOX_W_PT / PT
        # the column block is centred, and the conclusion band is then flush
        # with its outer edges
        bx_0 = X_LEFT + ((X_RIGHT - X_LEFT) - block / PT) / 2.0
        CX = [bx_0 + (BOX_W_PT / 2 + i * pitch) / PT for i in range(NCOL)]
        BAND_X0 = CX[0] - BOX_W / 2
        BAND_X1 = CX[-1] + BOX_W / 2

        # ---- bullets wrap to the band, so this must follow the block ----
        band_w_pt = (BAND_X1 - BAND_X0) * PT
        avail = band_w_pt - (2 * BAND_PAD + BAND_HANG) * PT
        tmp = plt.figure(figsize=(FIG_W, 4.0)); tmp.canvas.draw()
        WRAPPED = [_wrap(b, avail, FS_BULLET, tmp) for b in BULLETS if b.strip()]
        _r = tmp.canvas.get_renderer(); _p = tmp.text(0, 0, "", fontsize=FS_BULLET)
        _sc = 72.0 / tmp.dpi
        TEXT_W = 0.0
        for blk in WRAPPED:
            for line in blk:
                _p.set_text(line)
                TEXT_W = max(TEXT_W, _p.get_window_extent(renderer=_r).width * _sc)
        _p.remove(); plt.close(tmp)

        n_title = len(BAND_TITLE.split("\n"))
        title_h = n_title * FS_BAND_TITLE * 1.30
        lead = FS_BULLET * 1.45
        H_BAND = (7.0 + title_h + 7.0 + sum(len(w) for w in WRAPPED) * lead
                  + (len(WRAPPED) - 1) * lead * 0.45 + 8.0)

        FIG_H_PT = (M_TOP + H_BAN + G_BC + H_CARD + G_CR + H_RES + G_RB
                    + H_BAND + M_BOT)
        FIG_H = FIG_H_PT / 72.0

        def y(pt):
            return 1.0 - pt / FIG_H_PT

        B1, B0 = y(M_TOP), y(M_TOP + H_BAN)
        C1 = y(M_TOP + H_BAN + G_BC)
        C0 = y(M_TOP + H_BAN + G_BC + H_CARD)
        F1 = y(M_TOP + H_BAN + G_BC + H_CARD + G_CR)
        F0 = y(M_TOP + H_BAN + G_BC + H_CARD + G_CR + H_RES)
        V1 = y(M_TOP + H_BAN + G_BC + H_CARD + G_CR + H_RES + G_RB)
        V0 = y(FIG_H_PT - M_BOT)
        TIERS = list(zip(TIER_LABELS,
                         [(B0 + B1) / 2, (C0 + C1) / 2,
                          (F0 + F1) / 2, (V0 + V1) / 2]))

        fig = plt.figure(figsize=(FIG_W, FIG_H))
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
        checks = []
        pairs = []
        stubs = []

        def rbox(cx, y0, y1, fc, ec, w=BOX_W, lw=0.9):
            ax.add_patch(FancyBboxPatch((cx - w / 2, y0), w, y1 - y0,
                         boxstyle="round,pad=0.002,rounding_size=0.010",
                         linewidth=lw, facecolor=fc, edgecolor=ec,
                         mutation_aspect=1.0))

        def vconn(cx, ytop, ybot):
            """Dashed shaft with a solid head, one per column between rows.

            Drawn as two elements because a dash pattern applied to a
            FancyArrowPatch also breaks up the arrowhead outline.
            """
            if not SHOW_ARROWS:
                return
            head_h = 0.4 * ARROW_MS / FIG_H_PT      # head length, axes units
            ax.add_line(Line2D([cx, cx], [ytop, ybot + head_h], lw=0.6,
                               color=ARROW_COLOR, linestyle=ARROW_DASH,
                               solid_capstyle="butt"))
            ax.add_patch(FancyArrowPatch((cx, ybot + head_h * 1.2),
                         (cx, ybot), arrowstyle="->",
                         mutation_scale=ARROW_MS, lw=0.6,
                         color=ARROW_COLOR, shrinkA=0, shrinkB=0,
                         joinstyle="miter", capstyle="butt"))

        def funnel(xc, half, ytop, ybot, shaft_pt, head_pt, head_f):
            """Concave-sided funnel with a block head, gradient-filled.

            Used for both the row connectors and the band arrow, so the two
            are the same shape at different scales.
            """
            ws = shaft_pt / (FIG_W * 72.0) / 2.0
            wh = head_pt / (FIG_W * 72.0) / 2.0
            y_sh = ybot + (ytop - ybot) * head_f
            x0, x1 = xc - half, xc + half
            # control points sit above the shaft corners, so each side leaves
            # the top edge flaring outward and arrives at the shaft vertically
            verts = [(x0, ytop),
                     (xc - ws, ytop), (xc - ws, y_sh),
                     (xc - wh, y_sh), (xc, ybot), (xc + wh, y_sh),
                     (xc + ws, y_sh),
                     (xc + ws, ytop), (x1, ytop),
                     (x0, ytop)]
            codes = [Path.MOVETO, Path.CURVE3, Path.CURVE3,
                     Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO,
                     Path.CURVE3, Path.CURVE3, Path.CLOSEPOLY]
            patch = PathPatch(Path(verts, codes), facecolor="none",
                              edgecolor="none", linewidth=0)
            ax.add_patch(patch)
            tip = np.array(to_rgb(ARROW_COLOR))
            t = np.linspace(ARROW_TOP, 1.0, 256)[:, None]
            img = (1.0 - t) + t * tip[None, :]
            im = ax.imshow(img[:, None, :], extent=[x0, x1, ybot, ytop],
                           aspect="auto", interpolation="bilinear",
                           zorder=patch.get_zorder())
            im.set_clip_path(patch)


        for lab, yy in TIERS:
            t = ax.text(LBL_X, yy, lab, ha="center", va="center",
                        fontsize=FS_LBL, color="#555555", weight="bold",
                        linespacing=1.2)
            # the label column runs from the canvas edge to the first box
            checks.append((t, (0.0, 0.0, X_LEFT, 1.0)))

        for cx, title, card, nrun, find, (ban, tint, ec) in zip(
                CX, TITLE, CARD, NRUN, FIND, COLS):
            rbox(cx, B0, B1, ban, ban)
            t = ax.text(cx, (B0 + B1) / 2, title, ha="center", va="center",
                        color="white", fontsize=FS_BAN, weight="bold")
            checks.append((t, (cx - BOX_W / 2, B0, cx + BOX_W / 2, B1)))

            vconn(cx, B0 - ARROW_INSET / FIG_H_PT,
                  C1 + ARROW_INSET / FIG_H_PT)
            rbox(cx, C0, C1, "white", ec)
            t_card = ax.text(cx, C1 - CARD_PAD_T / FIG_H_PT, card,
                             ha="center", va="top", color=INK,
                             fontsize=FS_CARD, linespacing=1.4)
            checks.append((t_card, (cx - BOX_W / 2, C0, cx + BOX_W / 2, C1)))
            t_run = ax.text(cx, C0 + CARD_PAD_B / FIG_H_PT, nrun,
                            ha="center", va="bottom", color=ec,
                            fontsize=FS_RUN, weight="bold", style="italic")
            checks.append((t_run, (cx - BOX_W / 2, C0, cx + BOX_W / 2, C1)))
            pairs.append((t_card, t_run))

            vconn(cx, C0 - ARROW_INSET / FIG_H_PT,
                  F1 + ARROW_INSET / FIG_H_PT)
            rbox(cx, F0, F1, tint, ec)
            t = ax.text(cx, (F0 + F1) / 2, find, ha="center", va="center",
                        color=ec, fontsize=FS_RES, weight="bold",
                        linespacing=1.3)
            checks.append((t, (cx - BOX_W / 2, F0, cx + BOX_W / 2, F1)))

            # the five results converge on one conclusion, so the connector
            # is a manifold: a stub per column into a shared rule, then a
            # single arrow into the band. Drawn after the loop.
            stubs.append(cx)

        # ---- spreading funnel from the result row into the band ----
        if SHOW_ARROWS:
            funnel((CX[0] + CX[-1]) / 2.0,
                   (CX[-1] - CX[0]) * ARROW_SPAN / 2.0,
                   F0 - ARROW_INSET / FIG_H_PT,
                   V1 + ARROW_INSET / FIG_H_PT,
                   ARROW_SHAFT, ARROW_HEAD, ARROW_HEAD_F)

        bx0, bx1 = BAND_X0, BAND_X1
        ax.add_patch(FancyBboxPatch((bx0, V0), bx1 - bx0, V1 - V0,
                     boxstyle="round,pad=0.002,rounding_size=0.010",
                     linewidth=0.9, facecolor=GREY_F, edgecolor=GREY_E,
                     mutation_aspect=1.0))
        t = ax.text((bx0 + bx1) / 2, V1 - (7.0 + title_h / 2) / FIG_H_PT,
                    BAND_TITLE, ha="center", va="center",
                    fontsize=FS_BAND_TITLE, color=INK, weight="bold",
                    linespacing=1.3)
        checks.append((t, (bx0, V0, bx1, V1)))

        block_w = BAND_HANG + TEXT_W / (FIG_W * 72.0)
        x0 = bx0 + ((bx1 - bx0) - block_w) / 2.0
        ypt = 7.0 + title_h + 7.0 + lead * 0.75
        for blk in WRAPPED:
            for k, line in enumerate(blk):
                if k == 0:
                    t = ax.text(x0, V1 - ypt / FIG_H_PT, "\u2022", ha="left",
                                va="center", fontsize=FS_BULLET, color=INK)
                    checks.append((t, (bx0, V0, bx1, V1)))
                t = ax.text(x0 + BAND_HANG, V1 - ypt / FIG_H_PT, line,
                            ha="left", va="center", fontsize=FS_BULLET,
                            color=INK)
                checks.append((t, (bx0, V0, bx1, V1)))
                ypt += lead
            ypt += lead * 0.45

        fig.canvas.draw()
        rr = fig.canvas.get_renderer(); inv = ax.transAxes.inverted()
        bad = []
        for t, (qx0, qy0, qx1, qy1) in checks:
            bb = t.get_window_extent(renderer=rr)
            (a0, b0), (a1, b1) = inv.transform([(bb.x0, bb.y0),
                                                (bb.x1, bb.y1)])
            if (a0 < qx0 - 2e-3 or a1 > qx1 + 2e-3
                    or b0 < qy0 - 2e-3 or b1 > qy1 + 2e-3):
                bad.append(t.get_text().split("\n")[0][:52])

        # text-against-text: box containment alone will not catch a collision
        # between the card body and the run count inside the same box
        for ta, tb in pairs:
            ba = ta.get_window_extent(renderer=rr)
            bb2 = tb.get_window_extent(renderer=rr)
            gap_pt = (ba.y0 - bb2.y1) * 72.0 / fig.dpi
            if gap_pt < 0:
                bad.append(f"collision ({gap_pt:.1f} pt): "
                           f"{tb.get_text()[:28]!r}")

        os.makedirs(outdir, exist_ok=True)
        fig.savefig(os.path.join(outdir, "fig_study_design.pdf"))
        fig.savefig(os.path.join(outdir, "fig_study_design.png"), dpi=200)
        plt.close(fig)

    if bad:
        print(f"    WARNING: {len(bad)} text block(s) overflow their box:")
        for b in bad:
            print(f"      - {b!r}")
    else:
        print(f"    wrote fig_study_design.pdf and .png  "
              f"(all {len(checks)} text blocks fit, band {H_BAND:.0f} pt, "
              f"figure {FIG_W:.2f} x {FIG_H_PT:.0f} pt, "
              f"font x{FONT_SCALE:.2f}, box {BOX_W_PT:.1f} pt)")
        if squeezed:
            print("    NOTE: boxes hit the usable width; reduce FONT_SCALE "
                  "or widen FIG_W")
    return len(bad)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="figures")
    a = ap.parse_args()
    raise SystemExit(1 if build(a.out) else 0)
