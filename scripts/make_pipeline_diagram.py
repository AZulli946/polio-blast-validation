#!/usr/bin/env python3
"""Slide-ready diagram of the BLAST-DB build + BLAST validation process.

Output is sized for a 16:9 PowerPoint slide (13.33 x 7.5 in).
"""

import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import (
    Circle, Ellipse, FancyBboxPatch, Polygon, Rectangle,
)


# Slide is drawn in inches so axis units == slide inches.
SLIDE_W, SLIDE_H = 13.33, 7.5


def card(ax, cx, cy, w, h, fc):
    ax.add_patch(FancyBboxPatch(
        (cx - w / 2, cy - h / 2), w, h,
        boxstyle="round,pad=0.02,rounding_size=0.18",
        linewidth=0, facecolor=fc,
    ))


def chevron(ax, cx, cy, size=0.40, color="#3F4B66"):
    # (cx, cy) is the visual center; tail is size*0.55 to the left, tip to the right.
    x = cx - size * 0.55
    pts = [
        (x,              cy + size),
        (x + size * 1.1, cy),
        (x,              cy - size),
        (x + size * 0.35, cy),
    ]
    ax.add_patch(Polygon(pts, closed=True, facecolor=color,
                         edgecolor="none", zorder=10))


def database_icon(ax, cx, cy, w=1.2, h=0.85, color="#1F5AA3"):
    ax.add_patch(Rectangle((cx - w / 2, cy - h / 2), w, h,
                           facecolor=color, edgecolor="none", zorder=1))
    ax.add_patch(Ellipse((cx, cy - h / 2), w, h * 0.32,
                         facecolor=color, edgecolor="white",
                         linewidth=2, zorder=3))
    for dy in (-h * 0.15, h * 0.15):
        ax.add_patch(Ellipse((cx, cy + dy), w, h * 0.32,
                             facecolor="none", edgecolor="white",
                             linewidth=1.6, zorder=2))
    ax.add_patch(Ellipse((cx, cy + h / 2), w, h * 0.32,
                         facecolor=color, edgecolor="white",
                         linewidth=2, zorder=4))


def magnifier_icon(ax, cx, cy, r=0.38, color="#6A3F9D"):
    # Lens centered slightly up-and-left so the handle points down-right.
    lens_cx, lens_cy = cx - 0.08, cy + 0.08
    lens_lw = 5
    ax.add_patch(Circle((lens_cx, lens_cy), r,
                        facecolor="white", edgecolor=color,
                        linewidth=lens_lw, zorder=2))

    # Handle starts exactly on the lens edge (at -45 deg) so there is no gap.
    angle = math.radians(-45)
    # Offset by half the lens stroke width (in data units) so the handle
    # starts at the *outer* edge of the lens rim.
    edge_offset = r + (lens_lw / 200.0)  # fig is 13.33 in wide / 100-ish dpu
    sx = lens_cx + edge_offset * math.cos(angle)
    sy = lens_cy + edge_offset * math.sin(angle)
    handle_len = 0.55
    ex = sx + handle_len * math.cos(angle)
    ey = sy + handle_len * math.sin(angle)
    ax.plot([sx, ex], [sy, ey],
            color=color, linewidth=7, solid_capstyle="round", zorder=1)


def check_icon(ax, cx, cy, size=0.45, color="#2E7D32"):
    xs = [cx - size * 0.95, cx - size * 0.20, cx + size * 0.95]
    ys = [cy - size * 0.05, cy - size * 0.65, cy + size * 0.70]
    ax.plot(xs, ys, color=color, linewidth=8,
            solid_capstyle="round", solid_joinstyle="round")


def draw_card(ax, cx, cy, w, h, bg, step_num, icon_fn, icon_kwargs,
              title, title_color, body_lines):
    card(ax, cx, cy, w, h, bg)

    ax.text(cx - w / 2 + 0.40, cy + h / 2 - 0.50,
            f"STEP {step_num}", fontsize=13, fontweight="bold",
            color="#3F4B66", alpha=0.80)

    icon_y = cy + h / 2 - 1.30
    icon_fn(ax, cx, icon_y, **icon_kwargs)

    ax.text(cx, cy + 0.30, title,
            ha="center", va="center", fontsize=16,
            fontweight="bold", color=title_color)

    for dy, txt in zip([-0.55, -1.15, -1.75], body_lines):
        ax.text(cx, cy + dy, txt, ha="center", va="center",
                fontsize=12, color="#222")


def main():
    out_path = Path(__file__).resolve().parent.parent / "results" / "pipeline_diagram.png"
    out_path.parent.mkdir(exist_ok=True)

    fig, ax = plt.subplots(figsize=(SLIDE_W, SLIDE_H), dpi=1200)
    ax.set_xlim(0, SLIDE_W)
    ax.set_ylim(0, SLIDE_H)
    ax.axis("off")
    fig.patch.set_facecolor("white")

    # ---- Subtitle ----
    ax.text(SLIDE_W / 2, 6.95,
            "Curated reference genomes  →  BLAST  →  confirmed poliovirus reads",
            ha="center", va="center", fontsize=15, color="#5C6780", style="italic")

    # ---- Cards ----
    card_w, card_h = 3.7, 6.2
    cy = 3.35
    cx1, cx2, cx3 = 2.15, SLIDE_W / 2, SLIDE_W - 2.15

    card_font_kwargs = dict(
        w=card_w, h=card_h,
    )

    # Chevrons sit centered in the gaps between cards.
    chevron(ax, (cx1 + cx2) / 2, cy, size=0.50)
    chevron(ax, (cx2 + cx3) / 2, cy, size=0.50)

    # Step 1 — Build BLAST database
    draw_card(
        ax, cx1, cy, **card_font_kwargs,
        bg="#E8F1FB", step_num=1,
        icon_fn=database_icon,
        icon_kwargs=dict(w=1.5, h=1.05, color="#1F5AA3"),
        title="Build BLAST database",
        title_color="#1F3B6E",
        body_lines=[
            "Fetch complete poliovirus",
            "genomes from NCBI Entrez",
            "makeblastdb  →  nucleotide DB",
        ],
    )

    # Step 2 — Run BLAST
    draw_card(
        ax, cx2, cy, **card_font_kwargs,
        bg="#EEE5F9", step_num=2,
        icon_fn=magnifier_icon,
        icon_kwargs=dict(r=0.42, color="#6A3F9D"),
        title="Run BLAST",
        title_color="#563D7C",
        body_lines=[
            "Query: polio-aligned reads",
            "blastn against curated DB",
            "keep top hit per read",
        ],
    )

    # Step 3 — Validated hits
    draw_card(
        ax, cx3, cy, **card_font_kwargs,
        bg="#E3F2E1", step_num=3,
        icon_fn=check_icon,
        icon_kwargs=dict(size=0.60, color="#2E7D32"),
        title="Validated hits",
        title_color="#1B5E20",
        body_lines=[
            "Filter by % identity",
            "and genome region",
            "",
        ],
    )

    plt.savefig(out_path, bbox_inches="tight", pad_inches=0.05,
                dpi=1200, facecolor="white")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
