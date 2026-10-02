"""
Display helpers for GroupEvaluator's outputs: native's fd-gamma stretch, the geological theme fields, and the
native-style colour key for the theme composites.
"""
import math
from pathlib import Path

import matplotlib
import numpy as np
from PIL import Image, ImageDraw, ImageFont


def fd_stretch(dn):
    """Gamma stretch native applies to fd DN before display (davinci.image.to.gif -gamma), 0-255 float."""
    d = np.asarray(dn, dtype=np.float64)
    return np.clip(d * 7.5 * 0.08 ** np.sqrt(d / 400.0), 0, 255)


GEOLOGICAL_THEME_FIELDS = {"Classification": ("geologic-classifications", "-class"),
                           "material_group": ("geologic-groups", "-group")}


def display_legend(spec, width, single_column=False):
    """
    Native-style colour key (as color.keys/key_*.png) for a COLOUR_THEMES spec, as uint8 (height, width, 3): white bold
    text on black, the theme title, section headings over columns, a dark-to-full gradient swatch per class (flat for
    bins), "python-tetracorder" bottom left. Scaled to the width; native keys are drawn about 1200 wide.
    single_column: the sections stacked one under another in a single column, drawn at the scale of a 600 wide key,
    for composites with the key to the right (tall, narrow images such as core boxes).
    """
    font_file = str(Path(matplotlib.get_data_path()) / "fonts" / "ttf" / "DejaVuSans-Bold.ttf")
    s = (600 if single_column else width) / 1200.0
    font = ImageFont.truetype(font_file, max(9, round(15 * s)))
    small = ImageFont.truetype(font_file, max(8, round(12 * s)))
    sw, pad, gap = round(30 * s), round(14 * s), round(10 * s)
    line_h = round(font.size * 1.25)
 
    # key rows (section, label, [rgb], gradient) in key order, grouped by section
    method = spec["method"]
    rows = []
    if method == "classes":
        for cls in spec["classes"]:
            row = (cls["section"], cls["label"], [cls["rgb"]], True)
            if not rows or rows[-1] != row:                                  # coarse + medium goethite: one row
                rows.append(row)
    elif method == "acid_buffer":
        overlap = spec["overlap"]
        rows += [(c["section"], c["label"], [c["rgb"]], True) for c in spec["generating_classes"]]
        rows += [(overlap["label"], label, [overlap[key]], True) for label, key in
                 (("AGM / ABM < 0.33", "ratio_lt_0.33"), ("AGM / ABM 0.33 - 0.66", "otherwise"),
                  ("AGM / ABM > 0.66", "ratio_gt_0.66"))]
        rows += [(c["section"], c["label"], [c["rgb"]], True) for c in spec["buffering_classes"]]
    elif method == "rgb":
        primaries = {"R": [255, 0, 0], "G": [0, 255, 0], "B": [0, 0, 255]}
        rows += [("", spec["channels"][c]["label"], [primaries[c]], True) for c in "RGB"]
    elif method == "binned":                                                 # bins are full colour: flat swatches
        labels = [b.get("label", "") for b in spec["layers"][0]["bins"]]
        for layer in spec["layers"]:
            rows += [(layer["rule"], b.get("label") or labels[i], [b["rgb"]], False)
                     for i, b in enumerate(layer["bins"])]
    sections = {}
    for section, label, colours, gradient in rows:
        sections.setdefault(section, []).append((label, colours, gradient))
    # fewest rows per column that still gives every column >= 240 px (scaled)
    usable = width - 2 * pad

    per_col = next((r for r in range(3, 100)
                    if sum(math.ceil(len(v) / r) for v in sections.values()) * 240 * s <= usable),
                   max(len(v) for v in sections.values()))
    if single_column:                                                        # one section per block, no splitting
        per_col = max(len(v) for v in sections.values())
    columns = [(name if i == 0 else "", entries[i:i + per_col])
               for name, entries in sections.items() for i in range(0, len(entries), per_col)]
    col_w = usable if single_column else usable / len(columns)

    text_w = col_w - sw - 3 * gap
    spans = {name: math.ceil(len(entries) / per_col) for name, entries in sections.items()}
 
    measure = ImageDraw.Draw(Image.new("RGB", (1, 1)))
    def wrap(text, max_w=text_w):
        lines, line = [], ""
        for word in text.split():
            trial = f"{line} {word}".strip()
            if line and measure.textlength(trial, font=font) > max_w:
                lines.append(line)
                line = word
            else:
                line = trial
        return lines + [line] if line else lines
 
    title = wrap(spec.get("title", ""), usable)
    title_h = line_h * len(title) + gap if title else 0
    headings = [wrap(name, spans[name] * col_w - gap) if name else [] for name, _ in columns]
    head_h = line_h * max(map(len, headings)) + gap if any(headings) else 0
    laid, col_heights = [], []
    for name, entries in columns:
        rows, y = [], 0
        for label, colours, gradient in entries:
            lines = wrap(label)
            h = max(sw, line_h * len(lines))
            rows.append((y, h, lines, colours[0], gradient))
            y += h + gap
        laid.append((name, rows))
        col_heights.append(y)
    
    foot_h = round(small.size * 1.6)
    if single_column:                                                        # blocks stacked, each its own heading
        block_heights = [(line_h * len(h) + gap if h else 0) + y for h, y in zip(headings, col_heights)]
        height = pad + title_h + sum(block_heights) + foot_h + pad
    else:
        height = pad + title_h + head_h + max(col_heights) + foot_h + pad
 
    key = Image.new("RGB", (width, height), "black")
    draw = ImageDraw.Draw(key)
    for i, line in enumerate(title):
        draw.text((width / 2, pad + i * line_h), line, font=font, fill="white", anchor="ma")
    ramp = np.linspace(0.15, 1.0, sw)[None, :, None]                        # dark -> full, left to right
    for c, (name, rows) in enumerate(laid):
        if single_column:
            x0, y0 = pad, pad + title_h + sum(block_heights[:c])
            head_h = line_h * len(headings[c]) + gap if headings[c] else 0
        else:
            x0, y0 = pad + c * col_w, pad + title_h
        for i, line in enumerate(headings[c]):
            draw.text((x0, y0 + i * line_h), line, font=font, fill="white")
        for y, h, lines, rgb, gradient in rows:
            top = y0 + head_h + y
            swatch = np.asarray(rgb, dtype=np.float64)[None, None, :] * (ramp if gradient else 1.0)
            swatch = np.broadcast_to(swatch, (sw, sw, 3)).astype(np.uint8)
            key.paste(Image.fromarray(swatch), (round(x0 + gap), round(top + (h - sw) / 2)))
            text_y = top + (h - line_h * len(lines)) / 2
            for i, line in enumerate(lines):
                draw.text((x0 + sw + 2 * gap, text_y + i * line_h), line, font=font, fill="white")
    draw.text((pad, height - pad - small.size), "python-tetracorder", font=small, fill="white")
    return np.asarray(key)