#!/usr/bin/env python
#> in   output/figures/_fig3_layout.csv
#> in   output/figures/fig3a_feature_sample_heatmap.pdf
#> in   output/figures/fig3b_subtype_signature_forest.pdf
#> in   output/figures/fig3c_olink_heatmap.pdf
#> in   output/figures/fig3d_nmr_heatmap.pdf
#> in   output/figures/fig3e_disease_heatmap.pdf
#> out  output/figures/Figure3.pdf
import csv
import os
import sys

import fitz  # PyMuPDF

PROJECT_DIR = os.environ.get("PROJECT_DIR", ".")
if not os.path.isdir(PROJECT_DIR):
    sys.exit("PROJECT_DIR does not exist: %s" % PROJECT_DIR)
OUT = os.path.join(PROJECT_DIR, "output", "figures")
CSV = os.path.join(OUT, "_fig3_layout.csv")
DST = os.path.join(OUT, "Figure3.pdf")

MM = 72.0 / 25.4  # mm -> PDF points
ARIAL_BOLD = os.environ.get("ARIAL_BOLD_TTF") or next(
    (c for c in ("/usr/share/fonts/truetype/msttcorefonts/Arial_Bold.ttf",
                 "/Library/Fonts/Arial Bold.ttf",
                 os.path.join(os.environ.get("WINDIR", ""), "Fonts", "arialbd.ttf"))
     if c and os.path.exists(c)), "")

def main():
    if not os.path.exists(CSV):
        sys.exit("missing %s -- source _common/fig3_layout.R first" % CSV)
    with open(CSV, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        sys.exit("empty layout CSV")

    cw = float(rows[0]["canvas_w"])
    ch = float(rows[0]["canvas_h"])
    pt = float(rows[0]["letter_pt"])
    print("[canvas] %.1f x %.1f mm" % (cw, ch))
    if ch > 207.0 + 1e-6:
        sys.exit("canvas %.1f mm exceeds the 207 mm print frame" % ch)

    doc = fitz.open()
    page = doc.new_page(width=cw * MM, height=ch * MM)

    for r in rows:
        src_path = os.path.join(OUT, r["stem"] + ".pdf")
        if not os.path.exists(src_path):
            sys.exit("missing panel PDF: %s" % src_path)
        src = fitz.open(src_path)
        sp = src[0]
        sw, sh = sp.rect.width / MM, sp.rect.height / MM
        ew, eh = float(r["w"]), float(r["h"])
        if abs(sw - ew) > 0.6 or abs(sh - eh) > 0.6:
            sys.exit("%s is %.1f x %.1f mm but the layout expects %.1f x %.1f"
                     % (r["stem"], sw, sh, ew, eh))
        rect = fitz.Rect(float(r["x"]) * MM, float(r["y"]) * MM,
                         (float(r["x"]) + ew) * MM, (float(r["y"]) + eh) * MM)
        page.show_pdf_page(rect, src, 0, keep_proportion=True)
        print("[place ] %-32s %6.1f x %6.1f mm at (%5.1f, %5.1f)"
              % (r["stem"], sw, sh, float(r["x"]), float(r["y"])))
        src.close()

    # panel letters, drawn once, here -- never inside a panel
    if os.path.exists(ARIAL_BOLD):
        page.insert_font(fontname="ArialBd", fontfile=ARIAL_BOLD)
        fname = "ArialBd"
    else:
        print("[warn  ] Arial Bold not found; falling back to Helvetica-Bold")
        fname = "hebo"
    for r in rows:
        page.insert_text(
            fitz.Point((float(r["x"]) + 1.2) * MM,
                       (float(r["y"]) + 1.2) * MM + pt),  # baseline
            r["letter"], fontname=fname, fontsize=pt, color=(0, 0, 0))

    doc.save(DST, garbage=4, deflate=True)
    doc.close()

    # verify what we just wrote
    chk = fitz.open(DST)
    p0 = chk[0]
    n_txt = sum(len(s["text"].strip())
                for b in p0.get_text("dict")["blocks"]
                for l in b.get("lines", [])
                for s in l["spans"] if s["text"].strip())
    imgs = p0.get_images(full=True)
    dpis = []
    for im in imgs:
        rr = p0.get_image_rects(im[0])[0]
        if im[3] > 10:
            dpis.append(im[2] / (rr.width / 72))
    print("\n[assembled] %s  (%.1f x %.1f mm, %.1f MB)"
          % (DST, p0.rect.width / MM, p0.rect.height / MM,
             os.path.getsize(DST) / 1024 ** 2))
    print("[verify ] live text %d chars | embedded rasters %d | raster dpi %.0f-%.0f (need >=800)"
          % (n_txt, len(imgs), min(dpis) if dpis else 0, max(dpis) if dpis else 0))
    chk.close()

if __name__ == "__main__":
    main()
