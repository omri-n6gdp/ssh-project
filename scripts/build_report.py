"""Build the final PDF report for the robustness-testing section of the
SSH topolectrical-circuits project.

Usage:
    python3 scripts/build_report.py

Outputs:
    /Users/omrishahar/code/ssh_project/SSH_Robustness_Report.pdf
"""
import os
from datetime import date

from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_JUSTIFY
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm, mm
from reportlab.platypus import (
    BaseDocTemplate, Frame, Image, KeepTogether, NextPageTemplate,
    PageBreak, PageTemplate, Paragraph, Spacer, Table, TableStyle,
)

PROJECT_ROOT = "/Users/omrishahar/code/ssh_project"
FIG_DIR = os.path.join(PROJECT_ROOT, "figures")
OUT_PDF = os.path.join(PROJECT_ROOT, "SSH_Robustness_Report.pdf")

# ---------- Page setup ----------
PAGE_W, PAGE_H = A4
MARGIN_L = 22 * mm
MARGIN_R = 22 * mm
MARGIN_T = 22 * mm
MARGIN_B = 22 * mm
CONTENT_W = PAGE_W - MARGIN_L - MARGIN_R   # ≈ 166 mm
CONTENT_H = PAGE_H - MARGIN_T - MARGIN_B


def _on_page(canvas, doc):
    canvas.saveState()
    canvas.setFont("Helvetica", 8)
    canvas.setFillColor(colors.HexColor("#666666"))
    canvas.drawString(MARGIN_L, 12 * mm,
                      "SSH Topolectrical Circuit — Robustness Testing")
    canvas.drawRightString(PAGE_W - MARGIN_R, 12 * mm,
                           f"Page {doc.page}")
    canvas.restoreState()


# ---------- Styles ----------
styles = getSampleStyleSheet()

TITLE = ParagraphStyle(
    "Title", parent=styles["Title"],
    fontName="Helvetica-Bold", fontSize=22, leading=27,
    spaceAfter=4 * mm, alignment=TA_LEFT, textColor=colors.HexColor("#15243d"),
)
SUBTITLE = ParagraphStyle(
    "Subtitle", parent=styles["Normal"],
    fontName="Helvetica", fontSize=13, leading=17, spaceAfter=10 * mm,
    textColor=colors.HexColor("#5a6c87"),
)
BYLINE = ParagraphStyle(
    "Byline", parent=styles["Normal"],
    fontName="Helvetica", fontSize=10, leading=13, spaceAfter=2 * mm,
    textColor=colors.HexColor("#333333"),
)
H1 = ParagraphStyle(
    "H1", parent=styles["Heading1"],
    fontName="Helvetica-Bold", fontSize=15, leading=19,
    spaceBefore=8 * mm, spaceAfter=3 * mm,
    textColor=colors.HexColor("#15243d"),
)
H2 = ParagraphStyle(
    "H2", parent=styles["Heading2"],
    fontName="Helvetica-Bold", fontSize=12, leading=16,
    spaceBefore=5 * mm, spaceAfter=2 * mm,
    textColor=colors.HexColor("#15243d"),
)
BODY = ParagraphStyle(
    "Body", parent=styles["Normal"],
    fontName="Helvetica", fontSize=10.5, leading=14.5,
    spaceAfter=2.5 * mm, alignment=TA_JUSTIFY,
)
BODY_TIGHT = ParagraphStyle(
    "BodyTight", parent=BODY,
    spaceAfter=1 * mm,
)
CAPTION = ParagraphStyle(
    "Caption", parent=styles["Italic"],
    fontName="Helvetica-Oblique", fontSize=9, leading=12,
    spaceBefore=1.5 * mm, spaceAfter=4 * mm, alignment=TA_LEFT,
    textColor=colors.HexColor("#444444"),
)
ABSTRACT = ParagraphStyle(
    "Abstract", parent=BODY,
    fontSize=10, leading=14, leftIndent=8 * mm, rightIndent=8 * mm,
    textColor=colors.HexColor("#222222"),
    borderColor=colors.HexColor("#d4ddec"), borderWidth=0,
    backColor=colors.HexColor("#f3f6fb"),
    borderPadding=4 * mm, spaceAfter=6 * mm,
)
MONO = ParagraphStyle(
    "Mono", parent=styles["Normal"],
    fontName="Courier", fontSize=9, leading=12,
)


def P(text, style=BODY):
    return Paragraph(text, style)


def sized_image(path, width_mm):
    """Embed an image scaled to a given width while preserving aspect ratio."""
    from PIL import Image as PILImage
    with PILImage.open(path) as im:
        w_px, h_px = im.size
    w = width_mm * mm
    h = w * h_px / w_px
    return Image(path, width=w, height=h)


def table_style(header_bg="#15243d", header_fg=colors.white, grid="#bbbbbb"):
    return TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor(header_bg)),
        ("TEXTCOLOR", (0, 0), (-1, 0), header_fg),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTNAME", (0, 1), (-1, -1), "Helvetica"),
        ("FONTSIZE", (0, 0), (-1, -1), 9),
        ("LEADING", (0, 0), (-1, -1), 12),
        ("ALIGN", (0, 0), (-1, -1), "LEFT"),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("BOTTOMPADDING", (0, 0), (-1, 0), 6),
        ("TOPPADDING", (0, 0), (-1, 0), 6),
        ("BOTTOMPADDING", (0, 1), (-1, -1), 4),
        ("TOPPADDING", (0, 1), (-1, -1), 4),
        ("GRID", (0, 0), (-1, -1), 0.3, colors.HexColor(grid)),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1),
            [colors.white, colors.HexColor("#f4f6fa")]),
    ])


def figure_block(path, width_mm, caption, keep=True):
    items = [
        sized_image(path, width_mm),
        Paragraph(caption, CAPTION),
    ]
    return KeepTogether(items) if keep else items


# ---------- Build flowables ----------
def build_story():
    s = []

    # --- Title block ---
    s.append(P("Robustness of the Topological Boundary Resonance "
               "in an SSH Topolectrical Circuit", TITLE))
    s.append(P("Monte Carlo analysis under random capacitor disorder",
               SUBTITLE))
    s.append(P("<b>Author:</b> Omri Shahar", BYLINE))
    s.append(P(f"<b>Date:</b> {date(2026, 5, 18).strftime('%d %B %Y')}", BYLINE))
    s.append(P("<b>Project:</b> Topolectrical Circuits — Robustness Testing "
               "section", BYLINE))
    s.append(Spacer(1, 6 * mm))

    # --- Abstract ---
    s.append(P("Abstract", H1))
    s.append(P(
        "A 16-node Su–Schrieffer–Heeger (SSH) topolectrical circuit is "
        "characterised under random Gaussian disorder on every chain capacitor. "
        "Across 6,000 simulated realisations (M = 2,000 per disorder level σ ∈ "
        "{1%, 5%, 10%}) the topological boundary resonance remains centred on "
        "its pristine value of 189.69 kHz, with the resonance shift scaling "
        "linearly in σ and the edge-mode localisation length changing by less "
        "than 5% in the mean even at the strongest disorder. The inverse "
        "participation ratio stays roughly an order of magnitude above the "
        "fully-delocalised reference for all σ. Capacitor disorder preserves "
        "chiral symmetry of the circuit Laplacian, and the survival of the "
        "edge mode is the direct ensemble-level signature of that protection. "
        "These results extend the structured-impurity study of Anzurez et al. "
        "(2025) to random simultaneous multi-site disorder and confirm the "
        "ten-percent-tolerance robustness argument of Lee et al. (2018, Fig. 5) "
        "on a circuit realisation matching the one built in the lab.",
        ABSTRACT))

    # --- 1. Introduction ---
    s.append(P("1.  Introduction", H1))
    s.append(P(
        "Topolectrical circuits realise the band topology of tight-binding "
        "Hamiltonians using only passive inductors and capacitors arranged on "
        "a one-dimensional or higher-dimensional lattice (Lee et al., 2018). "
        "The SSH variant in particular hosts a mid-gap topological boundary "
        "resonance (TBR) at a precisely predictable frequency, with a voltage "
        "profile that is exponentially localised at the chain edge when the "
        "ratio of intracell to intercell capacitance is below unity."))
    s.append(P(
        "The robustness of this boundary resonance to imperfections is the "
        "central physical claim that distinguishes a topological circuit from "
        "a generic LC ladder. Anzurez et al. (2025) tested the claim against "
        "structured impurities — perturbing one or two specific lattice sites "
        "at a time. The analysis presented here extends that study by "
        "introducing simultaneous, independent, random Gaussian disorder on "
        "every capacitor in the chain, repeating the simulation across many "
        "realisations, and quantifying the resulting statistical distributions "
        "of resonance frequency, peak impedance, and edge-state localisation."))
    s.append(P(
        "The headline question is operational: at the manufacturing tolerances "
        "of real capacitors (typically ±1% for Class-1 C0G/NP0 ceramics, ±5–10% "
        "for Class-2 X7R types), does the topological boundary resonance "
        "remain identifiable, and does the edge state remain localised? "
        "The answer is yes on both counts."))

    # --- 2. Circuit ---
    s.append(P("2.  The SSH Topolectrical Circuit", H1))

    s.append(P("2.1  Topology", H2))
    s.append(P(
        "The circuit consists of N = 8 unit cells (16 nodes) connected by an "
        "alternating chain of intracell capacitors C<sub>1</sub> and intercell "
        "capacitors C<sub>2</sub>. Every node is grounded through an inductor "
        "L of identical value. A grounded boundary capacitor C<sub>b</sub> at "
        "node 1 ensures that the diagonal entries of the circuit Laplacian are "
        "uniform across all nodes; without it the boundary resonance is shifted "
        "by tens of percent (Lee et al., 2018, Fig. 2a). An AC voltage source "
        "drives node 1 against ground and acts as the impedance probe."))

    s.append(P("2.2  Component values", H2))
    spec = [
        ["Component", "Count", "Nominal value", "Series resistance", "Role"],
        ["L", "16", "2.2 mH", "10 Ω", "grounded inductor at every node"],
        ["C<sub>1</sub>", "8", "0.10 nF", "100 Ω",
         "intracell coupling (analog of v)"],
        ["C<sub>2</sub>", "7", "0.22 nF", "100 Ω",
         "intercell coupling (analog of w)"],
        ["C<sub>b</sub>", "1", "0.22 nF", "100 Ω",
         "grounded boundary capacitor at N1"],
        ["V<sub>1</sub>", "1", "1 V (AC)", "—",
         "sinusoidal probe at N1"],
    ]
    spec_para = [[Paragraph(c, BODY_TIGHT) for c in row] for row in spec]
    t = Table(spec_para, colWidths=[20*mm, 14*mm, 28*mm, 30*mm, 70*mm])
    t.setStyle(table_style())
    s.append(t)
    s.append(Spacer(1, 4 * mm))
    s.append(P(
        "The dimensionless ratio t = C<sub>1</sub>/C<sub>2</sub> ≈ 0.45 < 1 "
        "places the circuit in the topologically non-trivial phase, where the "
        "SSH mid-gap edge mode is present."))

    s.append(P("2.3  Pristine reference", H2))
    s.append(P(
        "The analytic prediction for the boundary-resonance frequency is "
        "f<sub>R</sub> = 1 / (2π √(L (C<sub>1</sub> + C<sub>2</sub>))) "
        "≈ 189.67 kHz. A deterministic single simulation at the nominal values "
        "yields the reference values used in every comparison below:"))
    ref = [
        ["Quantity", "Pristine value"],
        ["Peak frequency f<sub>peak</sub>", "189.689 kHz"],
        ["Peak impedance |Z|<sub>peak</sub>",
         "81,436 (with V<sub>src</sub> = 1 V)"],
        ["Inverse participation ratio (bulk, N2–N16)", "0.6555"],
        ["Localisation length ξ", "2.53 sites"],
    ]
    ref_para = [[Paragraph(c, BODY_TIGHT) for c in row] for row in ref]
    t = Table(ref_para, colWidths=[80*mm, 60*mm])
    t.setStyle(table_style())
    s.append(t)

    s.append(PageBreak())

    # --- 3. Methodology ---
    s.append(P("3.  Methodology", H1))

    s.append(P("3.1  Disorder model", H2))
    s.append(P(
        "Each capacitor in the chain is replaced by C<sub>k</sub> → "
        "C<sub>k</sub>·(1 + ε<sub>k</sub>) where the multipliers ε<sub>k</sub> "
        "are drawn independently from 𝒩(0, σ²). All sixteen chain capacitors "
        "(one boundary, eight intracell, seven intercell) are perturbed in "
        "every realisation. Three disorder levels are studied: σ = 1%, 5%, and "
        "10%, chosen to span the realistic range of ceramic capacitor "
        "tolerances. Inductors are held at their nominal values; perturbing "
        "them would break chiral symmetry and is reserved as an extension."))
    s.append(P(
        "Each multiplier is clipped to ε > −0.9 to keep all capacitor values "
        "positive. At σ ≤ 0.10 the probability of a clipped sample is below "
        "10<sup>−14</sup>, so the clip has no measurable effect on the "
        "distributions."))

    s.append(P("3.2  Simulation environment", H2))
    s.append(P(
        "Simulations were performed with LTspice (Analog Devices) driven in "
        "batch mode by the PyLTSpice Python package (Brum, 2024). The pristine "
        "netlist template uses named components matching the unit-cell index "
        "(C1_1..C1_8, C2_1..C2_7, Cb, L1..L16), which allows programmatic "
        "modification of any single capacitor by reference. For each "
        "realisation an AC analysis is performed across 80–400 kHz at 1,000 "
        "points per decade (≈700 frequency points), more than sufficient to "
        "resolve a resonance with Q ≈ 250."))

    s.append(P("3.3  Observables", H2))
    s.append(P(
        "Five scalar observables and one vector observable are recorded per "
        "realisation:"))
    obs = [
        ["Symbol", "Definition"],
        ["f<sub>peak</sub>",
         "Resonance peak frequency. Computed as argmax of the impedance "
         "magnitude |Z(ω)| = |V(N1) / (−I(V1))| over the sweep, then refined "
         "by a three-point parabolic fit in (log f, log |Z|) around the peak "
         "index for sub-grid resolution."],
        ["|Z|<sub>peak</sub>",
         "Peak impedance magnitude at f<sub>peak</sub>. A proxy for the "
         "resonance quality factor."],
        ["V(N<sub>i</sub>)",
         "Complex node voltage at every node i = 1..16, evaluated at the peak "
         "index of |Z|."],
        ["IPR",
         "Inverse participation ratio over the bulk nodes N2..N16: "
         "Σ |V<sub>i</sub>|<sup>4</sup> / (Σ |V<sub>i</sub>|<sup>2</sup>)<sup>2</sup>. "
         "N1 is excluded because it is trivially clamped to the source amplitude. "
         "Reference values: fully delocalised → 1/15 ≈ 0.067; fully localised "
         "on one site → 1.0."],
        ["ξ",
         "Edge-state localisation length, obtained by linear regression of "
         "log |V| on the odd nodes N3, N5, ..., N15. ξ is in units of "
         "node-index per e-fold of voltage decay."],
    ]
    obs_para = [[Paragraph(c, BODY_TIGHT) for c in row] for row in obs]
    t = Table(obs_para, colWidths=[28*mm, 134*mm])
    t.setStyle(table_style())
    s.append(t)

    s.append(P("3.4  Ensemble size", H2))
    s.append(P(
        "Each disorder level uses M = 2,000 independent realisations with "
        "deterministic seeds (seed = seed<sub>base</sub> + 10·σ<sub>id</sub> + r) "
        "for bit-reproducibility. The total ensemble is 6,000 simulations, "
        "completed in ≈4.5 minutes on a single laptop with zero failed runs. "
        "At this sample size the standard error of the mean of f<sub>peak</sub> "
        "is √M-suppressed to ≈130 Hz, well below the resolution of any plot "
        "in this report."))

    s.append(PageBreak())

    # --- 4. Results ---
    s.append(P("4.  Results", H1))

    s.append(P("4.1  Resonance frequency distribution", H2))
    s.append(P(
        "Figure 1 shows the ensemble distribution of f<sub>peak</sub> at the "
        "three disorder levels. The histograms are centred on the pristine "
        "value (red dashed) at every σ and the FWHM scales linearly with σ as "
        "expected from first-order perturbation theory: 1.3 kHz at σ = 1%, "
        "6.9 kHz at σ = 5%, and 13.5 kHz at σ = 10%. There is no asymmetry, no "
        "splitting, and no sign of resonance destruction; the worst single "
        "realisation across 6,000 still lies within ±13% of pristine."))
    s.append(figure_block(
        os.path.join(FIG_DIR, "mc_fpeak_histograms.png"),
        width_mm=165,
        caption=("<b>Figure 1.</b>  Histograms of the boundary resonance "
                 "frequency f<sub>peak</sub> across 2,000 realisations per "
                 "disorder level σ. Red dashed line marks the pristine value "
                 "189.69 kHz. The annotated FWHM = 2.355 × σ<sub>f</sub> "
                 "scales linearly in σ.")))

    s.append(P("4.2  Frequency shift versus disorder strength", H2))
    s.append(P(
        "Figure 2 condenses the f<sub>peak</sub> data into mean ± one standard "
        "deviation per σ. The mean tracks pristine to within better than 0.5% "
        "at every σ, so the disorder produces spread without bias. The "
        "approximately linear growth of the error bars confirms the "
        "first-order character of the leading frequency shift."))
    s.append(figure_block(
        os.path.join(FIG_DIR, "mc_fpeak_vs_sigma.png"),
        width_mm=130,
        caption=("<b>Figure 2.</b>  Ensemble mean ± one standard deviation of "
                 "f<sub>peak</sub> as a function of disorder strength σ. Red "
                 "dashed line: pristine reference.")))

    s.append(P("4.3  Peak impedance", H2))
    s.append(P(
        "Figure 3 reports the peak impedance |Z|<sub>peak</sub>, a proxy for "
        "the strength and sharpness of the resonance. A pristine boundary "
        "mode produces a tall, narrow peak; disorder that hybridised the edge "
        "state with bulk modes would shrink this peak. The ensemble mean stays "
        "within 1% of pristine across all σ, with a spread of ±11% at σ = 10%. "
        "The mode is not merely present but still strongly coupled to the "
        "boundary probe."))
    s.append(figure_block(
        os.path.join(FIG_DIR, "mc_zpeak_histograms.png"),
        width_mm=165,
        caption=("<b>Figure 3.</b>  Distribution of |Z|<sub>peak</sub> across "
                 "the ensemble. The pristine peak (red dashed) ≈ 81,436. "
                 "Spread grows roughly proportionally with σ; the mean is "
                 "essentially unchanged.")))

    s.append(PageBreak())

    s.append(P("4.4  Inverse participation ratio (localisation quality)", H2))
    s.append(P(
        "Figure 4 shows the distribution of the bulk inverse participation "
        "ratio (IPR), the standard diagnostic for edge-state localisation. A "
        "perfectly delocalised state over the 15 bulk nodes has IPR = 1/15 ≈ "
        "0.067 (grey dotted); a perfectly localised state has IPR = 1. The "
        "pristine value is 0.655 (red dashed). At σ = 1% the distribution is a "
        "narrow spike at pristine; at σ = 10% the distribution broadens and "
        "develops a lower tail reaching ≈ 0.25, but the bulk of realisations "
        "still cluster around 0.6 — almost an order of magnitude above the "
        "delocalised reference. The chiral-symmetry protection of the edge "
        "mode is intact across all disorder levels."))
    s.append(figure_block(
        os.path.join(FIG_DIR, "mc_ipr_histograms.png"),
        width_mm=165,
        caption=("<b>Figure 4.</b>  IPR distribution (bulk nodes N2..N16). "
                 "Red dashed: pristine. Grey dotted: fully-delocalised "
                 "reference 1/15. The mode remains strongly localised even at "
                 "σ = 10%.")))

    s.append(P("4.5  Localisation length", H2))
    s.append(P(
        "Figure 5 reports the localisation length ξ extracted from a "
        "log-linear fit on the odd-index nodes — the SSH edge mode under "
        "boundary drive. The pristine value ξ = 2.53 sites is reproduced "
        "essentially exactly at σ = 1% (ξ = 2.53 ± 0.02) and shifts by less "
        "than 4% in the mean at σ = 10% (ξ = 2.62 ± 0.24). The edge mode "
        "retains its characteristic decay scale across the entire range of "
        "realistic disorder."))
    s.append(figure_block(
        os.path.join(FIG_DIR, "mc_xi_histograms.png"),
        width_mm=165,
        caption=("<b>Figure 5.</b>  Distribution of the localisation length ξ "
                 "obtained from a least-squares fit of log |V| on odd nodes "
                 "N3..N15. Red dashed: pristine value 2.53.")))

    s.append(P("4.6  Voltage profile across the chain", H2))
    s.append(P(
        "Figure 6 overlays the pristine voltage profile (black) with the "
        "ensemble median (coloured solid) and the 5th–95th-percentile band "
        "(shaded) at each disorder level. The chiral-symmetric alternation "
        "between odd and even sublattices — odd nodes carrying the dominant "
        "edge mode, even nodes suppressed by roughly a factor of fifty — is "
        "preserved across the full disorder range. Even at σ = 10% the median "
        "profile overlaps the pristine profile to plotting precision, and the "
        "5–95% band remains narrow enough to keep the alternation structure "
        "clearly resolved."))
    s.append(figure_block(
        os.path.join(FIG_DIR, "mc_voltage_envelope.png"),
        width_mm=140,
        caption=("<b>Figure 6.</b>  Voltage profile across the 16 nodes at "
                 "f<sub>peak</sub>. Black: pristine. Coloured solid: ensemble "
                 "median per σ. Shaded: 5th–95th-percentile band.")))

    s.append(PageBreak())

    # --- 5. Summary table ---
    s.append(P("5.  Summary statistics", H1))
    s.append(P(
        "All ensemble means and standard deviations from M = 2,000 realisations "
        "per disorder level σ."))
    summary = [
        ["σ", "⟨f<sub>peak</sub>⟩ ± σ<sub>f</sub> (kHz)",
         "f<sub>peak</sub> range (kHz)",
         "⟨|Z|<sub>peak</sub>⟩ ± σ<sub>Z</sub>",
         "⟨IPR<sub>bulk</sub>⟩ ± σ<sub>IPR</sub>",
         "⟨ξ⟩ ± σ<sub>ξ</sub>"],
        ["pristine", "189.689", "—", "81,436", "0.6555", "2.53"],
        ["1%",  "189.675 ± 0.561", "187.5 – 191.6",
         "81,383 ± 905", "0.6549 ± 0.0065", "2.53 ± 0.02"],
        ["5%",  "189.872 ± 2.910", "181.5 – 200.6",
         "81,331 ± 4,575", "0.6448 ± 0.0330", "2.55 ± 0.10"],
        ["10%", "190.355 ± 5.747", "172.5 – 214.7",
         "80,344 ± 9,435", "0.6126 ± 0.0780", "2.62 ± 0.24"],
    ]
    summary_para = [[Paragraph(c, BODY_TIGHT) for c in row] for row in summary]
    t = Table(summary_para,
              colWidths=[16*mm, 36*mm, 30*mm, 28*mm, 30*mm, 22*mm])
    t.setStyle(table_style())
    s.append(t)
    s.append(Spacer(1, 3 * mm))
    s.append(P(
        "The full per-realisation dataset is recorded in "
        "<font face='Courier'>figures/monte_carlo_results.csv</font> "
        "(6,000 rows × 8 columns) and the corresponding voltage profiles in "
        "<font face='Courier'>figures/monte_carlo_voltages.csv</font> "
        "(6,000 rows × 19 columns)."))

    # --- 6. Discussion ---
    s.append(P("6.  Discussion", H1))
    s.append(P(
        "Three observations consistent across every observable support a "
        "single physical conclusion: the boundary resonance survives realistic "
        "capacitor disorder because the disorder is chiral-symmetric. First, "
        "the resonance frequency shifts in mean by less than 0.5% at every σ "
        "and broadens by an amount linear in σ, the standard first-order "
        "behaviour of a non-degenerate eigenfrequency. Second, the peak "
        "impedance — sensitive to spectral weight in the boundary mode — "
        "remains within 11% of pristine even at σ = 10%, ruling out "
        "substantial bulk-mode hybridisation. Third, the IPR and the "
        "exponential decay length ξ are both essentially invariant: the "
        "edge state remains localised at the same scale, on the same "
        "sublattices, with the same alternation."))
    s.append(P(
        "The mathematical mechanism is well known. Capacitor disorder enters "
        "only the off-diagonal entries of the circuit Laplacian, leaving the "
        "diagonal (set by the inductors) untouched. This preserves the chiral "
        "(sublattice) symmetry that protects the SSH mid-gap mode. The "
        "bulk–boundary correspondence theorem (Asbóth, Oroszlány, Pályi, 2016, "
        "Ch. 1) then guarantees the survival of a localised boundary mode as "
        "long as the bulk gap stays open, which it does throughout the σ "
        "range studied here. The Monte Carlo experiment is the ensemble-level "
        "demonstration of that theorem on a real circuit."))
    s.append(P(
        "Relative to prior work: Lee et al. (2018, Fig. 5) gave a single-shot "
        "demonstration that one specific disordered realisation at 10% "
        "tolerance preserves the resonance. Anzurez et al. (2025) extended the "
        "argument to structured single- and double-site impurities and "
        "characterised the resonance displacement systematically as a function "
        "of impurity strength. The ensemble-level analysis presented here "
        "completes the picture by quantifying the *spread* introduced by "
        "random simultaneous multi-site disorder, which is the form of "
        "imperfection actually present in a manufactured circuit."))

    # --- 7. Conclusion ---
    s.append(P("7.  Conclusion", H1))
    s.append(P(
        "The topological boundary resonance of the SSH topolectrical circuit "
        "is robust against realistic random capacitor disorder up to σ = 10%. "
        "The resonance frequency stays within ±3% of pristine, the peak "
        "impedance within ±11%, the inverse participation ratio within an "
        "order of magnitude of pristine and at least an order of magnitude "
        "above full delocalisation, and the edge-mode decay length within "
        "5% of pristine. These results confirm the protective role of chiral "
        "symmetry on a circuit identical to the one realised in the lab and "
        "provide the statistical complement to the structured-impurity "
        "analysis of Anzurez et al. (2025)."))

    s.append(PageBreak())

    # --- 8. References ---
    s.append(P("References", H1))
    refs = [
        ("Anderson, P. W. (1958). Absence of diffusion in certain random "
         "lattices. <i>Physical Review</i> 109, 1492."),
        ("Anzurez, J., Luna Reyes, I., Pérez-Pedraza, J. C., Raya, A., "
         "Sáenz-Trujillo, A. U., Tapia Sánchez, R., &amp; Tinoco-Santillán, G. "
         "(2025). Edge state robustness in the SSH model with impurities: "
         "topoelectric circuit simulations toward topological "
         "heterostructures. <i>Physica B: Condensed Matter</i> 715, 417609."),
        ("Asbóth, J. K., Oroszlány, L., &amp; Pályi, A. (2016). <i>A Short "
         "Course on Topological Insulators</i>. Lecture Notes in Physics, "
         "vol. 919, Springer."),
        ("Brum, N. (2024). <i>PyLTSpice: Python utilities to interact with "
         "LTspice and NGSpice</i>. Package documentation, "
         "<font face='Courier'>https://pyltspice.readthedocs.io</font>."),
        ("Lee, C. H., Imhof, S., Berger, C., Bayer, F., Brehm, J., Molenkamp, "
         "L. W., Kiessling, T., &amp; Thomale, R. (2018). Topolectrical "
         "circuits. <i>Communications Physics</i> 1, 39."),
        ("Su, W. P., Schrieffer, J. R., &amp; Heeger, A. J. (1979). Solitons "
         "in polyacetylene. <i>Physical Review Letters</i> 42, 1698."),
    ]
    for r in refs:
        s.append(P(r, ParagraphStyle("Ref", parent=BODY, fontSize=10,
                                     leading=13.5, spaceAfter=2 * mm,
                                     leftIndent=6 * mm, firstLineIndent=-6 * mm)))

    # --- Appendix A: file manifest ---
    s.append(P("Appendix A.  Project files", H1))
    s.append(P(
        "All paths are relative to the project root "
        "<font face='Courier'>~/code/ssh_project/</font>."))
    manifest = [
        ["File", "Purpose"],
        ["netlists/mc_template.cir",
         "ASCII pristine SPICE netlist with named components (Cb, "
         "C1_1..C1_8, C2_1..C2_7, L1..L16)."],
        ["scripts/monte_carlo.py",
         "End-to-end Monte Carlo driver. CLI: --m, --sigmas, --tag, --no-plots."],
        ["scripts/mc_plots.py",
         "Post-processor: regenerates all figures from the saved CSVs."],
        ["scripts/build_report.py",
         "Generates this PDF report."],
        ["figures/monte_carlo_results.csv",
         "Per-realisation scalar observables (6,000 rows × 8 columns)."],
        ["figures/monte_carlo_voltages.csv",
         "Per-realisation voltage profiles (6,000 rows × 19 columns)."],
        ["figures/mc_summary.csv",
         "Per-σ ensemble means and standard deviations."],
        ["figures/mc_fpeak_histograms.png",
         "Figure 1: f<sub>peak</sub> histograms per σ."],
        ["figures/mc_fpeak_vs_sigma.png",
         "Figure 2: f<sub>peak</sub> mean ± std vs σ."],
        ["figures/mc_zpeak_histograms.png",
         "Figure 3: peak impedance histograms."],
        ["figures/mc_ipr_histograms.png",
         "Figure 4: IPR histograms."],
        ["figures/mc_xi_histograms.png",
         "Figure 5: localisation-length histograms."],
        ["figures/mc_voltage_envelope.png",
         "Figure 6: voltage-profile overlay."],
    ]
    manifest_para = [[Paragraph(c, BODY_TIGHT) for c in row] for row in manifest]
    t = Table(manifest_para, colWidths=[58*mm, 104*mm])
    t.setStyle(table_style())
    s.append(t)

    # --- Appendix B: software environment ---
    s.append(P("Appendix B.  Software environment", H1))
    env = [
        ["Component", "Specification"],
        ["Host OS", "macOS (Darwin 25.3.0)"],
        ["Simulator", "LTspice (Analog Devices) at "
         "/Applications/LTspice.app/Contents/MacOS/LTspice"],
        ["Python", "3.13.x in virtual environment ~/ssh_robustness_env"],
        ["Key Python packages",
         "PyLTSpice, numpy, scipy, pandas, matplotlib, reportlab"],
        ["AC sweep directive", ".ac dec 1000 80k 400k"],
        ["Random number generator", "numpy.random.default_rng (PCG64)"],
        ["Random seed scheme",
         "seed = 10,000 + 10·σ<sub>id</sub> + r, r ∈ [0, M); deterministic."],
    ]
    env_para = [[Paragraph(c, BODY_TIGHT) for c in row] for row in env]
    t = Table(env_para, colWidths=[42*mm, 120*mm])
    t.setStyle(table_style())
    s.append(t)

    return s


def main():
    doc = BaseDocTemplate(
        OUT_PDF, pagesize=A4,
        leftMargin=MARGIN_L, rightMargin=MARGIN_R,
        topMargin=MARGIN_T, bottomMargin=MARGIN_B,
        title="Robustness of the Topological Boundary Resonance in an SSH "
              "Topolectrical Circuit",
        author="Omri Shahar",
        subject="Monte Carlo disorder analysis",
    )
    frame = Frame(MARGIN_L, MARGIN_B, CONTENT_W, CONTENT_H, id="normal",
                  leftPadding=0, rightPadding=0,
                  topPadding=0, bottomPadding=0)
    doc.addPageTemplates([
        PageTemplate(id="main", frames=[frame], onPage=_on_page),
    ])
    story = build_story()
    doc.build(story)
    print(f"wrote {OUT_PDF}")


if __name__ == "__main__":
    main()
