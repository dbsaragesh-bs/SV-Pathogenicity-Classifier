from pptx import Presentation
from pptx.util import Inches, Pt
F = "reports/figures/"
prs = Presentation(); prs.slide_width = Inches(13.33); prs.slide_height = Inches(7.5)
def slide(title, bullets=None, img=None, img2=None):
    s = prs.slides.add_slide(prs.slide_layouts[5]); s.shapes.title.text = title
    s.shapes.title.text_frame.paragraphs[0].font.size = Pt(30)
    if bullets:
        w = Inches(6.2) if img else Inches(12.3)
        tb = s.shapes.add_textbox(Inches(0.5), Inches(1.5), w, Inches(5.5)); tf = tb.text_frame; tf.word_wrap = True
        for i, b in enumerate(bullets):
            p = tf.paragraphs[0] if i == 0 else tf.add_paragraph(); p.text = "• " + b; p.font.size = Pt(18); p.space_after = Pt(8)
    if img: s.shapes.add_picture(F + img, Inches(6.9) if bullets else Inches(1.5), Inches(1.5), height=Inches(5.2) if not img2 else Inches(4.2))
    return s
slide("Do long reads help classify pathogenic SVs?", ["Matched simulation: same SV, same genome, two technologies", "Arm A short only vs Arm B short + long vs Oracle", "Primary metric: AUPRC with cluster-bootstrap CI", "Result: no gain once an SV is detected; gain is in detection"])
slide("Why simulate", ["Pathogenic SVs have no reads", "Real benign short and long reads come from different material", "Plant every SV (both classes) into a GRCh38 window; simulate both read types with identical label-blind code", "Cost: result depends on simulator realism (checked on HG002)"])
slide("Data: 4,030 SVs, 2,015 per class", ["ClinVar P/LP vs B/LB (>=1 star) + gnomAD-SV v4.1 common PASS", "Conflicts and near-duplicates removed; local-density cap", "Matched on type x size bins (KS p = 0.42)", "Size, type, chromosome, position baselines: AUROC 0.48-0.53"])
slide("Pipeline", ["ART 2x150 + in-house HiFi, 30x, het SV, background variants", "bwa-mem / minimap2 to full GRCh38", "Manta+Delly | Sniffles2+cuteSV", "Detection 0.916 short, 0.945 HiFi; median breakpoint error 1 bp both", "29 pseudo-sample batches; CPU for genomics, GPU for models"])
slide("Leakage audit: pipeline is label-blind", ["Nuisance-only, shuffle, adversarial, shortcut models: AUROC 0.49-0.55", "Source swap (gnomAD-benign -> ClinVar-benign): AUROC 0.89 / 0.96", "Simulation parameters identical across classes (KS p > 0.27)", "High AUPRC (~0.98) comes from gene context, not leakage"])
slide("Main result", ["XGBoost AUPRC: A 0.984, B 0.985, Oracle 0.986", "B - A = +0.0004 [-0.0021, +0.0038]", "LR, RF, MLP agree", "Exon-overlap status changed in 0/3,510 SVs between A and B coordinates"], "pr_curves.png")
slide("Where would a gain be? Ladder and controls", ["Ladder steps all within CI of 0", "Leave-one-family-out: max drop 0.005", "All SVs incl. undetected: +0.0036 [+0.0021, +0.0054]", "Exon-hit only: +0.0019 [+0.0004, +0.0035] (exploratory)"], "ablation_ladder.png")
slide("Why: SHAP and feature families", ["Annotation dominates (exons hit, CDS overlap, LOEUF)", "Long-read evidence adds redundant information", "GPU TreeSHAP on XGBoost"], "shap_families.png")
slide("Stratified gain", ["No stratum excludes 0", "Small SVs: +0.002 [-0.002, +0.010]", "Segdup: -0.017 [-0.040, +0.005] (n = 263)", "Power: MDE ~0.02 for ~250 SVs"], "stratified_forest.png")
slide("Depth sweep", ["Long detection 97.8 / 97.4 / 94.8 / 76.4 % at 30/20/10/5x", "Classifier AUPRC flat", "Detection benefit plateaus at ~10-20x", "60x short: detection 92.9 % vs HiFi 97.8 %"], "depth_sweep.png")
slide("Real-data check (HG002)", ["104 GIAB deletions through the same code", "Adversarial real-vs-sim AUROC 0.97 / 0.98", "Real HiFi breakpoint dispersion 5-10 bp vs 0 in simulation", "Real detection: short(Delly) 0.11, HiFi 0.64", "Simulation is too clean -> optimistic for breakpoints"])
slide("Limitations", ["Simulated errors; one DUP model; het only; DEL-dominant", "Benign mostly gnomAD; easy benchmark (gene context)", "No ONT, INS/INV, phasing, 60x short arm; one replicate", "No real pathogenic reads"])
slide("Conclusion", ["After detection, long reads add no measurable classification value here (CI +/-0.004)", "Their value is detecting SVs short reads miss (94.5 vs 91.6 %)", "Next: real long-read call sets with real breakpoint noise"])
prs.save("reports/SV_ablation_slides.pptx"); print("saved")
