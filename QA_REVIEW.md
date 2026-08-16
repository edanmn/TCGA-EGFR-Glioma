# Skeptical peer review — round 5 (20 passes)

**Scope.** Targeted the three things flagged as untested and never closed: the split-variance
combination promised in limitation twelve, the arbitrary choices in the corruption experiment,
and whether the published figures match their captions. Verification now **383 checks**:
28 assertions, 14 recomputations, 209 revision-figure audits, 137 supplementary audits, minus
overlap. One finding forced a substantial revision of §6.9 for the second time.

## CRITICAL

### C11. Split variability dominates, and it breaks the "all three" claim
`R/38` runs K=6 independent splits per setting and combines within- and between-split variance
by Rubin's rule. Result: **between-split variance is 77–84% of the total.** Which patients fall
in which half matters far more than which genes are sampled. Combined intervals are **3.1–8.2×
wider** (median 6.5×) than the single-split gene bootstrap, and split-level point estimates for
the difference form range **+0.030 to +0.283**.

| Setting | difference form, combined | verdict |
|---|---|---|
| CGGA-693 | +0.063 (−0.036, +0.162) | spans 0 |
| CGGA-325 | **+0.198 (+0.091, +0.306)** | excludes 0 |
| array-301 | +0.170 (−0.008, +0.348) | boundary |

The manuscript claimed the difference form "exceeds its null in all three glioma pairs." With
both variance sources counted, **only CGGA-325 is conclusive**. The direction is consistent —
positive in all three settings and in all 18 split-level estimates — but the paper now claims a
consistent direction rather than three significant results. Abstract, intro, §6.9, Discussion,
Conclusion, limitation twelve and Table 5's caption all revised. **Confidence.** High.

### C12. I mislabelled the variance share, and my own audit caught it
Wrote "between-split variance accounts for 87–92%". That is the **SD** share; the **variance**
share is 77–84%. Caught only because `R/32` asserts paper figures against source CSVs and
failed on two checks. Fixed in four places in the paper and the ambiguous label removed from
`R/38`'s output. This is the audit paying for itself.

## MAJOR

### M15. All three main figures contradicted their captions
Reviewed visually for the first time:
- **Figure 1** subtitle said CGGA is "flat". It is not: it *dips* II→III (−0.015 → −0.055) then
  *rises* III→IV (→ +0.072), with overlapping standard errors. The near-zero Spearman r reflects
  a **non-monotonic** pattern, not a level line.
- **Figure 2**'s caption directs the reader to "the diagonal" — **no diagonal was drawn**. Its
  subtitle asserted flagged genes are "measured inconsistently", the causal reading §6.4 has
  since qualified.
- **Figure 3** drew a dashed line at 0.5 labelled as chance, visually asserting exactly the
  interpretation §6.9 retracts. It showed no matched null.

All three redrawn (`R/39`): Figure 1 with an accurate subtitle, Figure 2 with the identity line
and the fitted attenuation slope (0.55), Figure 3 with per-cohort matched-null levels (0.900,
0.893, 0.793) so the correct comparison is the visible one. Captions rewritten to match.

## Passes returning clean
Threshold grid re-verified; determinism; all 6 tabular column specs; no dangling cross-refs;
supplementary audit 137/137; rounding; provenance; reference completeness.

## Confidence after round 5
**High on arithmetic; the central claim is now materially weaker but honestly stated.**
Rounds 3, 4 and 5 each found that an inference in the revision was too confident — replicate
t-tests, then an anticonservative gene bootstrap, then an interval omitting the dominant
variance component. Each correction moved the same direction: the evidence is weaker than the
previous version claimed. That the direction is consistent across all 18 split-level estimates
is the durable part; the significance of any individual cohort pair is not.

The honest summary is that **cohorts of a few hundred patients per arm cannot resolve
differences of this size with a matched-null design.** That is a statement about the data, not
the method, and no further re-analysis changes it.

---

# Skeptical peer review — round 4 (20 passes)

**Scope.** Rounds 1–3 are clean on re-run, so these 20 passes went after the largest
untested surface: whether the conclusions survive the analyst choices made arbitrarily,
plus the supplementary document (never machine-audited) and the validity of the
revision's own intervals. Verification totals now **362 checks**: 28 assertions,
14 recomputations, 183 revision-figure audits, 137 supplementary-table audits.

## CRITICAL

### C9. The conclusions are anchor-dependent, and only one anchor was calibrated
`R/36` varies FDR ∈ {0.01, 0.05, 0.10} × replication p ∈ {0.01, 0.05, 0.10} × anchor ∈
{grade, IDH} × metric, across the three glioma settings: 108 cells.

- **Thresholds: fully robust.** All **54 of 54** grade-anchor cells reproduce the reported
  sign. The conclusions are not an artifact of our cutoffs.
- **Anchor: not robust.** Only **39 of 54** IDH-anchor cells agree, and CGGA-325's product
  form inverts in **all 9** (0% agreement).

§6.9's calibration used the grade anchor exclusively, so Table 4's IDH column is reported
**uncalibrated**. Added to §6.9 and as limitation fourteen. **Confidence.** High.

### C10. The gene bootstrap was anticonservative *(round-4 carryover, fixed)*
Genes are co-expressed; resampling them independently treats correlated observations as
independent. `R/35` re-runs as a cluster bootstrap over 200 co-expression clusters:
intervals widen by a median **2.2×** (range 1.6–2.7×). Every conclusion survives. Reported
in §6.9 and limitation twelve rather than silently substituted.

## MAJOR

- **M12. The supplement had never been machine-checked.** Tables S1–S8 carry ~60 hand-
  transcribed numbers with no audit (R/21 covers only S9; R/22 only the main tables).
  `R/37_supp_audit.R` added, gated in `run_all.R`: **137/137 pass.** No errors found — but
  the exposure was real for four review rounds.
- **M13. Split-to-split variability is not in the intervals.** Two independent splits gave
  array-301 difference-form estimates of +0.148 and +0.227. Limitation twelve now says so.
- **M14. README was stale for 13 scripts** and is the repo's entry point. Rewritten with the
  full pipeline, all four verification suites and their expected totals.

## MODERATE
- **N13. Duplicate result files for the same quantity.** `pathway_screen.csv` (per-unit) vs
  `pathway_screen_std.csv` (per-SD, the one Table S2 uses); `cgga_positive_controls.csv` vs
  `panel_correlations.csv`. The manuscript uses the correct file each time, but a reader
  opening the wrong one gets different numbers. `R/37` now prints both discrepancies
  explicitly. Recommend deleting the superseded files before release. *Not applied.*
- **N14.** Audit-script bug found and fixed during this round: the pathway-screen lookup
  matched three cohorts and returned NA. Caught because the check failed loudly.

## Passes returning clean
Determinism (R/33 reproduces **bit-identically**); rounding uniformly 3 d.p. across all six
main tables; provenance strings intact (GDC 45.0, CGGA 20200506, R 4.6.0, Zenodo DOI); every
§6.9 figure traces to a source result file; no dangling § or table refs; all 29 references
cited; balanced math delimiters; 9 supplementary tables and 5 supplementary figures all
defined and cited; paper–supplement numeric agreement on shared quantities.

## Confidence after round 4
**High on arithmetic and internal consistency; moderate on the central claim, now with a
sharper boundary.** The threshold grid materially raises confidence — 54/54 is a strong
robustness result that no prior round had established. The anchor sensitivity lowers it in a
specific, stated way: the calibration is a grade-anchor result. Residual risk remains
concentrated in generalization, not correctness.

---

# Skeptical peer review — EGFR/glioma positive-control paper (round 3)

**Round 3 scope.** Round 2's ten passes are now clean, so repeating them would be theatre.
These ten passes target the *revision's own new experiments*, which no one had reviewed,
plus a literature check of the novelty claim. Verification totals: 28/28 assertions,
14/14 recomputations, **168/168 revision-figure audits** (after the N12 follow-up).

---

## CRITICAL (round 3)

### C6. The p-values in Table 5 were invalid *(self-inflicted, now fixed)*
**Passage (prior).** "the real cohort significantly exceeds the null (0.942 vs 0.904, p=0.010)."
**Issue.** Those were Welch t-tests across 8 *resampling replicates* of one fixed cohort.
Replicate spread is Monte-Carlo error of the estimate, shrinking as 1/√reps, so p tracks how
long the simulation was run. Demonstrated directly: re-running with 3,4,…,8 replicates moves
p from **0.32 → 0.010** while the effect moves only +0.030 → +0.037.
**Fix.** `R/33_gene_level_inference.R` — paired bootstrap over **genes** (the AUC's actual
sampling unit, and both arms are scored on the same genes), 2,000 resamples, BH within metric.
Table 5 rebuilt on that basis; `R/28` marked superseded.
**Confidence.** High.

### C7. The corrected inference changes a headline claim — then N12 changed it back
With valid intervals at the 2,000-gene universe, the difference form exceeded its null in **two
of three** glioma pairs: CGGA-325 +0.268 (+0.219,+0.319) and array-301 +0.218 (+0.162,+0.273),
but CGGA-693 +0.051 (−0.019,+0.121). **Superseded by the N12 run below.** On Table 4's own
universe (`R/34`) CGGA-693 becomes **+0.138** (+0.102,+0.176) and the difference form separates
in **all three** glioma pairs. The manuscript now reports the full universe as primary and the
2,000-gene run as the robustness check, flagging that CGGA-693's result is universe-dependent.

### C8. Four nearest bodies of prior art were uncited
Found by literature search, not introspection: **sigQC** (signature-level cross-dataset QC),
**RUV** (negative control genes — the mirror construct), the **replication-prediction
literature** (effect size predicts replication at AUC≈0.77 in psych/econ — i.e. the baseline
we omitted is the field's *established* one), and the **TCGA→CGGA validation** design. Added as
refs 26–29 with a positioning paragraph in §2.

---

## MAJOR (round 3)

- **M7. Breast BRCA→METABRIC is not evaluable** under a matched null: splitting 622 patients
  with 67 events leaves **26** prognostic genes. Previously reported as if evaluable. Now stated
  in Table 5 and the Discussion; the scope claim rests on **four** cohort pairs, not five.
- **M8. No multiplicity correction** across the 30 metric×setting cells. BH now applied; 17 of
  18 nominally significant cells survive.
- **M9. Self-correlation leakage** in the co-expression metric — hub genes were drawn from the
  evaluated set, so 100/2,000 genes included their own r=1 in both cohorts. Fixed in `R/33`.
- **M10. "The exact signature EGFR shows"** was an overstatement. Three modes (compression 0.50,
  permutation 0.50, noise k=2) land inside EGFR's 0.23–0.42 variance-ratio band. The
  discriminator is *within*-stratum variance, which compression preserves and the other two
  inflate; EGFR's is unchanged (1.84 vs 1.88). Claim rewritten as an explicit two-statistic
  match, flagged as a signature match rather than a known cause.
- **M11. `R/26` bootstrap CIs were unseeded** — not reproducible. Seed added.

## MODERATE (round 3)
- **N7.** Table 5 (2,000 genes, discovery n=368–1,066) vs Table 4 (6,512 genes, n=793) were
  compared as if interchangeable. Both now state their design; text flags the within-row
  comparison as the valid one.
- **N8.** §4.6 said "steps 00–29" when scripts run to 33. Fixed.
- **N9.** Figure S1's pipeline diagram omitted the revision analyses. Box added.
- **N10.** Four generated figures orphaned; documented as repo artifacts (page budget).
- **N11.** The effect-size increment (+0.035 to +0.199) was described in Methods and asserted by
  the audit but **absent from Results** after the §6.9 rewrite. Restored.
- **N12. RESOLVED.** The matched null was re-run on Table 4's universe (`R/34`; 8,000-gene rule,
  5,495–6,381 prognostic after the split). Signs held in **15 of 16** comparisons; the single flip
  is the co-expression metric in CGGA-693, not a headline claim. The difference form strengthens
  to 3 of 3 glioma pairs, and the product form is shown **inverted** in CGGA-693 (−0.075) where
  the difference form is clearly positive (+0.138) — the sharpest form of the functional-form
  argument. BRCA→METABRIC re-tested and confirmed unevaluable (**0** prognostic genes of 8,000).
  Table 5 rebuilt on the full universe; this also retires N7, since Tables 4 and 5 now share a
  gene universe.

## Passes that came back clean
Dependency-order integrity across all 32 pipeline steps (0 violations); null/discovery
disjointness; all 15 abstract numerics traceable to Results; all 29 references cited; no dangling
§ or table refs; all tabular column specs; no duplicated text; balanced math delimiters; no
surviving unsupported universals.

## Confidence after round 3
**High on arithmetic and internal consistency; moderate on the central claim.**
The claim is now supported by valid inference rather than by an artifact, which is the main
reason to trust it over round 2's version. It rests on four evaluable cohort pairs, three from
one consortium, with one glioma pair showing no effect and one breast direction unevaluable.
The ground-truth experiment remains the strongest component; its link to EGFR is a
two-statistic signature match, not a demonstrated cause.

---

# Round 2 (superseded where round 3 revises it)

# Skeptical peer review — EGFR/glioma positive-control paper (round 2)

**Scope.** Every numeric claim in the manuscript was re-derived from the cached data,
either by the pre-existing suites (`R/21`, `R/22`) or by scripts written for this review
(`R/23`–`R/32`). Nothing was accepted because it appeared in the paper.

**Verification totals.** 28/28 cohort assertions, 14/14 statistical recomputations,
**120/120 revision-figure audits** (`R/32_paper_audit.R`, now gated in `run_all.R`).
Power calculations checked by hand against the Schoenfeld relation: √(7.849/253)→HR 1.19,
√(7.849/127)→1.28, binary at p=0.55→1.65. All correct.

> **Round-1 self-correction.** The previous review concluded that "a same-population null
> scores higher than a real cohort, therefore the gate detects nothing." That conclusion
> rested on **one** cohort pair and did not survive testing in five. It was an
> overgeneralization of exactly the kind this review is meant to catch, and it has been
> retracted and replaced throughout. See C1.

---

## CRITICAL

### C1. The round-1 calibration claim was itself overgeneralized *(new, self-correcting)*
**Section.** 6.9, Abstract, Discussion, Conclusion.
**Prior claim.** "Against a matched null the same score reaches AUC 0.94, *higher* than
against the real cohort … a correction that applies to any cross-cohort quality metric."
**Issue.** Tested across five cohort pairs, two tumor types, three platforms, six metrics
(`R/30`, 8 replicates each). The picture is mixed, not uniform:

| Setting | product null | product real | diff |
|---|---|---|---|
| TCGA→CGGA-693 | 0.931 | 0.874 | **−0.058** (p=8e-4) |
| TCGA→CGGA-325 | 0.904 | 0.942 | **+0.037** (p=0.010) |
| TCGA→array-301 | 0.888 | 0.935 | **+0.047** (p=6e-5) |
| BRCA→METABRIC | 0.739 | 0.714 | −0.025 (ns) |
| METABRIC→BRCA | 0.717 | 0.582 | −0.135 (p=2e-4) |

The real cohort significantly **exceeds** its null in two of three glioma pairs. Rate-matching
the null arm to the real arm's replication rate changes nothing (0.931→0.921), so this is
not a base-rate artifact.
**Correction.** Applied: §6.9 rewritten, Table 5 added, abstract/intro/discussion/conclusion
and limitation six restated. What survives is that the *product-form* score is heavily
confounded (null AUC 0.89–0.93) — not that the gate detects nothing.
**Confidence.** High.

### C2. The confound is in the functional form, and is fixable *(new, positive)*
**Finding.** The product `r_D × r_R` is large when *both* correlations are large, so it
rewards anchor magnitude — which tracks effect size — as much as agreement. The difference
form `−|r_D − r_R|` scores only disagreement:

| | product | difference |
|---|---|---|
| null AUC (glioma) | 0.89–0.93 | **0.60–0.62** |
| real − null, CGGA-693 | −0.058 | **+0.048** |
| real − null, CGGA-325 | +0.037 | **+0.189** |
| real − null, array-301 | +0.047 | **+0.154** |

The difference form is nearly unconfounded and beats its null in **all three** glioma pairs.
**Why it matters.** This converts the round-1 negative result into an actionable
recommendation: change the statistic, don't abandon the method.
**Correction.** Applied — §6.9, Table 5, and the recommendation in the Conclusion.
**Confidence.** High — `R/30_multisetting_null.R`.

### C3. Negative control behaves correctly, licensing the above
**Finding.** The discovery effect size `|β_D|`, which by construction cannot carry
cross-cohort information, favors the null in **all five** settings (−0.061 to −0.167). That
it separates in the expected direction while the difference-form anchor separates in the
opposite one is what licenses reading the latter as genuine cross-cohort signal.
**Confidence.** High.

### C4. The effect-size baseline was missing for glioma *(carried from round 1)*
|β_TCGA| gives 0.738/0.770/0.783 (unadjusted) and 0.645/0.651/0.707 (adjusted) — far above
the "0.36–0.59" the reported baselines suggest. It is computed in the breast analysis
(`R/17:84`) but never in glioma. The anchor retains incremental value (+0.035 to +0.199,
LR p<10⁻³⁷). §6.8's rebuttal of the proxy objection was a non-sequitur (it tests noise, not
signal). **Applied.**

### C5. "Reproduced across three datasets" does not survive control *(carried)*
True of the marginal grade correlation, which is composition-confounded. Holding IDH fixed,
CGGA-693's EGFR–grade slope (+0.124) matches TCGA's (+0.115); within-IDH-wildtype the
correlation is positive in all three CGGA sets (+0.163/+0.065/+0.079). Only CGGA-693's
reversed IDH slope is composition-immune (p=8.8e-4), and it does not reproduce in CGGA-325
(p=0.99) or array-301 (p=0.29). **Applied.**

---

## MAJOR

### M1. Detectability is strongly mode-dependent, with a real blind spot *(new)*
Five mechanistically distinct failure modes injected at known doses (`R/31`, ground truth =
genes actually corrupted):

| Mode | AUC range | Note |
|---|---|---|
| Compression (shrink between-stratum signal) | 0.639–0.848 | **best detected; EGFR's signature** |
| Permutation | 0.630–0.858 | |
| Additive noise | 0.567–0.779 | |
| Contamination (paralog mixing) | 0.543–0.734 | |
| **Flooring** (detection limit / saturation) | **0.498–0.599** | **blind spot** |

Replicate SDs 0.004–0.026. Flooring is near-invisible even where corrupted genes' replication
falls to 0.82 vs 0.95 for clean genes.
**Why it matters.** Round 1 tested permutation only — which is *not* what EGFR shows. This
closes that gap: compression at λ=0.50 leaves a variance share of 0.136 against a 0.339
baseline (ratio 0.40, inside EGFR's observed 0.23–0.42 band) and is detected at AUC 0.730.
The gate could have caught the EGFR failure specifically. It also documents a failure mode it
cannot catch. **Applied** — Table 6, §6.9, limitation ten.

### M2. Internal inconsistency between Table 6's caption and text *(new, found by audit)*
The caption states an uncorrupted baseline of 0.339; the text's ratio (0.393) was computed
against a condition-specific baseline. A reader dividing as the caption instructs gets 0.401.
**Corrected** to 0.40 with both numbers stated explicitly. Caught only because `R/32` asserts
transcribed figures against source CSVs. **Confidence.** High.

### M3. Table 6 column-spec error *(new)*
`\begin{tabular}{llccccc}` declared 7 columns for a 6-field table. LaTeX compiled it silently
with a stray empty column. **Fixed** to `{llcccc}`; a column-count check now covers all six
tables (all pass).

### M4. Two §6.10 numbers were unreproducible *(carried)*
"basal ≤0.54" and "r=0.84" came from no committed script. Recomputed: 0.542 and 0.839 — both
**correct**, not hallucinated. `R/29_breast_verify.R` added. Also found: in the reverse
direction the effect-size baseline (0.574) **beats** the ER anchor (0.556). **Applied.**

### M5. Composition explains a minority genome-wide *(carried)*
Transporting TCGA's within-stratum structure onto CGGA's stratum weights: composition accounts
for only 8–25% of mean anchor discordance, and the transported prediction does not beat TCGA's
raw correlations (r=0.815–0.868 vs 0.812–0.882). EGFR's variance share collapses 0.106 →
0.045/0.030/0.024 with within-stratum SD unchanged (1.88 vs 1.84). **Applied.**

### M6. Table S7 confounds two changes *(carried, NOT applied)*
Its "subtype-adjusted" row changes discovery set *and* replication criterion at once on 183
genes. The 2×2 decomposition (`R/23`) isolates them: holding the gene set at 6,512 and changing
only the criterion gives 0.700/0.782/0.741 with tight CIs. **Recommend replacing Table S7.**

---

## MODERATE

- **N1. Spearman/Pearson split.** Screens use Spearman, transport/baseline use Pearson
  (needed for the closed-form transport). CGGA-693 anchor AUC 0.815 vs 0.808. Internally
  consistent within each analysis; must not be pooled. **Applied** (limitation eleven).
- **N2. Stale outputs contradict the paper.** `results/cgga_positive_controls.csv` (from
  `R/09`) gives EGFR r=0.183, PTEN=−0.384; `panel_correlations.csv` (from `R/20`, which the
  paper uses) gives 0.179, −0.357. `_cgga_diag.log` reports grade HR 8.50 vs Table 2's 8.0.
  Paper is self-consistent; the superseded artifacts should be deleted. **Not applied.**
- **N3. `n=422` overstates the CGGA Cox sample** in `R/23`/`R/26` logs (effective n=403 after
  `coxph` drops missing OS). Affects no AUC. **Not applied.**
- **N4. Breast rate-matched arm is degenerate** in the BRCA→METABRIC direction (only 59
  prognostic genes); reported as NA rather than interpreted. Correct handling, worth a note.
- **N5. Page budget.** Now **14 pages** (was 10). The 10-page target was the author's, not a
  venue's. Levers: move Table 5 or Table 6 to the supplement, or Figure 2 (content is in
  Table S6). **Not applied** — depends on venue.
- **N6. Corruption doses are not calibrated to any real assay**, and batch/annotation failure
  modes are untested. **Applied** as limitation ten.

## Verified correct (no action)
§6.7's ten screen numbers exact (88.0% QC-pass, 92.8% vs 24.3%, Fisher p=4.87e-21, OR=39.2,
no-gating 87.3%, 2/2 sign flips, three threshold rows). §6.8 de-duplication (264−67=197, AUC
0.91 vs 0.92, 0.684→0.578). Tables 1, 3, S9, all cohort sizes. All Table 1/3 HRs, p-values,
C-indices. CGGA-693 EGFR means 4.17 vs 3.63. recount3 8.92/7.99 and 9.20/7.63. All six tables'
column specs. All section cross-references resolve (§4.6, §6.2, §6.4, §6.6–§6.10 all exist).

## Potential hallucinated content
**None.** Every number resolved to a real computation. M4's two figures were unreproducible
but proved correct on recomputation.

## Claims not verifiable from the paper alone
1. The mechanism of EGFR's attenuation in CGGA — raw reads are controlled-access. M5 narrows
   it to signal compression; M1 shows that mode is detectable; neither identifies the cause.
2. Generalization beyond the five cohort pairs tested. The calibration problem is now shown to
   be *setting-dependent*, which is itself the finding.
3. Patient overlap between CGGA datasets rests on CGGA's own identifiers.

## Prioritized fix list
1. **C1** — retract the round-1 overgeneralization. *Applied.*
2. **C2/C3** — adopt the difference form; report the negative control. *Applied.*
3. **C4/C5, M1, M2, M3, M4, M5, N1, N6** — *Applied.*
4. **M6** — replace Table S7 with the 2×2. *Not applied* (supplement untouched).
5. **N2, N3** — delete superseded outputs; fix printed cohort sizes. *Not applied.*
6. **N5** — page budget, after venue selection. *Not applied.*

## Confidence that the paper is scientifically sound after revision

**High for the descriptive results; moderate-to-high for the methodological claim.**

162 independent checks pass (28 + 14 + 120). The EGFR case study (Tables 1, 3) is solid and
confirmatory. The revised methodological claim is now *tested* rather than asserted, and the
test partly contradicted the previous round's conclusion — which is the main reason to trust
the current version more than either predecessor.

Two things I would still expect a reviewer to press. First, the difference-form result rests
on three glioma pairs drawn from one consortium; the breast pair shows nothing for any metric,
so "works in glioma, silent in breast" is two data points about tumor types, not a pattern.
Second, the ground-truth experiment establishes detectability against *simulated* corruption;
the inference that real CGGA EGFR resembles the `compress` arm rests on a variance-signature
match, which is suggestive rather than conclusive.

Residual risk is concentrated in generalization, not in correctness or arithmetic.
