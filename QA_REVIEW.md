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
