# Manuscript result edit checklist

This checklist compares:

- `Main Text without saved changes.docx`
- `Supporting Information.docx`

with the corrected publication run in `manuscript_outputs/`. The corrected run
uses Camelot camera-local timestamps (Asia/Jerusalem) and passed all validation
checks.

## Main text

### Results: detections and naive occupancy

The independent-detection totals and land-use proportions are unchanged.

Change:

- European badger naive occupancy: `0.85` to `0.83` (50/60 cameras).
- Human naive occupancy: `0.88` to `0.95` (57/60 cameras).

Source: `supplementary_tables/Table_A2_detection_summary_and_naive_occupancy.csv`.

### Results: diel activity shifts among land uses

Retain the conclusion that land use significantly affected five of seven
species. Update the detailed counts and low-sample result:

- Golden jackal: 6/6 Holm-significant contrasts.
- Wild boar: 6/6 Holm-significant contrasts.
- Red fox: 5/6 Holm-significant contrasts.
- Mountain gazelle: 4/6 Holm-significant contrasts.
- Indian crested porcupine: 3/6, not 4/6, Holm-significant contrasts.
- Grey wolf: High culling versus Nature reserves, Delta = 0.883,
  p = 0.0394, adjusted p = 0.236; High culling versus Minefields,
  Delta = 0.886, p = 0.0404, adjusted p = 0.236. Neither is significant
  after Holm adjustment.
- European badger: no nominal or Holm-significant contrast. All raw
  permutation p-values are at least 0.19 and all adjusted p-values are 1.

Replace the submitted Table 1 with
`main_tables/Table_1_land_use_activity_tests.xlsx` or its CSV equivalent.
The replacement table identifies whether Delta-hat 1 or Delta-hat 4 was used
and identifies the uncertainty source for every row.

### Results: nocturnality

The direction, effect sizes, and biological conclusion are unchanged. Update:

- Community slope p-value: `9.5 x 10^-19` to `1.90 x 10^-18`.
- Describe the badger slope as negative but imprecise, not positive:
  beta = -0.102, SE = 0.100, p = 0.312.
- The gazelle response is statistically supported rather than only a
  "tendency": beta = +0.119, SE = 0.054, p = 0.0271. With the manuscript's
  land-use coding, this means gazelles were less nocturnal under greater
  management disturbance.

Source:
`additional_results/part2_nocturnality_gradient_slopes.csv`.

### Results: pairwise temporal overlap

The directions and biological conclusions are unchanged, but the submitted
rank-model coefficients were approximately tenfold too large. Replace them
and the associated tests:

- Non-gazelle pairs: beta = +0.136, SE = 0.0359, p = 0.000149.
- Gazelle-other pairs: beta = -0.132, SE = 0.0346, p = 0.000136.

Correct unranked estimated marginal means:

- Non-gazelle pairs: Minefields = 0.821, Nature reserves = 0.810,
  Low culling = 0.821, High culling = 0.879.
- Gazelle-other pairs: Minefields = 0.567, Nature reserves = 0.525,
  Low culling = 0.471, High culling = 0.476.

Correct the two gazelle contrasts quoted in the Results:

- Minefields versus Low culling: p = 0.000830.
- Minefields versus High culling: p = 0.001765.

Source: `supplementary_tables/Table_A5_pairwise_overlap_beta_models.csv`.

### Results: gazelle response to local human activity

The submitted primary result is correct and should be retained:

- Q4 minus Q1 nocturnality = +1.0678 h.
- 95% CI = 0.6900 to 1.4362 h.
- p = 2.03 x 10^-8.
- Temporal overlap Delta = 0.8096; kernel randomization p < 0.001.
- Q1 n = 960 and Q4 n = 483 independent detections.

Source: `additional_results/gazelle_human_activity_summary.csv`.

In the opening sentence of this subsection, replace the reference to
`Table A.4` with `Tables A.1 and A.3`: Table A.1 contains the human-activity
contrast and Table A.3 contains the land-use nocturnality summaries.

### Results: gazelle reserve comparison

The overlap and sample sizes are unchanged after rounding. Change:

- Kernel randomization p-value: `0.0007` to `0.0014`.

Retain Delta = 0.85 and n = 992 versus 214 independent detections.

Source: `additional_results/part6_gazelle_reserve_overlap_test.csv`.

## Figure captions and displayed figures

### Figure 3

Replace the submitted figure with
`main_figures/Figure_3_land_use_activity.png` or the PDF.

The corrected wolf panel is High culling versus Nature reserves. Remove the
caption statement about asterisks because the corrected figure has no panel
asterisks. Describe panels as selected land-use contrasts; the badger panel
retains the specified High culling versus Low culling contrast.

### Figure 6

Replace the submitted figure with
`main_figures/Figure_6_gazelle_human_and_reserve_activity.png` or the PDF.

In the Panel B caption, replace "peaking near midday" for the temporally zoned
reserve with "peaking in mid-afternoon and declining before sunset." Human
activity in the open-access reserve peaks in late afternoon near sunset.

## Supporting information

Replace these tables with the generated versions:

- Table A.1: `supplementary_tables/Table_A1_species_responses_to_human_activity.csv`.
  Key changes include Q4 wild-boar n = 886, total wild-boar detections =
  4,333, and porcupine kernel p = 0.109
  (not significant), and badger kernel p = 0.841.
- Table A.2: `supplementary_tables/Table_A2_detection_summary_and_naive_occupancy.csv`.
  Key changes are badger occupancy = 0.83 and human occupancy = 0.95.
- Table A.3: `supplementary_tables/Table_A3_nocturnality_by_species_land_use.csv`.
  Replace all submitted confidence intervals; the point means are unchanged.
- Table A.4: `supplementary_tables/Table_A4_pairwise_overlap_by_land_use.csv`.
  The replacement records the estimator and sample sizes and uses Delta-hat 1
  for low-event badger comparisons.
- Table A.5: `supplementary_tables/Table_A5_pairwise_overlap_beta_models.csv`.
  Use Delta, not Delta-hat 4, in its caption because the input includes both
  Delta-hat 1 and Delta-hat 4 coefficients.

Replace Figures A.1 and A.2 with the generated files in
`supplementary_figures/`.

For the Figure A.2 caption:

- use Delta rather than Delta-hat 4;
- state that overlap and kernel-test results are in Table A.1, not Table A.2;
- do not mark the porcupine full-distribution comparison significant
  (kernel p = 0.109).

Appendix B, Table B.1, and the Fisher-Freeman-Halton result are unchanged:
p = 0.781 (exact output 0.781252).

## Methods and table-caption consistency

These are not changes to biological results but are required for consistency
with the corrected analysis:

- State that Camelot timestamps were already camera-local
  (Asia/Jerusalem) before conversion to solar time; no UTC+2/UTC+3 offset was
  added.
- Remove the claim that the displayed activity curves contain bootstrap
  envelopes unless uncertainty ribbons are restored.
- Do not describe every Table 1 uncertainty value as a bootstrap SE.
  Delta-hat 4 rows use the standard deviation of the `compareCkern`
  randomization-null distribution; Delta-hat 1 rows use within-group
  bootstrap SEs.
- State that standard `overlapEst()` coefficients use independent detection
  events without group-size weights because that function does not support
  weights. The weighted curve displays and weighted Delta-hat 4
  `compareCkern()` tests retain `sighting_quantity` as the analysis weight.
- For nocturnality confidence intervals, say that independent detection
  events were resampled uniformly and `sighting_quantity` was retained as the
  analysis weight. Events were not also sampled in proportion to group size.

## Conclusions that remain supported

- Community nocturnality increases toward the high-culling end of the
  land-use gradient.
- Wolves, jackals, and wild boar show the clearest movement toward
  nocturnality under greater disturbance.
- Gazelles are least nocturnal in High culling and most nocturnal in
  Minefields.
- Gazelles become more nocturnal at high-human-activity cameras.
- Non-gazelle temporal overlap increases with disturbance, whereas
  gazelle-other overlap decreases.
- Temporal zoning retains the gazelle evening activity peak while human
  activity declines before sunset.
