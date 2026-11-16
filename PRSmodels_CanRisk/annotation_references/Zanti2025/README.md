# Published BRCA1 and BRCA2 variant odds ratios

`BRCA1_BRCA2_variant_odds_ratios.xlsx` is an unchanged copy of Supplementary
Data 1–8 from Zanti et al. (2025), *Analysis of more than 400,000 women provides
case-control evidence for BRCA1 and BRCA2 variant classification*.

- Article: https://doi.org/10.1038/s41467-025-59979-6
- Original workbook: https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-025-59979-6/MediaObjects/41467_2025_59979_MOESM4_ESM.xlsx
- Downloaded: 2026-11-12

Use the **Supplementary Data 4** sheet. The lowest header row is Excel row 5.
Filter column J (`GENE`) to `BRCA1`. There are 4,291 BRCA1 rows, of which 1,095
have a reported combined-dataset odds ratio in column AO. This includes variants
of multiple classifications, not only pathogenic variants.

Key columns:

| Excel column | Meaning |
|---|---|
| A, B | Variant identifiers in hg19 and hg38, respectively |
| C–G | Chromosome, positions in both assemblies, REF and ALT |
| M, N | HGVS cDNA and protein annotations |
| P, Q | ClinVar review status and classification used by the publication |
| AC, AF, AI | Study-specific OR and 95% CI for BRIDGES, CARRIERS and UKB |
| AO, AP | Combined-dataset OR and 95% CI, and p-value |
| AL | Combined case-control likelihood ratio for pathogenicity (not an OR) |

Keep missing ORs missing. The likelihood-ratio columns assess evidence for
pathogenicity and must not be substituted for disease odds ratios. Some ORs have
wide confidence intervals because variants are rare. The publication's ClinVar
classifications are a historical snapshot.

This workbook is an association-study reference, not a CanRisk coefficient file
or a validated PRS scoring file. CanRisk accepts gene-level genetic test results
and models risk using its own parameters. Its official file specification
explicitly lists BRCA1 as pedigree field 18:
https://canrisk.atlassian.net/wiki/spaces/FAQS/pages/2749661185/CanRisk+Pedigree+Data+v4+File+Format+Specification
