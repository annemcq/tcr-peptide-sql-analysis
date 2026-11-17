# TCR–Peptide Recognition Through Relational SQL

**SQLite · molecular contact maps · structure-based scoring · scientific validation**

![SQL validation](https://github.com/annemcq/tcr-peptide-sql-analysis/actions/workflows/tests.yml/badge.svg)

This project asks whether a structure-based T-cell receptor (TCR) peptide score
can be reconstructed in SQL from its component residue contacts and statistical
interaction potentials. It then examines how the **cognate peptide's rank**
changes between two benchmark designs: uniformly generated negative peptides
and composition-matched shuffled negatives.

The project is a *relational analysis and reproducibility exercise*, not a new
TCR-recognition model. It complements the source repository
[`annemcq/tcr-peptide-ranking`](https://github.com/annemcq/tcr-peptide-ranking)
by expressing its core contact-scoring logic as auditable SQL joins.

## Research motivation

A TCR–peptide scoring workflow draws on several kinds of information: the TCR
system, its residue contact map, the candidate peptide sequence, and a lookup
potential for the interacting amino-acid pair. In a flat table, the connection
between these data types is obscured. A relational database makes the joins
explicit and enables questions such as:

1. Can we reproduce every original Python TCRen score, without retraining any ML model?
2. Which contact pairs contribute to the scores?
3. How sensitive are *within-TCR* positive peptide ranks to the design of negative candidates?
4. Which records would disappear silently from an inner join, and how do we catch them?

The source benchmark contains **25 TCR systems**, each with one cognate and 100
synthetic negative peptides per scenario. These synthetic negatives are not
experimentally validated nonbinders. All published numeric results in the source
repository are properties of that benchmark, not validated binding probabilities.

## Database model

| Table/view | Granularity | Key fields |
| --- | --- | --- |
| `tcr_systems` | TCR system | `pdb_id` |
| `scenarios` | Negative-sampling design | `scenario` |
| `candidates` | One peptide in one TCR and scenario | `candidate_id`, `label`, `peptide` |
| `candidate_residues` | One position of a candidate peptide | `(candidate_id, peptide_position)` |
| `native_contacts` | One TCR–peptide contact record | `contact_id`, `pdb_id`, `tcr_aa`, `peptide_position` |
| `tcren_potentials` | One amino-acid interaction pair | `(tcr_aa, peptide_aa)` |
| `candidate_scores` | SQL view with one reconstructed score per candidate | `sql_score`, `n_contacts` |

The `peptide_position` uses **zero-based indexing**, matching the source Python
scorer's `peptide[pos]`. Input sequences are normalized to uppercase and system
IDs to lowercase. Duplicated candidates, invalid residue positions, missing
interaction potentials, and unmatched candidates cause explicit validation
failures rather than misleading scores.

### Scoring in SQL

The central calculation is:

```sql
SELECT c.candidate_id,
       SUM(p.interaction_energy) AS sql_score
FROM candidates c
JOIN native_contacts nc ON nc.pdb_id = c.pdb_id
JOIN candidate_residues r
  ON r.candidate_id = c.candidate_id
 AND r.peptide_position = nc.peptide_position
JOIN tcren_potentials p
  ON p.tcr_aa = nc.tcr_aa AND p.peptide_aa = r.amino_acid
GROUP BY c.candidate_id;
```

Each source contact record contributes one interaction potential to the sum.
The stored source contact geometry is **reused for every candidate** belonging
to that TCR; alternative peptide binding structures are not modeled here.
The resulting number is a *statistical potential score*, **not** an absolute
thermodynamic binding free energy.

### Tie-aware ranking

`RANK()` does not reproduce the source project's `pandas.rank(method='average')`.
For an equal-score group occupying ranks `k` through `k + m - 1`, use
`k + (m - 1)/2.0`. The ranking query explicitly computes this correction with
window functions. It then calculates each cognate peptide's reciprocal rank,
Top-1 and Top-5. All comparisons are evaluated **within each TCR**.

## Ten reproducible SQL analyses

| SQL file | Analysis questions |
| --- | --- |
| `sql/01_data_quality.sql` | 1. How many candidate peptides and cognates are there per system? 2. Which contact/potential lookups are missing? 3. Are any candidates silently lost by joins? Plus a consistency check of cognate scores across scenarios. |
| `sql/02_contact_statistics.sql` | 4. How does native contact count vary by TCR? 5. Which residue-pair interactions contribute most often? |
| `sql/03_tcren_scoring.sql` | 6. Can SQL reconstruct the full scores? 7. How accurately do they match Python? 8. Are any scores outside a numerical tolerance? |
| `sql/04_peptide_ranking.sql` | 9. Where does the cognate rank in each TCR, accounting for ties? |
| `sql/05_negative_controls.sql` | 10. How do ranks and aggregate MRR/Top-k change when negative candidates are composition-matched? |

These are **10 scientific analysis questions**, implemented as 12 SQL statements,
not 10 arbitrary database examples. Methods include relational joins, `GROUP BY`,
`LEFT JOIN` validation, `CASE`, CTEs, window functions and conditional aggregation.

## Reproduce locally

Requires Python 3.10+ and SQLite (included in Python's standard library).
`pytest` is the only required third-party dependency for tests; `matplotlib`
plus `numpy` are optional for figures.

```bash
python -m pip install -r requirements.txt
python -m pytest -q
```

Option A: fetch the five public CSV inputs at runtime (the source data are **not** redistributed here), then build the database. For strict reproducibility, pass a source commit SHA with `--ref` rather than relying on the moving `main` branch:

```bash
python -m src.fetch_public_data
python -m src.build_database \
  --original data/raw/original.csv \
  --matched data/raw/matched.csv \
  --contacts data/raw/contacts.csv \
  --potential data/raw/potential.csv \
  --reference-original data/raw/reference_original.csv
python -m src.run_analysis
python -m pip install numpy matplotlib
python -m src.plot_results
```

Option B: obtain the original files from a local clone of the source repository (see [`data/README.md`](data/README.md)) and use their corresponding paths:

```bash
python -m src.build_database \
  --original ../tcr-peptide-ranking/processed/first_dataset_pairs.csv \
  --matched ../tcr-peptide-ranking/processed/matched_composition_pairs.csv \
  --contacts ../tcr-peptide-ranking/data/contact_maps_PDB.csv \
  --potential ../tcr-peptide-ranking/data/TCRen_potential.csv \
  --reference-original ../tcr-peptide-ranking/processed/data_tcren_score.csv \
  --reference-matched ../tcr-peptide-ranking/processed/data_tcren_features_matched_composition.csv

python -m src.run_analysis
python -m pip install numpy matplotlib
python -m src.plot_results
```

If the optional `--reference-*` files are unavailable, omit those flags;
scoring and rank analyses still run, but the Python-parity report cannot
be assessed for those rows. The builder *refuses to overwrite* an existing
SQLite database: delete `results/tcr_analysis.sqlite` deliberately before
rerunning the import. Output CSV reports appear under `results/reports/` and
figures under `figures/`.

### Independently checked upstream benchmark

An independent contact-loop calculation using the public GitHub inputs (before executing this project's full SQLite workflow) verified the following reference values on 17 November 2026:

| Upstream check | Observed value |
| --- | ---: |
| Candidate peptides per scenario | 2,525 |
| TCR systems per scenario | 25 |
| Original-score maximum absolute difference from public Python scores | 2.66 × 10⁻¹⁵ |
| TCRen full-dataset MRR, original negative set | 0.3349 |
| TCRen full-dataset MRR, composition-matched set | 0.2677 |

These are **independent checks of the original public inputs**, not claims that this repository's SQLite pipeline has already been run on the full dataset. They are also full-dataset baseline summaries; the source README reports **averages over repeated held-out grouped splits**, which need not equal these numbers.

### Expected sanity checks

- The public benchmark's 25 systems and 101 candidates per system should
  yield 2,525 candidate rows **per scenario** if the source CSVs are unchanged.
- Every `(scenario,pdb_id)` group should contain **exactly one** positive peptide.
- SQL reconstruction should produce one score per candidate, with all contact
  lookup records accounted for.
- With source reference scores provided, SQL and Python should agree within an
  appropriate floating-point tolerance. The default discrepancy report shows
  absolute differences larger than `1e-7`; precision of the inputs matters.

**Results policy:** This repository template has a tested small synthetic
fixture, but the full public CSV analysis must be run locally. Do not treat
these expected counts as newly measured results, and do not claim parity or
new rank summaries until they are verified from the real data.

## Figures (generated only from the real imported data)

`src/plot_results.py` can produce:

- `rank_sensitivity.png` — within-TCR cognate ranks for original versus matched negatives
- `sql_python_parity.png` — SQL versus Python scores (requires reference CSV)

## Tests and automation

Run `python -m pytest -q`. Tests use **small synthetic inputs with known
answers**, exercising positive labels, contact sums, rank ties, missing lookup
errors, duplicate candidates and report export. CI runs pytest on every push
and pull request. Tests demonstrate software correctness on fixtures; they do
not substitute for quantitative validation against the full public dataset.

## Limitations and scientific interpretation

- The contact map for the cognate complex is reused to score decoy peptides.
  Different candidate binding geometries are not considered.
- Synthetic negative peptides are not confirmed biological nonbinders. The
  original repository shows that ML ranking results are sensitive to the
  negative-generation strategy; this project only re-evaluates the **TCRen
  baseline**, without rerunning logistic regression or random forests.
- The two negative sets change the candidate competition; changes in positive
  rank do not imply that TCRen itself changes for the same cognate peptide.
- Reproducing the original scores establishes computational consistency, not
  biophysical validity or transfer to new TCRs.
- SQL grouping is exact for the supplied schema, while real deployments may
  require careful indexing, performance benchmarking and transaction handling.

## Provenance and credit

The scientific scoring method and source dataset are from the original TCRen
work (Karnaukhov et al., *Nature Computational Science*, 2024,
DOI: [10.1038/s43588-024-00653-0](https://doi.org/10.1038/s43588-024-00653-0)).
The dataset construction, Python scorer and the matched-negative experiment
come from [`annemcq/tcr-peptide-ranking`](https://github.com/annemcq/tcr-peptide-ranking).
The contribution of this repository is a **normalized SQLite representation,
independently expressed SQL scoring, tie-aware ranking, reproducibility checks,
and an automated comparison pipeline**. The work was developed with coding
assistance; what matters is being able to explain and modify the implementation.

## Technology

Python · SQLite · SQL (CTEs, window functions, joins) · pytest · GitHub Actions
