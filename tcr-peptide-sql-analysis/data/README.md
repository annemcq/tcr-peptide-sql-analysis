# Data provenance

This repository intentionally does **not** re-host the source TCRen datasets,
contact maps, or other scientific inputs. Follow the source project's
license and original authors' conditions before redistributing data.

Expected inputs from `annemcq/tcr-peptide-ranking` (default branch):

- `processed/first_dataset_pairs.csv` — original synthetic negative candidates
- `processed/matched_composition_pairs.csv` — composition-matched candidate set
- `data/contact_maps_PDB.csv` — source contact records
- `data/TCRen_potential.csv` — residue-pair potential lookup
- `processed/data_tcren_score.csv` — optional original Python score reference
- `processed/data_tcren_features_matched_composition.csv` — optional matched Python score reference

The builder interprets `residue.index.to` as **0-based**, consistent with the
source `compute_tcren_score` Python implementation.

Dataset facts are described in the source repository's README and may change
in future commits. Always inspect the inputs before computing scientific results.
