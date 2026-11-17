# Publish this project on GitHub

1. Sign in to https://github.com/new.
2. Repository name: `tcr-peptide-sql-analysis`.
3. Description: `Reproducible SQLite analysis of structure-based TCR–peptide scoring, tie-aware ranking, and matched-negative controls.`
4. Visibility: Public. **Do not tick** Add README, .gitignore, or license when creating the empty repository; those items are already provided or can be added deliberately later.
5. Create repository.
6. Unzip the downloaded package on your computer. Open the extracted folder. On your empty GitHub repository click **uploading an existing file** (or **Add file → Upload files**) and drag in the **contents** of the extracted folder, preserving `src/`, `sql/`, `tests/`, `.github/workflows/`, etc. If GitHub's browser uploader does not preserve the folder hierarchy, use the Git command-line method below instead.
7. Commit message: `Initial reproducible SQL biophysics analysis` → Commit changes.
8. Check **Actions → SQL validation**. It should run on your push (tests do not need the real source data).
9. For local reproduction, follow README.md. Actual full-data analysis requires downloading the public source CSVs first and is not included in this archive.

## Reliable alternative: Git terminal commands

Extract this package and open a terminal **inside the extracted `tcr-peptide-sql-analysis` directory**:

```bash
git init
git add .
git commit -m "Initial reproducible SQL biophysics analysis"
git branch -M main
git remote add origin https://github.com/annemcq/tcr-peptide-sql-analysis.git
git push -u origin main
```

If Git requests authentication, sign in using Git Credential Manager, `gh auth login`, or an approved personal access token. Do not paste tokens into chats or source files.

## Running it yourself

```bash
python -m pip install -r requirements.txt
python -m pytest -q
python -m src.fetch_public_data
python -m src.build_database --original data/raw/original.csv --matched data/raw/matched.csv --contacts data/raw/contacts.csv --potential data/raw/potential.csv --reference-original data/raw/reference_original.csv
python -m src.run_analysis
```

Optional figures: `python -m pip install numpy matplotlib && python -m src.plot_results`.

### Important
- The public source dataset is fetched at run time, **not redistributed** in this project.
- The database builder refuses to overwrite `results/tcr_analysis.sqlite`. Delete that file intentionally before rebuilding.
- This project complements the existing TCR repository and does not claim a new ML model, new biological measurements or experimentally validated negative peptides.
- Review and be able to explain every SQL query before including it in job applications.
