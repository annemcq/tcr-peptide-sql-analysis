PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS tcr_systems (
    pdb_id TEXT PRIMARY KEY
);

CREATE TABLE IF NOT EXISTS scenarios (
    scenario TEXT PRIMARY KEY,
    description TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS candidates (
    candidate_id INTEGER PRIMARY KEY,
    scenario TEXT NOT NULL REFERENCES scenarios(scenario),
    pdb_id TEXT NOT NULL REFERENCES tcr_systems(pdb_id),
    peptide TEXT NOT NULL,
    label INTEGER NOT NULL CHECK (label IN (0, 1)),
    python_score REAL,
    UNIQUE (scenario, pdb_id, peptide)
);
CREATE INDEX IF NOT EXISTS idx_candidate_group ON candidates(scenario, pdb_id);

CREATE TABLE IF NOT EXISTS candidate_residues (
    candidate_id INTEGER NOT NULL REFERENCES candidates(candidate_id) ON DELETE CASCADE,
    peptide_position INTEGER NOT NULL CHECK (peptide_position >= 0),
    amino_acid TEXT NOT NULL CHECK (length(amino_acid) = 1),
    PRIMARY KEY (candidate_id, peptide_position)
);

CREATE TABLE IF NOT EXISTS native_contacts (
    contact_id INTEGER PRIMARY KEY,
    pdb_id TEXT NOT NULL REFERENCES tcr_systems(pdb_id),
    peptide_position INTEGER NOT NULL CHECK (peptide_position >= 0),
    tcr_aa TEXT NOT NULL CHECK (length(tcr_aa) = 1)
);
CREATE INDEX IF NOT EXISTS idx_contacts_pdb ON native_contacts(pdb_id, peptide_position);

CREATE TABLE IF NOT EXISTS tcren_potentials (
    tcr_aa TEXT NOT NULL CHECK (length(tcr_aa) = 1),
    peptide_aa TEXT NOT NULL CHECK (length(peptide_aa) = 1),
    interaction_energy REAL NOT NULL,
    PRIMARY KEY (tcr_aa, peptide_aa)
);

-- Scores can be recomputed from raw amino-acid contacts without ML inference.
CREATE VIEW IF NOT EXISTS candidate_scores AS
SELECT c.candidate_id, c.scenario, c.pdb_id, c.peptide, c.label,
       c.python_score, COUNT(*) AS n_contacts,
       SUM(p.interaction_energy) AS sql_score
FROM candidates AS c
JOIN native_contacts AS nc ON nc.pdb_id = c.pdb_id
JOIN candidate_residues AS r
  ON r.candidate_id = c.candidate_id
 AND r.peptide_position = nc.peptide_position
JOIN tcren_potentials AS p
  ON p.tcr_aa = nc.tcr_aa
 AND p.peptide_aa = r.amino_acid
GROUP BY c.candidate_id;
