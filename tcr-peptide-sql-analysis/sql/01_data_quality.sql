-- 01.1 One row per group: each TCR/scenario should contain exactly one cognate.
SELECT scenario, pdb_id, COUNT(*) AS n_candidates,
       SUM(label) AS n_positives,
       SUM(CASE WHEN label = 0 THEN 1 ELSE 0 END) AS n_negatives
FROM candidates
GROUP BY scenario, pdb_id
ORDER BY scenario, pdb_id;

-- 01.2 Unmatched contacts: missing peptide position or potential lookup.
-- Each contact should generate one scored contribution for every candidate.
SELECT c.scenario, c.pdb_id, c.candidate_id, nc.contact_id,
       nc.peptide_position, nc.tcr_aa, r.amino_acid
FROM candidates c
JOIN native_contacts nc ON nc.pdb_id = c.pdb_id
LEFT JOIN candidate_residues r
  ON r.candidate_id = c.candidate_id
 AND r.peptide_position = nc.peptide_position
LEFT JOIN tcren_potentials p
  ON p.tcr_aa = nc.tcr_aa AND p.peptide_aa = r.amino_acid
WHERE r.candidate_id IS NULL OR p.tcr_aa IS NULL;

-- 01.3 One row per peptide must appear in candidate_scores, or joins lost candidates.
SELECT c.candidate_id, c.scenario, c.pdb_id
FROM candidates c
LEFT JOIN candidate_scores s ON s.candidate_id = c.candidate_id
WHERE s.candidate_id IS NULL;

-- 01.4 The same cognate peptide should have the same TCRen score in both scenarios.
-- A nonzero delta indicates a data/import/scoring problem, not a sampling effect.
SELECT a.pdb_id, a.peptide,
       a.sql_score AS original_cognate_score,
       b.sql_score AS matched_cognate_score,
       ABS(a.sql_score - b.sql_score) AS abs_difference
FROM candidate_scores a
JOIN candidate_scores b
  ON a.pdb_id = b.pdb_id AND a.peptide = b.peptide
WHERE a.scenario = 'original' AND b.scenario = 'composition_matched'
  AND a.label = 1 AND b.label = 1
ORDER BY abs_difference DESC, a.pdb_id;
