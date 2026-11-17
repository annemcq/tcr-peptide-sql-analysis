-- 02.1 Contact count per TCR (source contact records, not unique atom pairs).
SELECT pdb_id, COUNT(*) AS n_contacts
FROM native_contacts
GROUP BY pdb_id
ORDER BY n_contacts DESC, pdb_id;

-- 02.2 Highest-scoring residue contributions (mean signed potential per type).
SELECT nc.tcr_aa, r.amino_acid AS peptide_aa,
       COUNT(*) AS n_contributions,
       ROUND(AVG(p.interaction_energy), 5) AS mean_potential
FROM candidates c
JOIN native_contacts nc ON nc.pdb_id = c.pdb_id
JOIN candidate_residues r
 ON r.candidate_id = c.candidate_id AND r.peptide_position = nc.peptide_position
JOIN tcren_potentials p ON p.tcr_aa = nc.tcr_aa AND p.peptide_aa = r.amino_acid
WHERE c.scenario = 'original'
GROUP BY nc.tcr_aa, r.amino_acid
ORDER BY n_contributions DESC, nc.tcr_aa, r.amino_acid
LIMIT 15;
