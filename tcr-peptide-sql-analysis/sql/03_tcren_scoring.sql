-- 03.1 Molecular scoring from contact-level joins and SUM (see view DDL).
SELECT scenario, pdb_id, candidate_id, label, n_contacts, sql_score
FROM candidate_scores
ORDER BY scenario, pdb_id, sql_score, candidate_id;

-- 03.2 Validate SQL scores against the Python reference, when supplied.
SELECT scenario,
       COUNT(python_score) AS n_reference_scores,
       ROUND(MAX(ABS(sql_score - python_score)), 10) AS max_abs_error,
       ROUND(AVG(ABS(sql_score - python_score)), 10) AS mean_abs_error
FROM candidate_scores
GROUP BY scenario
ORDER BY scenario;

-- 03.3 Any discrepancy exceeding tolerance (floating-point summation may differ).
SELECT scenario, pdb_id, peptide, sql_score, python_score,
       ABS(sql_score - python_score) AS abs_error
FROM candidate_scores
WHERE python_score IS NOT NULL
  AND ABS(sql_score - python_score) > 1e-7
ORDER BY abs_error DESC;
