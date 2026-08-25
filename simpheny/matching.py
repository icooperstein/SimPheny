"""Identify candidate-gene overlaps with reference diagnoses."""
from typing import Dict, List
import pandas as pd
from .similarity import phenosim_jaccard

def identify_matches(query: dict, reference_df: pd.DataFrame, source: str) -> pd.DataFrame:
    candidate_genes = set(query["candidate_genes"])
    rows: List[Dict] = []
    for _, ref in reference_df.iterrows():
        sim = phenosim_jaccard(query["hpo_terms"], ref["HPO_Terms_List"])
        overlaps = candidate_genes.intersection(ref["Diagnostic_Genes_List"])
        for gene in sorted(overlaps):
            rows.append({
                "Query_ID": query["id"],
                "Reference_Source": source,
                "Reference_ID": str(ref["ID"]),
                "Gene_Hit": gene,
                "PhenoSimJaccard": sim,
            })
    return pd.DataFrame(rows)
