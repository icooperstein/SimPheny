"""Reference-specific EBM scoring."""
from math import log, log10
from scipy.stats import chi2
from .config import (
    EBM_PARAMETERS, HIGH_CONFIDENCE_THRESHOLD, MEDIUM_CONFIDENCE_THRESHOLD
)

def combine_fixed_ebm(pheno_p: float, gene_p: float, source: str) -> float:
    p = EBM_PARAMETERS[source]
    stat = -2.0 * (log(pheno_p) + log(gene_p))
    adjusted_stat = stat / p["scale"]
    combined = float(chi2.sf(adjusted_stat, p["dof"]))
    return max(combined, 1e-300)

def simpheny_score(p: float) -> float:
    return -log10(p)

def confidence_tier(score: float) -> str:
    if score >= HIGH_CONFIDENCE_THRESHOLD:
        return "High"
    if score >= MEDIUM_CONFIDENCE_THRESHOLD:
        return "Medium"
    return "Low"
