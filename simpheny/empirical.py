"""Empirical phenotype and gene p-values."""
import random
from typing import Iterable, List, Optional, Tuple
from pyhpo import HPOSet
from .config import (
    DEFAULT_COMBINE, DEFAULT_KIND, DEFAULT_METHOD,
    DEFAULT_MAX_RANDOM_HPO_TERMS
)

def _sample_unique_terms(corpus: List[str], n: int, rng: random.Random) -> List[str]:
    if n > len(set(corpus)):
        raise ValueError("Not enough unique HPO terms in phenotype background.")
    sample = rng.sample(corpus, n)
    while len(set(sample)) != n:
        sample += rng.sample(corpus, n - len(set(sample)))
    # preserve first occurrence order
    return list(dict.fromkeys(sample))

def empirical_p_values(
    reference_terms: Iterable[str],
    observed_similarity: float,
    hit_gene: str,
    query_hpo_count: int,
    query_candidate_count: int,
    phenotype_term_corpus: List[str],
    gene_corpus: List[str],
    iterations: int,
    seed: Optional[int] = None,
) -> Tuple[float, float]:
    if query_candidate_count > len(gene_corpus):
        raise ValueError("Candidate list is larger than gene background corpus.")
    rng = random.Random(seed)
    n_hpo = min(query_hpo_count, DEFAULT_MAX_RANDOM_HPO_TERMS)
    ref_hpo = HPOSet.from_queries(list(reference_terms))

    pheno_hits = 0
    gene_hits = 0
    for _ in range(iterations):
        random_terms = _sample_unique_terms(phenotype_term_corpus, n_hpo, rng)
        random_hpo = HPOSet.from_queries(random_terms)
        score = random_hpo.similarity(
            ref_hpo, method=DEFAULT_METHOD, kind=DEFAULT_KIND, combine=DEFAULT_COMBINE
        )
        if score >= observed_similarity:
            pheno_hits += 1

        sampled_genes = rng.sample(gene_corpus, query_candidate_count)
        if hit_gene in set(sampled_genes):
            gene_hits += 1

    return (
        (pheno_hits + 1) / (iterations + 1),
        (gene_hits + 1) / (iterations + 1),
    )
