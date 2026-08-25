"""PhenoSimJaccard implementation."""
from typing import Iterable, List, Tuple
import pyhpo
from pyhpo import Ontology, HPOSet
from pyhpo.similarity.base import SimScore, SimilarityBase
from .config import DEFAULT_METHOD, DEFAULT_KIND, DEFAULT_COMBINE

class CustomJaccardIC(SimilarityBase):
    def __call__(self, term1: "pyhpo.HPOTerm", term2: "pyhpo.HPOTerm",
                 kind: str, dependencies: List[float]) -> float:
        if term1 == term2:
            return 1.0
        common = sum(x.information_content[kind] for x in term1.common_ancestors(term2))
        union = sum(
            x.information_content[kind]
            for x in (term1.all_parents | term2.all_parents)
        )
        if term1 in term2.all_parents and term2 in term1.all_parents:
            pass
        elif term1 in term2.all_parents:
            union += term2.information_content[kind]
        elif term2 in term1.all_parents:
            union += term1.information_content[kind]
        else:
            union += (
                term1.information_content[kind]
                + term2.information_content[kind]
            )
        return common / union if union else 0.0

_registered = False

def initialize_ontology():
    global _registered
    _ = Ontology()
    if not _registered:
        try:
            SimScore.register(DEFAULT_METHOD, CustomJaccardIC)
        except Exception:
            pass
        _registered = True

def validate_hpo_terms(terms: Iterable[str]) -> Tuple[List[str], List[str]]:
    valid, invalid = [], []
    for term in terms:
        term = str(term).strip()
        if not term:
            continue
        try:
            Ontology.get_hpo_object(term)
            valid.append(term)
        except Exception:
            invalid.append(term)
    return valid, invalid

def hpo_set(terms: Iterable[str]) -> HPOSet:
    valid, invalid = validate_hpo_terms(terms)
    if invalid:
        raise ValueError("Invalid HPO term(s): " + ", ".join(invalid))
    if not valid:
        raise ValueError("No valid HPO terms were provided.")
    return HPOSet.from_queries(valid)

def phenosim_jaccard(terms1: Iterable[str], terms2: Iterable[str]) -> float:
    a, b = hpo_set(terms1), hpo_set(terms2)
    return float(a.similarity(
        b, method=DEFAULT_METHOD, kind=DEFAULT_KIND, combine=DEFAULT_COMBINE
    ))
