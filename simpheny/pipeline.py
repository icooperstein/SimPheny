"""End-to-end SimPheny workflow for one or more query patients."""
from typing import List, Optional, Tuple
import time
import pandas as pd

from .config import REFERENCE_CONFIG, GENE_BACKGROUND_FILE, normalize_reference_names
from .io import (
    read_queries,
    read_reference,
    read_gene_background,
    read_phenotype_background,
)
from .similarity import initialize_ontology, validate_hpo_terms
from .matching import identify_matches
from .empirical import empirical_p_values
from .scoring import combine_fixed_ebm, simpheny_score, confidence_tier
from .privacy import mask_reference_ids
from .ranking import add_gene_rankings


def _format_duration(seconds):
    seconds = max(0, int(round(seconds)))
    minutes, seconds = divmod(seconds, 60)
    hours, minutes = divmod(minutes, 60)
    if hours:
        return f"{hours}h {minutes}m {seconds}s"
    if minutes:
        return f"{minutes}m {seconds}s"
    return f"{seconds}s"


def _prepare_reference(source):
    cfg = REFERENCE_CONFIG[source]
    label = cfg["label"]

    print(f"\n[{label}] Loading reference data...", flush=True)
    reference = read_reference(cfg["reference_file"])
    print(f"[{label}] Loaded {len(reference):,} reference records.", flush=True)

    print(f"[{label}] Validating reference HPO terms...", flush=True)
    cleaned = []
    invalid_count = 0

    for terms in reference["HPO_Terms_List"]:
        valid, invalid = validate_hpo_terms(terms)
        cleaned.append(valid)
        invalid_count += len(invalid)

    reference["HPO_Terms_List"] = cleaned
    reference = reference[
        reference["HPO_Terms_List"].map(bool)
    ].reset_index(drop=True)

    if invalid_count:
        print(
            f"[{label}] Ignored {invalid_count:,} invalid reference HPO annotations.",
            flush=True,
        )

    phenotype_background_file = cfg.get("phenotype_background_file")

    if phenotype_background_file is not None:
        print(f"[{label}] Loading phenotype background...", flush=True)
        phenotype_corpus = read_phenotype_background(phenotype_background_file)
        phenotype_corpus, invalid_background = validate_hpo_terms(phenotype_corpus)

        if invalid_background:
            print(
                f"[{label}] Ignored {len(invalid_background):,} invalid "
                "phenotype-background annotations.",
                flush=True,
            )
    else:
        print(
            f"[{label}] Building phenotype background from reference cohort...",
            flush=True,
        )
        phenotype_corpus = [
            term
            for terms in reference["HPO_Terms_List"]
            for term in terms
        ]

    print(
        f"[{label}] Phenotype background contains "
        f"{len(phenotype_corpus):,} HPO annotations.",
        flush=True,
    )

    return reference, phenotype_corpus


def _score_query_against_reference(
    query,
    source,
    reference,
    phenotype_corpus,
    gene_corpus,
    iterations,
    seed,
    mask_ids=True,
):
    cfg = REFERENCE_CONFIG[source]
    label = cfg["label"]

    print(
        f"[{label}] Calculating phenotype similarities and gene overlaps "
        f"for {query['id']}...",
        flush=True,
    )

    matches = identify_matches(query, reference, source=source).reset_index(drop=True)
    print(
        f"[{label}] {query['id']}: found {len(matches):,} candidate-gene matches.",
        flush=True,
    )

    if matches.empty:
        return matches

    ref_terms = dict(zip(reference["ID"], reference["HPO_Terms_List"]))
    scored = []
    scoring_start = time.perf_counter()
    total = len(matches)

    print(
        f"[{label}] {query['id']}: calculating empirical p-values "
        f"({iterations:,} iterations per match)...",
        flush=True,
    )

    for match_index in range(total):
        row = matches.iloc[match_index]
        display_index = match_index + 1
        match_start = time.perf_counter()

        match_seed = None if seed is None else seed + match_index

        pheno_p, gene_p = empirical_p_values(
            reference_terms=ref_terms[str(row["Reference_ID"])],
            observed_similarity=float(row["PhenoSimJaccard"]),
            hit_gene=str(row["Gene_Hit"]),
            query_hpo_count=len(query["hpo_terms"]),
            query_candidate_count=len(query["candidate_genes"]),
            phenotype_term_corpus=phenotype_corpus,
            gene_corpus=gene_corpus,
            iterations=iterations,
            seed=match_seed,
        )

        ebm = combine_fixed_ebm(pheno_p, gene_p, source)
        score = simpheny_score(ebm)

        out = row.to_dict()
        out.update(
            {
                "pheno_p": pheno_p,
                "gene_p": gene_p,
                "EBM_p": ebm,
                "SimPheny_Score": score,
                "Confidence": confidence_tier(score),
            }
        )
        scored.append(out)

        elapsed = time.perf_counter() - scoring_start
        avg_per_match = elapsed / display_index
        remaining = avg_per_match * (total - display_index)
        match_elapsed = time.perf_counter() - match_start

        print(
            f"[{label}] {query['id']} [{display_index}/{total}] "
            f"{row['Gene_Hit']} scored in {_format_duration(match_elapsed)} | "
            f"elapsed {_format_duration(elapsed)} | "
            f"ETA {_format_duration(remaining)}",
            flush=True,
        )

    df = pd.DataFrame(scored)

    if mask_ids and cfg["mask_reference_ids"]:
        df = mask_reference_ids(df, source)

    print(f"[{label}] {query['id']}: scoring complete.", flush=True)
    return df


def run_simpheny(
    query_file: str,
    references: List[str],
    iterations: int = 10_000,
    seed: Optional[int] = None,
    mask_restricted_ids: bool = True,
) -> Tuple[pd.DataFrame, pd.DataFrame, dict]:
    total_start = time.perf_counter()

    print("Initializing ontology...", flush=True)
    initialize_ontology()

    print(f"Loading query patients: {query_file}", flush=True)
    queries = read_queries(query_file)

    for query in queries:
        valid, invalid = validate_hpo_terms(query["hpo_terms"])
        if invalid:
            raise ValueError(
                f"Query {query['id']} contains invalid HPO term(s): "
                + ", ".join(invalid)
            )
        query["hpo_terms"] = valid

    print(f"Loaded {len(queries):,} query patient(s).", flush=True)
    for query in queries:
        print(
            f"  {query['id']}: {len(query['hpo_terms'])} HPO terms, "
            f"{len(query['candidate_genes'])} candidate genes.",
            flush=True,
        )

    selected = normalize_reference_names(references)
    print("Reference dataset(s): " + ", ".join(selected), flush=True)

    print("Loading empirical gene background...", flush=True)
    gene_corpus = read_gene_background(str(GENE_BACKGROUND_FILE))
    print(
        f"Gene background contains {len(gene_corpus):,} gene observations.",
        flush=True,
    )

    # Load/validate each reference only once, even for multi-patient input.
    prepared_references = {}
    for source in selected:
        prepared_references[source] = _prepare_reference(source)

    all_matches = []
    matches_per_query = {}
    matches_per_reference = {source: 0 for source in selected}

    for query_number, query in enumerate(queries, start=1):
        print(
            f"\n=== Query {query_number}/{len(queries)}: {query['id']} ===",
            flush=True,
        )

        query_match_count = 0

        for source in selected:
            reference, phenotype_corpus = prepared_references[source]

            # Offset seed by query number so distinct patients do not receive
            # identical Monte Carlo streams when --seed is supplied.
            query_seed = None
            if seed is not None:
                query_seed = seed + ((query_number - 1) * 1_000_000)

            df = _score_query_against_reference(
                query=query,
                source=source,
                reference=reference,
                phenotype_corpus=phenotype_corpus,
                gene_corpus=gene_corpus,
                iterations=iterations,
                seed=query_seed,
                mask_ids=mask_restricted_ids,
            )

            query_match_count += len(df)
            matches_per_reference[source] += len(df)

            if not df.empty:
                all_matches.append(df)

        matches_per_query[query["id"]] = query_match_count

    matches = (
        pd.concat(all_matches, ignore_index=True)
        if all_matches
        else pd.DataFrame()
    )

    print("\nRanking candidate genes independently for each query patient...", flush=True)
    matches, genes = add_gene_rankings(matches)

    total_elapsed = time.perf_counter() - total_start

    metadata = {
        "query_file": query_file,
        "query_ids": [q["id"] for q in queries],
        "num_query_patients": len(queries),
        "references": selected,
        "iterations": iterations,
        "matches_per_query": matches_per_query,
        "matches_per_reference": matches_per_reference,
        "total_matches": len(matches),
        "total_ranked_query_gene_pairs": len(genes),
        "restricted_reference_ids_masked": mask_restricted_ids,
        "elapsed_seconds": total_elapsed,
    }

    print(f"SimPheny analysis finished in {_format_duration(total_elapsed)}.", flush=True)
    return matches, genes, metadata
