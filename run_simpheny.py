#!/usr/bin/env python3
"""Command-line interface for SimPheny."""
import argparse
import json
from datetime import datetime
from pathlib import Path

from simpheny.io import read_queries
from simpheny.pipeline import run_simpheny


def parser():
    p = argparse.ArgumentParser(
        description=(
            "Prioritize candidate genes for one or more undiagnosed patients "
            "using SimPheny patient-to-patient phenotypic matching."
        )
    )
    p.add_argument(
        "--input",
        required=True,
        help=(
            "Query TSV with one patient per row and columns "
            "ID, HPO_Terms, and Candidate_Genes."
        ),
    )
    p.add_argument(
        "--reference",
        nargs="+",
        required=True,
        help="Choose: udn, clinvar, phenopacket, decipher, orphanet, or all.",
    )
    p.add_argument(
        "--output-dir",
        default=None,
        help=(
            "Optional output directory. If omitted, SimPheny creates a unique "
            "timestamped directory under results/."
        ),
    )
    p.add_argument("--iterations", type=int, default=10000)
    p.add_argument("--seed", type=int, default=None)
    p.add_argument(
        "--show-restricted-reference-ids",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    return p


def _safe_name(value):
    value = str(value).strip()
    safe = "".join(
        c if c.isalnum() or c in ("-", "_") else "_"
        for c in value
    )
    return safe.strip("_") or "query"


def _default_output_dir(input_file, references):
    queries = read_queries(input_file)

    if len(queries) == 1:
        query_label = _safe_name(queries[0]["id"])
    else:
        query_label = f"{_safe_name(Path(input_file).stem)}_{len(queries)}patients"

    refs = [r.lower() for r in references]
    if "all" in refs:
        ref_label = "all"
    else:
        ref_label = "-".join(_safe_name(r) for r in refs)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return Path("results") / f"{query_label}_{ref_label}_{timestamp}"


def main():
    args = parser().parse_args()

    if args.output_dir:
        outdir = Path(args.output_dir)
    else:
        outdir = _default_output_dir(args.input, args.reference)

    outdir.mkdir(parents=True, exist_ok=False)

    matches, genes, metadata = run_simpheny(
        query_file=args.input,
        references=args.reference,
        iterations=args.iterations,
        seed=args.seed,
        mask_restricted_ids=not args.show_restricted_reference_ids,
    )

    matches_path = outdir / "simpheny_matches.tsv"
    genes_path = outdir / "simpheny_genes.tsv"
    metadata_path = outdir / "run_metadata.json"

    matches.to_csv(matches_path, sep="\t", index=False)
    genes.to_csv(genes_path, sep="\t", index=False)

    metadata["output_directory"] = str(outdir)
    metadata_path.write_text(json.dumps(metadata, indent=2))

    print("\nSimPheny complete.")
    print(f"Query patients: {metadata['num_query_patients']}")
    print("References:", ", ".join(metadata["references"]))
    print("Total matches:", metadata["total_matches"])

    if not genes.empty:
        print("\nTop candidate genes:")

        for query_id in metadata["query_ids"]:
            query_genes = genes[genes["Query_ID"] == query_id].head(10)
            print(f"\n{query_id}")
            if query_genes.empty:
                print("  No candidate-gene matches found.")
            else:
                print(
                    query_genes[
                        [
                            "Candidate_Rank",
                            "Gene_Hit",
                            "SimPheny_Gene_Score",
                            "Num_Reference_Matches",
                            "Sources",
                        ]
                    ].to_string(index=False)
                )

    print("\nResults written to:")
    print(f"  {outdir}")
    print(f"  - {genes_path.name}")
    print(f"  - {matches_path.name}")
    print(f"  - {metadata_path.name}")


if __name__ == "__main__":
    main()
