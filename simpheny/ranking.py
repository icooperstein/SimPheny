"""Aggregate supporting matches into candidate-gene rankings."""
import pandas as pd


def add_gene_rankings(matches: pd.DataFrame):
    """Add SimPheny Gene Scores and candidate ranks independently per query."""
    gene_columns = [
        "Query_ID",
        "Gene_Hit",
        "SimPheny_Gene_Score",
        "Candidate_Rank",
        "Num_Reference_Matches",
        "Best_SimPheny_Score",
        "Best_Reference_Source",
        "Best_Reference_ID",
        "Sources",
    ]

    if matches.empty:
        return matches.copy(), pd.DataFrame(columns=gene_columns)

    gene_summaries = []

    for query_id, query_matches in matches.groupby("Query_ID", sort=False):
        query_rows = []

        for gene, group in query_matches.groupby("Gene_Hit", sort=False):
            ordered = group.sort_values("SimPheny_Score", ascending=False)
            top_scores = ordered["SimPheny_Score"].head(2).tolist()
            best = ordered.iloc[0]

            query_rows.append(
                {
                    "Query_ID": query_id,
                    "Gene_Hit": gene,
                    "SimPheny_Gene_Score": sum(top_scores) / len(top_scores),
                    "Num_Reference_Matches": len(group),
                    "Best_SimPheny_Score": float(best["SimPheny_Score"]),
                    "Best_Reference_Source": best["Reference_Source"],
                    "Best_Reference_ID": best["Reference_ID"],
                    "Sources": ";".join(sorted(set(group["Reference_Source"]))),
                }
            )

        query_genes = (
            pd.DataFrame(query_rows)
            .sort_values(
                ["SimPheny_Gene_Score", "Gene_Hit"],
                ascending=[False, True],
            )
            .reset_index(drop=True)
        )
        query_genes["Candidate_Rank"] = query_genes.index + 1
        gene_summaries.append(query_genes)

    genes = pd.concat(gene_summaries, ignore_index=True)

    matches_out = matches.merge(
        genes[["Query_ID", "Gene_Hit", "SimPheny_Gene_Score", "Candidate_Rank"]],
        on=["Query_ID", "Gene_Hit"],
        how="left",
    )

    matches_out = matches_out.sort_values(
        ["Query_ID", "Candidate_Rank", "SimPheny_Score"],
        ascending=[True, True, False],
    ).reset_index(drop=True)

    genes = genes[
        [
            "Query_ID",
            "Candidate_Rank",
            "Gene_Hit",
            "SimPheny_Gene_Score",
            "Num_Reference_Matches",
            "Best_SimPheny_Score",
            "Best_Reference_Source",
            "Best_Reference_ID",
            "Sources",
        ]
    ]

    return matches_out, genes
