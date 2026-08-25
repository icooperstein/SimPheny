"""Input/output helpers."""
from pathlib import Path
from typing import List
import pandas as pd


def split_semicolon(value) -> List[str]:
    if pd.isna(value):
        return []
    return [x.strip() for x in str(value).split(";") if x.strip()]


def read_queries(path: str) -> List[dict]:
    """Read one or more query patients from a TSV file.

    Required columns:
      ID, HPO_Terms, Candidate_Genes

    Each row is treated as one independent query patient.
    """
    df = pd.read_csv(path, sep="\t")
    required = {"ID", "HPO_Terms", "Candidate_Genes"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Query file missing columns: {sorted(missing)}")

    if df.empty:
        raise ValueError("Query file contains no patients.")

    if df["ID"].isna().any():
        raise ValueError("Every query row must contain an ID.")

    if df["ID"].astype(str).duplicated().any():
        duplicates = sorted(
            set(df.loc[df["ID"].astype(str).duplicated(keep=False), "ID"].astype(str))
        )
        raise ValueError(
            "Query patient IDs must be unique. Duplicate ID(s): "
            + ", ".join(duplicates)
        )

    queries = []
    for _, row in df.iterrows():
        hpo_terms = split_semicolon(row["HPO_Terms"])
        candidate_genes = split_semicolon(row["Candidate_Genes"])

        if not hpo_terms:
            raise ValueError(f"Query {row['ID']} has no HPO terms.")
        if not candidate_genes:
            raise ValueError(f"Query {row['ID']} has no candidate genes.")

        queries.append(
            {
                "id": str(row["ID"]),
                "hpo_terms": hpo_terms,
                "candidate_genes": candidate_genes,
            }
        )

    return queries


def read_reference(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(
            f"Reference data not found: {path}\n"
            "See reference_data/README.md for setup instructions."
        )

    df = pd.read_csv(path, sep="\t")
    required = {"ID", "HPO_Terms", "Diagnostic_Genes"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Reference file missing columns: {sorted(missing)}")

    out = df.copy()
    out["ID"] = out["ID"].astype(str)
    out["HPO_Terms_List"] = out["HPO_Terms"].apply(split_semicolon)
    out["Diagnostic_Genes_List"] = out["Diagnostic_Genes"].apply(split_semicolon)

    return out[
        out["HPO_Terms_List"].map(bool)
        & out["Diagnostic_Genes_List"].map(bool)
    ].reset_index(drop=True)


def read_gene_background(path: str) -> List[str]:
    df = pd.read_csv(path, sep="\t")
    if "Gene" not in df.columns:
        raise ValueError("Gene-background file must contain a Gene column.")

    genes = [str(x).strip() for x in df["Gene"].dropna() if str(x).strip()]
    if not genes:
        raise ValueError("Gene-background file contains no genes.")

    return genes


def read_phenotype_background(path: Path) -> List[str]:
    """Read the broader UDN phenotype-background corpus.

    Accepted formats:
      1. HPO_Terms column with semicolon-delimited terms per row
      2. HPO_Term column with one term per row

    Duplicate annotations are retained intentionally.
    """
    if not path.exists():
        raise FileNotFoundError(
            f"Phenotype background not found: {path}\n"
            "UDN requires the broader UDN phenotype corpus used in the manuscript analysis."
        )

    df = pd.read_csv(path, sep="\t")

    if "HPO_Terms" in df.columns:
        corpus = []
        for value in df["HPO_Terms"].dropna():
            corpus.extend(split_semicolon(value))
        if corpus:
            return corpus

    if "HPO_Term" in df.columns:
        corpus = [
            str(x).strip()
            for x in df["HPO_Term"].dropna()
            if str(x).strip()
        ]
        if corpus:
            return corpus

    raise ValueError(
        "Phenotype-background file must contain either an HPO_Terms column "
        "or an HPO_Term column."
    )
