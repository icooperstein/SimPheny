"""Reference registry and manuscript-derived constants."""
from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
REFERENCE_ROOT = PACKAGE_ROOT / "reference_data"
GENE_BACKGROUND_FILE = REFERENCE_ROOT / "gene_background.tsv"

EBM_PARAMETERS = {
    "udn": {"rho": 0.140654089465608, "scale": 1.070327044732804, "dof": 3.737175491999793},
    "clinvar": {"rho": 0.06694913859450671, "scale": 1.0334745692972533, "dof": 3.8704387305049393},
    "phenopacket": {"rho": 0.16390782006306615, "scale": 1.0819539100315332, "dof": 3.6970151527835613},
    "decipher": {"rho": 0.121718358199912691, "scale": 1.0608591790995634, "dof": 3.7705287174826747},
    "orphanet": {"rho": 0.12629563845713915, "scale": 1.0631478192285695, "dof": 3.7624118938629243},
}

REFERENCE_CONFIG = {
    "udn": {
        "label": "UDN",
        "reference_file": REFERENCE_ROOT / "udn" / "reference.tsv",
        # UDN is the one source whose phenotype null uses the broader UDN
        # phenotype corpus rather than only the diagnosed reference cohort.
        "phenotype_background_file": REFERENCE_ROOT / "udn" / "phenotype_background.tsv",
        "restricted": True,
        "mask_reference_ids": True,
    },
    "clinvar": {
        "label": "ClinVar",
        "reference_file": REFERENCE_ROOT / "clinvar" / "reference.tsv",
        "phenotype_background_file": None,
        "restricted": False,
        "mask_reference_ids": False,
    },
    "phenopacket": {
        "label": "Phenopacket Store",
        "reference_file": REFERENCE_ROOT / "phenopacket" / "reference.tsv",
        "phenotype_background_file": None,
        "restricted": False,
        "mask_reference_ids": False,
    },
    "decipher": {
        "label": "DECIPHER",
        "reference_file": REFERENCE_ROOT / "decipher" / "reference.tsv",
        "phenotype_background_file": None,
        "restricted": True,
        "mask_reference_ids": True,
    },
    "orphanet": {
        "label": "Orphanet",
        "reference_file": REFERENCE_ROOT / "orphanet" / "reference.tsv",
        "phenotype_background_file": None,
        "restricted": False,
        "mask_reference_ids": False,
    },
}

HIGH_CONFIDENCE_THRESHOLD = 4.5
MEDIUM_CONFIDENCE_THRESHOLD = 2.5
DEFAULT_ITERATIONS = 10_000
DEFAULT_MAX_RANDOM_HPO_TERMS = 10
DEFAULT_METHOD = "custom_jaccardIC"
DEFAULT_KIND = "omim"
DEFAULT_COMBINE = "funSimAvg"

def normalize_reference_names(names):
    names = [x.lower() for x in names]
    if "all" in names:
        return list(REFERENCE_CONFIG.keys())
    unknown = [x for x in names if x not in REFERENCE_CONFIG]
    if unknown:
        raise ValueError(
            f"Unknown reference dataset(s): {', '.join(unknown)}. "
            f"Valid choices: {', '.join(REFERENCE_CONFIG)}, all"
        )
    return list(dict.fromkeys(names))
