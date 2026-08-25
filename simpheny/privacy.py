"""Reference-ID masking for restricted datasets."""
import hashlib
import pandas as pd

def stable_mask(source: str, raw_id: str) -> str:
    digest = hashlib.sha256(f"{source}:{raw_id}".encode()).hexdigest()[:10].upper()
    return f"{source.upper()}_REF_{digest}"

def mask_reference_ids(df: pd.DataFrame, source: str) -> pd.DataFrame:
    if df.empty:
        return df.copy()
    out = df.copy()
    out["Reference_ID"] = out["Reference_ID"].map(lambda x: stable_mask(source, str(x)))
    return out
