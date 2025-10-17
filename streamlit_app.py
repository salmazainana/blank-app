import streamlit as st
import re
import pandas as pd
import requests
import os

st.set_page_config(page_title="Fine-mapping summary → FinnGen links", layout="wide")
st.title("Fine-mapping summary → FinnGen (R12) links")

RELEASE = "r12"
BASE_URL = f"https://{RELEASE}.finngen.fi"

@st.cache_data(show_spinner=False)
def get_variant_id_from_rsid(rsid: str) -> str | None:
    """
    Resolve rsID → chr:pos-ref-alt using Ensembl REST.
    Returns None if not resolvable quickly.
    """
    try:
        r = requests.get(
            f"https://rest.ensembl.org/variation/human/{rsid}",
            headers={"Content-Type": "application/json"},
            timeout=6,
        )
        r.raise_for_status()
        data = r.json()
        if "mappings" in data and data["mappings"]:
            mapping = data["mappings"][0]
            chr_num = str(mapping["seq_region_name"]).replace("23", "X").replace("24", "Y")
            pos = mapping["start"]
            allele_string = mapping.get("allele_string", "")
            ancestral = mapping.get("ancestral_allele")

            parts = allele_string.split("/")
            if len(parts) < 2:
                return None
            ref = parts[0]
            alts = parts[1:]

            # Choose alt: if ancestral equals an alt, prefer it; if ancestral==ref and multiple alts, pick second; else first
            if ancestral and ancestral in alts:
                alt = ancestral
            elif ancestral == ref and len(alts) > 1:
                alt = alts[1]
            else:
                alt = alts[0]
            return f"{chr_num}:{pos}-{ref}-{alt}"
    except Exception:
        pass
    return None

def finngen_link_for_rsid(rsid: str) -> str:
    """Direct variant page when resolvable; else fall back to search."""
    v_id = get_variant_id_from_rsid(rsid)
    if v_id:
        return f"{BASE_URL}/variant/{v_id}"
    return f"{BASE_URL}/?query={rsid}"

# === Load data ===
df = pd.read_csv("fine_mapping_summary.csv")

# === Parsers ===
rs_pat = re.compile(r"(rs\d+|Affx-\d+)", re.IGNORECASE)
pip_pat = re.compile(r"\(([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\)")

def parse_top_pip_snps(cell: str):
    """Parse 'Top_PIP_SNPs' like 'rs123 (0.91); rs456 (0.40)' into [{'rsid', 'pip'}, ...] sorted by PIP desc."""
    if not isinstance(cell, str) or not cell.strip():
        return []
    parts = [p.strip() for p in cell.split(";") if p.strip()]
    items = []
    for p in parts:
        m_rsid = rs_pat.search(p)
        m_pip = pip_pat.search(p)
        if m_rsid:
            rsid = m_rsid.group(1)
            pip_val = float(m_pip.group(1)) if m_pip else None
            items.append({"rsid": rsid, "pip": pip_val})
    items.sort(key=lambda d: (-d["pip"] if isinstance(d["pip"], (float, int)) else float("-inf")))
    return items

def parse_top_p_values(cell: str):
    """
    Optional: parse a sibling column 'Top_P_values' like 'rs123 (3.2e-10); rs456 (0.04)'
    Returns dict {rsid -> pval}
    """
    if not isinstance(cell, str) or not cell.strip():
        return {}
    parts = [p.strip() for p in cell.split(";") if p.strip()]
    out = {}
    for p in parts:
        m_rsid = rs_pat.search(p)
        m_val = pip_pat.search(p)  # same pattern handles scientific notation
        if m_rsid and m_val:
            try:
                out[m_rsid.group(1)] = float(m_val.group(1))
            except ValueError:
                continue
    return out

# === Build table of top SNPs ===
rows = []
for _, r in df.iterrows():
    parsed = parse_top_pip_snps(r["Top_PIP_SNPs"])
    top = parsed[0] if parsed else {"rsid": None, "pip": None}
    link = finngen_link_for_rsid(top["rsid"]) if top["rsid"] else None

    rows.append({
        "GWAS_File": r.get("GWAS_File"),
        "Lead_SNP": r.get("Lead_SNP"),
        "CS_Size": r.get("CS_Size"),
        "Top_SNP": top["rsid"],
        "Top_SNP_PIP": top["pip"],
        "FinnGen": link,
    })

out = pd.DataFrame(rows)

# === Add Top_SNP_P (p-value) from available sources ===
out["Top_SNP_P"] = pd.NA

# Source A: parse from a 'Top_P_values' column if present
if "Top_P_values" in df.columns:
    # Build a per-row dict rsid->p, then pick for Top_SNP
    p_dicts = df["Top_P_values"].apply(parse_top_p_values)
    out["Top_SNP_P"] = [
        d.get(rsid) if isinstance(d, dict) and isinstance(rsid, str) else pd.NA
        for d, rsid in zip(p_dicts, out["Top_SNP"])
    ]

# Source B: optional p-value map file (GWAS_File, rsid, pval)
if out["Top_SNP_P"].isna().any() and os.path.exists("pval_map.csv"):
    pmap = pd.read_csv("pval_map.csv")  # expects columns: GWAS_File, rsid, pval
    pmap = pmap.rename(columns={"rsid": "Top_SNP", "pval": "Top_SNP_P"})
    out = out.merge(pmap[["GWAS_File", "Top_SNP", "Top_SNP_P"]], how="left", on=["GWAS_File", "Top_SNP"])
    # Prefer existing non-null, else take merged
    out["Top_SNP_P"] = out["Top_SNP_P_x"].combine_first(out["Top_SNP_P_y"])
    out = out.drop(columns=[c for c in out.columns if c.endswith("_x") or c.endswith("_y")])

# Final niceties
# Clip any tiny numeric drift in PIP, and format
if "Top_SNP_PIP" in out.columns:
    out["Top_SNP_PIP"] = pd.to_numeric(out["Top_SNP_PIP"], errors="coerce").clip(lower=0, upper=1)

# === Display ===
st.subheader("Summary")
st.dataframe(
    out[["GWAS_File", "Lead_SNP", "CS_Size", "Top_SNP", "Top_SNP_PIP", "Top_SNP_P", "FinnGen"]],
    use_container_width=True,
    column_config={
        "Top_SNP_PIP": st.column_config.NumberColumn(format="%.3f", help="Posterior inclusion probability"),
        "Top_SNP_P": st.column_config.NumberColumn(format="%.2e", help="p-value of the top-PIP SNP"),
        "FinnGen": st.column_config.LinkColumn(display_text="Open in FinnGen (R12)"),
    }
)

st.caption(
    "P-values are taken from 'Top_P_values' when present; otherwise from optional 'pval_map.csv' "
    "(columns: GWAS_File, rsid, pval). Links prefer the canonical R12 variant page resolved via Ensembl; "
    "if resolution fails, they fall back to the search page."
)
