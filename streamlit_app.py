import streamlit as st
import re
import pandas as pd
import requests

st.set_page_config(page_title="Fine-mapping summary → FinnGen links", layout="wide")
st.title("Fine-mapping summary → FinnGen (R12) links")

# === Settings ===
RELEASE = "r12"
BASE_URL = f"https://{RELEASE}.finngen.fi"

@st.cache_data(show_spinner=False)
def get_variant_id_from_rsid(rsid: str) -> str | None:
    """
    Query Ensembl REST API to resolve rsID to chr:pos-ref-alt (FinnGen format).
    Handles multi-allelic variants by taking the first allele pair.
    """
    try:
        r = requests.get(f"https://rest.ensembl.org/variation/human/{rsid}", 
                        headers={"Content-Type": "application/json"}, timeout=6)
        r.raise_for_status()
        data = r.json()
        if "mappings" in data and data["mappings"]:
            mapping = data["mappings"][0]
            chr = mapping["seq_region_name"].replace("23", "X").replace("24", "Y")
            pos = mapping["start"]
            allele_string = mapping["allele_string"]
            
            # Handle multi-allelic: split A-G/T → take first pair A-G
            if "/" in allele_string:
                allele_string = allele_string.split("/")[0]  # A-G/T → A-G
            
            ref, alt = allele_string.split("-", 1)
            return f"{chr}:{pos}-{ref}-{alt}"  # FinnGen format: colon, not dash
    except Exception:
        pass
    return None

def finngen_link_for_rsid(rsid: str) -> str | None:
    """Prefer direct variant page; fall back to the search page if resolution fails."""
    v_id = get_variant_id_from_rsid(rsid)
    return f"{BASE_URL}/variant/{v_id}" if v_id else f"{BASE_URL}/?query={rsid}"

# === Load data ===
# Make sure your CSV file is committed in the same directory as this script.
df = pd.read_csv("fine_mapping_summary.csv")

# === Parse the Top_PIP_SNPs column ===
rs_pat = re.compile(r"(rs\d+|Affx-\d+)", re.IGNORECASE)
pip_pat = re.compile(r"\(([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\)")

def parse_top_pip_snps(cell: str):
    """Parse entries like 'rs123 (0.91); rs456 (0.40)' and return [{'rsid':..., 'pip':...}]"""
    if not isinstance(cell, str) or not cell.strip():
        return []
    parts = [p.strip() for p in cell.split(";") if p.strip()]
    items = []
    for p in parts:
        rs_match = rs_pat.search(p)
        pip_match = pip_pat.search(p)
        if rs_match:
            rsid = rs_match.group(1)
            pip_val = None
            if pip_match:
                try:
                    pip_val = float(pip_match.group(1))
                except ValueError:
                    pip_val = None
            items.append({"rsid": rsid, "pip": pip_val})
    items.sort(key=lambda d: (-d["pip"] if isinstance(d["pip"], (float, int)) else float("-inf")))
    return items

# === Build table ===
rows = []
for _, r in df.iterrows():
    parsed = parse_top_pip_snps(r["Top_PIP_SNPs"])
    top = parsed[0] if parsed else {"rsid": None, "pip": None}
    link = finngen_link_for_rsid(top["rsid"]) if top["rsid"] else None

    rows.append({
        "GWAS_File": r["GWAS_File"],
        "Lead_SNP": r["Lead_SNP"],
        "CS_Size": r["CS_Size"],
        "Top_SNP": top["rsid"],
        "Top_SNP_PIP": top["pip"],
        "FinnGen": link,
    })

out = pd.DataFrame(rows)

# === Display ===
st.dataframe(
    out,
    use_container_width=True,
    column_config={
        "Top_SNP_PIP": st.column_config.NumberColumn(format="%.3f"),
        "FinnGen": st.column_config.LinkColumn(display_text="Open in FinnGen"),
    }
)

st.markdown(f"""
This table shows each GWAS file, its lead SNP, number of credible sets,
and the **top SNP by PIP** with a direct link to the [FinnGen {RELEASE} browser]({BASE_URL}).
""")
