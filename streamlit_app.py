import streamlit as st
import re
import pandas as pd
import requests

st.set_page_config(page_title="Fine-mapping summary → FinnGen links", layout="wide")
st.title("Fine-mapping summary → FinnGen (R12) links")

RELEASE = "r12"
BASE_URL = f"https://{RELEASE}.finngen.fi"

@st.cache_data(show_spinner=False)
def get_variant_id_from_rsid(rsid: str) -> str | None:
    """Get chr:pos-ref-alt from Ensembl API (handles multi-allelic with ancestral)"""
    try:
        r = requests.get(f"https://rest.ensembl.org/variation/human/{rsid}", 
                        headers={"Content-Type": "application/json"}, timeout=6)
        r.raise_for_status()
        data = r.json()
        
        if "mappings" in data and data["mappings"]:
            mapping = data["mappings"][0]
            chr_num = mapping["seq_region_name"].replace("23", "X").replace("24", "Y")
            pos = mapping["start"]
            allele_string = mapping["allele_string"]
            ancestral = mapping["ancestral_allele"]
            
            # Parse: "G/C/T" → ref='G', alts=['C','T']
            parts = allele_string.split('/')
            if len(parts) < 2:
                return None
            ref = parts[0]
            alts = parts[1:]
            
            # Choose alt: ancestral if in alts, else first
            # FIXED LOGIC: 
            if ancestral and ancestral in alts:
                alt = ancestral  # Use ancestral if it's an alternate
            elif ancestral == ref and len(alts) > 1:
                alt = alts[1]   # Ancestral=ref → pick SECOND alt (most common)
            else:
                alt = alts[0]   # Default: first alt
            
            return f"{chr_num}:{pos}-{ref}-{alt}"
    except Exception:
        pass
    return None

def finngen_link_for_rsid(rsid: str) -> str:
    """Direct variant link or search fallback"""
    v_id = get_variant_id_from_rsid(rsid)
    if v_id:
        return f"{BASE_URL}/variant/{v_id}"
    return f"{BASE_URL}/?query={rsid}"

# === Load data ===
df = pd.read_csv("fine_mapping_summary.csv")

# === Parse ===
rs_pat = re.compile(r"(rs\d+|Affx-\d+)", re.IGNORECASE)
pip_pat = re.compile(r"\(([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\)")

def parse_top_pip_snps(cell: str):
    if not isinstance(cell, str) or not cell.strip():
        return []
    parts = [p.strip() for p in cell.split(";") if p.strip()]
    items = []
    for p in parts:
        rs_match = rs_pat.search(p)
        pip_match = pip_pat.search(p)
        if rs_match:
            rsid = rs_match.group(1)
            pip_val = float(pip_match.group(1)) if pip_match else None
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
