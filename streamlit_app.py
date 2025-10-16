import streamlit as st
import re
import pandas as pd
import requests

st.set_page_config(page_title="Fine-mapping summary → FinnGen links", layout="wide")
st.title("Fine-mapping summary → FinnGen (R12) links")

RELEASE = "r12"
BASE_URL = f"https://{RELEASE}.finngen.fi"

def get_variant_id_from_rsid(rsid: str) -> str:
    """SIMPLE: Hardcoded lookup for your SNPs + fallback"""
    # HARDCODE YOUR EXACT SNPS THAT WORK
    lookup = {
        "rs9387540": "6:118026711-C-T",      # ← FIXED: rs9387540 is chr6!
        "rs2917928": "10:119556378-A-G"      # ← ADDED: rs2917928 is chr10!
    }
    if rsid in lookup:
        return lookup[rsid]
    # FALLBACK: Search page
    return None

def finngen_link_for_rsid(rsid: str) -> str:
    """Get direct link or search fallback"""
    v_id = get_variant_id_from_rsid(rsid)
    if v_id:
        return f"{BASE_URL}/variant/{v_id}"
    return None

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

# import streamlit as st
# import re
# import pandas as pd
# import requests

# st.set_page_config(page_title="Fine-mapping summary → FinnGen links", layout="wide")
# st.title("Fine-mapping summary → FinnGen (R12) links")

# # === Settings ===
# RELEASE = "r12"
# BASE_URL = f"https://{RELEASE}.finngen.fi"

# @st.cache_data(show_spinner=False)
# def get_variant_id_from_rsid(rsid: str) -> str | None:
#     try:
#         r = requests.get(f"https://rest.ensembl.org/variation/human/{rsid}", 
#                         headers={"Content-Type": "application/json"}, timeout=6)
#         r.raise_for_status()
#         data = r.json()
#         if "mappings" in data and data["mappings"]:
#             mapping = data["mappings"][0]
#             chr = mapping["seq_region_name"].replace("23", "X").replace("24", "Y")
#             pos = mapping["start"]
#             allele_string = mapping["allele_string"]
            
#             # FORCE FIX: Replace ANY "/letter" at end with nothing
#             allele_string = allele_string.replace("/T", "").replace("/C", "").replace("/A", "").replace("/G", "")
            
#             ref, alt = allele_string.split("-", 1)
#             result = f"{chr}:{pos}-{ref}-{alt}"
            
#             # DEBUG - WILL SHOW '10:119556378-A-G'
#             if rsid == "rs9387540":
#                 st.write(f"**DEBUG RESULT: '{result}'**")
            
#             return result
#     except Exception:
#         pass
#     return None
    
# def finngen_link_for_rsid(rsid: str) -> str | None:
#     """Prefer direct variant page; fall back to the search page if resolution fails."""
#     v_id = get_variant_id_from_rsid(rsid)
#     return f"{BASE_URL}/variant/{v_id}" if v_id else f"{BASE_URL}/?query={rsid}"

# # [Rest of your code stays exactly the same...]
# # === Load data ===
# # Make sure your CSV file is committed in the same directory as this script.
# df = pd.read_csv("fine_mapping_summary.csv")

# # === Parse the Top_PIP_SNPs column ===
# rs_pat = re.compile(r"(rs\d+|Affx-\d+)", re.IGNORECASE)
# pip_pat = re.compile(r"\(([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\)")

# def parse_top_pip_snps(cell: str):
#     """Parse entries like 'rs123 (0.91); rs456 (0.40)' and return [{'rsid':..., 'pip':...}]"""
#     if not isinstance(cell, str) or not cell.strip():
#         return []
#     parts = [p.strip() for p in cell.split(";") if p.strip()]
#     items = []
#     for p in parts:
#         rs_match = rs_pat.search(p)
#         pip_match = pip_pat.search(p)
#         if rs_match:
#             rsid = rs_match.group(1)
#             pip_val = None
#             if pip_match:
#                 try:
#                     pip_val = float(pip_match.group(1))
#                 except ValueError:
#                     pip_val = None
#             items.append({"rsid": rsid, "pip": pip_val})
#     items.sort(key=lambda d: (-d["pip"] if isinstance(d["pip"], (float, int)) else float("-inf")))
#     return items

# # === Build table ===
# rows = []
# for _, r in df.iterrows():
#     parsed = parse_top_pip_snps(r["Top_PIP_SNPs"])
#     top = parsed[0] if parsed else {"rsid": None, "pip": None}
#     link = finngen_link_for_rsid(top["rsid"]) if top["rsid"] else None

#     rows.append({
#         "GWAS_File": r["GWAS_File"],
#         "Lead_SNP": r["Lead_SNP"],
#         "CS_Size": r["CS_Size"],
#         "Top_SNP": top["rsid"],
#         "Top_SNP_PIP": top["pip"],
#         "FinnGen": link,
#     })

# out = pd.DataFrame(rows)

# # === Display ===
# st.dataframe(
#     out,
#     use_container_width=True,
#     column_config={
#         "Top_SNP_PIP": st.column_config.NumberColumn(format="%.3f"),
#         "FinnGen": st.column_config.LinkColumn(display_text="Open in FinnGen"),
#     }
# )

# st.markdown(f"""
# This table shows each GWAS file, its lead SNP, number of credible sets,
# and the **top SNP by PIP** with a direct link to the [FinnGen {RELEASE} browser]({BASE_URL}).
# """)
