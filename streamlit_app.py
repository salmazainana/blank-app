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
    """Resolve rsID → chr:pos-ref-alt via Ensembl; return None if not resolvable quickly."""
    try:
        r = requests.get(
            f"https://rest.ensembl.org/variation/human/{rsid}",
            headers={"Content-Type": "application/json"},
            timeout=6,
        )
        r.raise_for_status()
        data = r.json()
        if "mappings" in data and data["mappings"]:
            m = data["mappings"][0]
            chrom = str(m["seq_region_name"]).replace("23", "X").replace("24", "Y")
            pos = m["start"]
            allele_string = m.get("allele_string", "")
            ancestral = m.get("ancestral_allele")

            parts = allele_string.split("/")
            if len(parts) < 2:
                return None
            ref = parts[0]
            alts = parts[1:]
            if ancestral and ancestral in alts:
                alt = ancestral
            elif ancestral == ref and len(alts) > 1:
                alt = alts[1]
            else:
                alt = alts[0]
            return f"{chrom}:{pos}-{ref}-{alt}"
    except Exception:
        pass
    return None

def finngen_link_for_rsid(rsid: str) -> str:
    v_id = get_variant_id_from_rsid(rsid)
    if v_id:
        return f"{BASE_URL}/variant/{v_id}"
    return f"{BASE_URL}/?query={rsid}"

# ===== Load data =====
df = pd.read_csv("fine_mapping_summary.csv")

# ===== Robust parser for Top_PIP_SNPs =====
# Accept rsIDs, Affx IDs, or already-normalized variant IDs (chr:pos-ref-alt)
RS_PATTERN = re.compile(r"(rs\d+|Affx-\d+|\d+:\d+-[ACGT]+-[ACGT]+)", re.IGNORECASE)
NUM_PATTERN = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")

def parse_top_pip_snps_to_dict(cell: str) -> dict:
    """
    Turn 'rs123 (0.91, 3.2e-10); rs456 (pip=0.40; p=0.05)' into
    {'rs123': [0.91, 3.2e-10], 'rs456': [0.40, 0.05]}
    Missing p-values become None. Missing PIP become None (row will be ignored when ranking).
    """
    d = {}
    if not isinstance(cell, str) or not cell.strip():
        return d
    # split by semicolons at top level
    for chunk in [c.strip() for c in cell.split(";") if c.strip()]:
        # find SNP id
        m_id = RS_PATTERN.search(chunk)
        if not m_id:
            continue
        snp = m_id.group(1)

        # extract content inside parentheses (if any)
        pip_val, p_val = None, None
        lpar = chunk.find("(")
        rpar = chunk.rfind(")")
        inside = chunk[lpar+1:rpar] if (lpar != -1 and rpar != -1 and rpar > lpar) else ""

        if inside:
            low = inside.lower()

            # named forms first
            m_pip_named = re.search(r"pip\s*=\s*(" + NUM_PATTERN.pattern + r")", low)
            m_p_named   = re.search(r"\b(pval|p)\s*=\s*(" + NUM_PATTERN.pattern + r")", low)

            if m_pip_named:
                try: pip_val = float(m_pip_named.group(1))
                except: pass
            if m_p_named:
                try: p_val = float(m_p_named.group(2))
                except: pass

            # positional fallback if still missing
            if pip_val is None or p_val is None:
                nums = [n.group(0) for n in NUM_PATTERN.finditer(inside)]
                if pip_val is None and len(nums) >= 1:
                    try: pip_val = float(nums[0])
                    except: pass
                if p_val is None and len(nums) >= 2:
                    try: p_val = float(nums[1])
                    except: pass

        d[snp] = [pip_val, p_val]
    return d

# ===== Build output rows =====
rows = []
for _, r in df.iterrows():
    d = parse_top_pip_snps_to_dict(r.get("Top_PIP_SNPs", ""))
    # pick SNP with highest PIP (ignore None); tie-breaker = smaller p-value
    best_snp, best_pip, best_p = None, None, None
    if d:
        # sort by (-pip, p) where None pip ranks last; None p ranks after numeric p
        def sort_key(item):
            snp, (pip, p) = item
            k1 = float("inf") if pip is None else -pip
            k2 = float("inf") if p is None else p
            return (k1, k2)
        best_snp, (best_pip, best_p) = sorted(d.items(), key=sort_key)[0]

    rows.append({
        "GWAS_File": r.get("GWAS_File"),
        "Lead_SNP": r.get("Lead_SNP"),
        "CS_Size": r.get("CS_Size"),
        "top_pip_snps": best_snp,     # requested name
        "pip": best_pip,              # requested name
        "pvalue": best_p,             # requested name
        "FinnGen": finngen_link_for_rsid(best_snp) if isinstance(best_snp, str) and best_snp.lower().startswith("rs") else None,
    })

out = pd.DataFrame(rows)

# Optional: clip PIP to [0,1] if you want to suppress >1 due to formatting noise
# out["pip"] = pd.to_numeric(out["pip"], errors="coerce").clip(lower=0, upper=1)

# ===== Display =====
st.subheader("Summary")
st.dataframe(
    out[["GWAS_File", "Lead_SNP", "CS_Size", "top_pip_snps", "pip", "pvalue", "FinnGen"]],
    use_container_width=True,
    column_config={
        "pip": st.column_config.NumberColumn(format="%.3f", help="Posterior inclusion probability (PIP)"),
        "pvalue": st.column_config.NumberColumn(format="%.2e", help="p-value for the selected top-PIP SNP"),
        "FinnGen": st.column_config.LinkColumn(display_text="Open in FinnGen (R12)"),
    }
)

st.caption(
    "Parsed Top_PIP_SNPs → {snp: [pip, pvalue]}; sorted by PIP desc (ties by smaller p). "
    "FinnGen links are only built for rsIDs; chr:pos-ref-alt links can be used directly if provided."
)

