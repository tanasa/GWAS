import requests
import time
import pandas as pd
import os

def get_variants_as_rows(chrom, pos):
    chrom_clean = chrom.replace("chr", "")
    url = f"https://rest.ensembl.org/overlap/region/human/{chrom_clean}:{pos}-{pos}?feature=variation"
    headers = {"Content-Type": "application/json"}
    r = requests.get(url, headers=headers)

    if not r.ok:
        return [{
            "input_chrom": chrom, "input_pos": pos,
            "id": None, "seq_region_name": None, "start": None, "end": None,
            "alleles_raw": None, "alleles_list": None, "n_alleles": None,
            "is_two_allele_snv": False, "consequence_type": None,
            "minor_af": None, "minor_allele": None, "clinical_significance": None,
            "strand": None, "source": None, "variant_class": None,
            "feature_type": None, "assembly_name": None,
            "http_status": r.status_code, "error": r.text[:200]
        }]

    variants = r.json()
    if not variants:
        return [{
            "input_chrom": chrom, "input_pos": pos,
            "id": None, "seq_region_name": chrom_clean, "start": pos, "end": pos,
            "alleles_raw": None, "alleles_list": None, "n_alleles": 0,
            "is_two_allele_snv": False, "consequence_type": None,
            "minor_af": None, "minor_allele": None, "clinical_significance": None,
            "strand": None, "source": None, "variant_class": None,
            "feature_type": None, "assembly_name": None,
            "http_status": 200, "error": "No variants found at locus"
        }]

    rows = []
    for v in variants:
        alleles_raw = v.get("alleles")
        if isinstance(alleles_raw, str):
            alleles_list = alleles_raw.split("/")
        elif isinstance(alleles_raw, list):
            alleles_list = alleles_raw
        else:
            alleles_list = None

        is_two_allele_snv = False
        if alleles_list and len(alleles_list) == 2:
            a1, a2 = alleles_list[0], alleles_list[1]
            if all(isinstance(x, str) and len(x) == 1 and x in "ACGT" for x in (a1, a2)):
                if v.get("variant_class") in (None, "SNV"):
                    is_two_allele_snv = True

        rows.append({
            "input_chrom": chrom, "input_pos": pos,
            "id": v.get("id"),
            "seq_region_name": v.get("seq_region_name"),
            "start": v.get("start"), "end": v.get("end"),
            "alleles_raw": alleles_raw,
            "alleles_list": "/".join(alleles_list) if alleles_list else None,
            "n_alleles": len(alleles_list) if alleles_list else None,
            "is_two_allele_snv": is_two_allele_snv,
            "consequence_type": v.get("consequence_type"),
            "minor_af": v.get("minor_af"),
            "minor_allele": v.get("minor_allele"),
            "clinical_significance": ",".join(v.get("clinical_significance")) if isinstance(v.get("clinical_significance"), list) else v.get("clinical_significance"),
            "strand": v.get("strand"),
            "source": v.get("source"),
            "variant_class": v.get("variant_class"),
            "feature_type": v.get("feature_type"),
            "assembly_name": v.get("assembly_name"),
            "http_status": 200, "error": None
        })
    return rows


# ---------- MAIN ----------
input_file = "41467_2025_64070_MOESM4_ESM.2a_NHCFV.snp.POS.txt"

# Auto-name outputs based on input filename
basename = os.path.splitext(input_file)[0]
out_all  = f"{basename}.ensembl_all.tsv"
out_snv2 = f"{basename}.ensembl_snv2.tsv"

print(f"Input:       {input_file}")
print(f"Output all:  {out_all}")
print(f"Output snv2: {out_snv2}")
print()

COLUMNS = [
    "input_chrom", "input_pos", "id", "seq_region_name", "start", "end",
    "alleles_raw", "alleles_list", "n_alleles", "is_two_allele_snv",
    "consequence_type", "minor_af", "minor_allele", "clinical_significance",
    "strand", "source", "variant_class", "feature_type", "assembly_name",
    "http_status", "error"
]

total_all  = 0
total_snv2 = 0

with open(input_file) as f_in, \
     open(out_all,  "w") as f_all, \
     open(out_snv2, "w") as f_snv2:

    # Write headers once
    header = "\t".join(COLUMNS) + "\n"
    f_all.write(header)
    f_snv2.write(header)

    for i, line in enumerate(f_in):
        line = line.strip()
        if not line:
            continue

        chrom, pos = line.split()
        pos = int(pos)

        print(f"Processing {i+1}: {chrom}:{pos}", end="\r", flush=True)

        rows = get_variants_as_rows(chrom, pos)

        for row in rows:
            tsv_line = "\t".join(str(row.get(col, "")) for col in COLUMNS) + "\n"

            # Write to all file immediately
            f_all.write(tsv_line)
            f_all.flush()
            total_all += 1

            # Write to snv2 file only if it passes the filter
            if row["is_two_allele_snv"]:
                f_snv2.write(tsv_line)
                f_snv2.flush()
                total_snv2 += 1

        time.sleep(0.1)

print(f"\nDone!")
print(f"All rows written:        {total_all}  → {out_all}")
print(f"Two-allele SNVs written: {total_snv2} → {out_snv2}")
