#!/usr/bin/env python3
"""Freeze the complete Wang A549 activated panel using externally defined HGNC mappings."""
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
SHA256 = "6f43d6ff43aa9fdfa5fb2f20a20a7cace66e6e02e2a0dcf19d9b726e2e248d20"
SOURCE = ROOT / "hgnc_complete_set.txt.gz"
raw = gzip.decompress(SOURCE.read_bytes())
assert hashlib.sha256(raw).hexdigest() == SHA256, "HGNC source changed; use the frozen snapshot"
hgnc = list(csv.DictReader(io.StringIO(raw.decode()), delimiter="\t"))
original = list(csv.DictReader((ROOT / "wang2004-activated-original.tsv").open(), delimiter="\t"))
assert len(original) == 50

# Original typography/compound labels require explicit source-independent expansion.
overrides = {
    "ENaCa": ("SCNN1A", "alias_symbol", "ENaCalpha"),
    "HIAP1": ("BIRC3", "alias_symbol", "hiap-1"),
    "MT-1I": ("MT1IP", "prev_symbol", "MT1I"),
    "Kip2/p57": ("CDKN1C", "alias_symbol", "KIP2|P57"),
    "Stomatin": ("STOM", "name", "stomatin"),
    "ET-2": ("EDN2", "alias_symbol", "ET2"),
    "LOC144100": ("PLEKHA7", "entrez_id", "144100"),
    "PPG/Serglycin": ("SRGN", "alias_symbol", "PPG"),
}
mapped = []
for gene in original:
    label = gene["original_label"]
    matches = []
    for field in ("symbol", "prev_symbol", "alias_symbol"):
        matches = [r for r in hgnc if label in r[field].split("|")]
        if matches:
            matched_value = label
            break
    if not matches and label in overrides:
        expected, field, matched_value = overrides[label]
        matches = [r for r in hgnc if set(matched_value.split("|")) <= set(r[field].split("|"))]
        assert len(matches) == 1 and matches[0]["symbol"] == expected, (label, matches)
    assert len(matches) == 1, (label, "ambiguous/unresolved", [r["symbol"] for r in matches])
    record = matches[0]
    assert record["status"] == "Approved"
    mapped.append(gene | {
        "symbol": record["symbol"], "hgnc_id": record["hgnc_id"],
        "entrez_id": record["entrez_id"], "ensembl_gene_id": record["ensembl_gene_id"], "hgnc_name": record["name"],
        "matching_field": field, "matching_value": matched_value,
        "mapping_policy": "explicit_typography_expansion" if label in overrides else "exact_field_match",
    })
assert len({r["symbol"] for r in mapped}) == 50
with (ROOT / "wang2004-activated-hgnc.tsv").open("w") as stream:
    writer = csv.DictWriter(stream, delimiter="\t", fieldnames=mapped[0], lineterminator="\n")
    writer.writeheader()
    writer.writerows(mapped)
provenance = {
    "article_doi": "10.1073/pnas.0407008101",
    "table_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC524211/",
    "selection": "Entire Table 1 GR-activated section; no target-study expression selection",
    "hgnc_source": URL, "hgnc_sha256": SHA256, "hgnc_bytes": len(raw),
    "hgnc_http_last_modified": "Fri, 04 Sep 2026 15:22:27 GMT",
    "declared_genes": 50, "mapped_genes": len(mapped), "ambiguous_genes": 0,
    "original_sha256": hashlib.sha256((ROOT / "wang2004-activated-original.tsv").read_bytes()).hexdigest(),
    "mapping_sha256": hashlib.sha256((ROOT / "wang2004-activated-hgnc.tsv").read_bytes()).hexdigest(),
    "source_caveats": [
        "Table repeats GenBank AA775091 for GILZ and S100P; map the gene labels, not accessions",
        "MT-1I maps to approved pseudogene MT1IP; retain in declared coverage",
        "Discovery used A549 with cycloheximide; transfer to primary airway smooth muscle is context dependent",
    ],
}
(ROOT / "hgnc-provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
print(json.dumps({k: provenance[k] for k in ("mapped_genes", "ambiguous_genes", "mapping_sha256")}))
