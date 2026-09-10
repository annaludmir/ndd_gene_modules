import re
import pandas as pd

# 1. Target terms list
TARGET_TERMS = ["intracellular non-membrane-bounded organelle (GO:0043232)","Regulation of sister chromatid separation at the metaphase-anaphase transition WP4240","regulation of cell cycle G2/M phase transition (GO:1902749)","regulation of G2/M transition of mitotic cell cycle (GO:0010389)","mitotic spindle organization (GO:0007052)","microtubule cytoskeleton organization involved in mitosis (GO:1902850)","chromosome (GO:0005694)","mitotic sister chromatid segregation (GO:0000070)","condensed chromosome (GO:0000793)","regulation of mitotic cell cycle phase transition (GO:1901990)"
]


def parse_genes(gene_str):
    """Splits gene string by ;, , or whitespace and returns a set of cleaned gene symbols."""
    if pd.isna(gene_str):
        return set()
    # Split by semicolon, comma, slash, or space
    genes = re.split(r"[;,/\s]+", str(gene_str))
    return {g.strip().upper() for g in genes if g.strip()}


def extract_genes_union(
    csv_paths,
    target_terms,
    term_col="Term",
    genes_col="Genes",
    output_txt="union_genes.txt",
):
    """Extracts the union of genes corresponding to specified GO terms across multiple CSV files."""
    # Normalize target terms for matching
    normalized_targets = {t.strip().lower(): t for t in target_terms}

    all_union_genes = set()
    term_gene_map = {t: set() for t in target_terms}
    file_gene_counts = {}

    for path in csv_paths:
        df = pd.read_csv(path)

        if term_col not in df.columns or genes_col not in df.columns:
            raise ValueError(
                f"File '{path}' must contain columns '{term_col}' and '{genes_col}'."
            )

        file_genes = set()

        for _, row in df.iterrows():
            raw_term = str(row[term_col]).strip()
            term_key = raw_term.lower()

            if term_key in normalized_targets:
                original_term_name = normalized_targets[term_key]
                genes = parse_genes(row[genes_col])

                all_union_genes.update(genes)
                file_genes.update(genes)
                term_gene_map[original_term_name].update(genes)

        file_gene_counts[path] = len(file_genes)

    # Convert sorted list for output
    sorted_genes = sorted(all_union_genes)

    # Save to TXT file (one gene per line)
    with open(output_txt, "w") as f:
        f.write("\n".join(sorted_genes))

    # Print Summary Reports
    print("=" * 60)
    print("EXTRACTION SUMMARY")
    print("=" * 60)
    for path, count in file_gene_counts.items():
        print(f"File '{path}': {count} matching genes extracted.")

    print(
        f"\nTotal unique genes in UNION across all files: {len(sorted_genes)}"
    )
    print(f"Saved gene list to: {output_txt}\n")

    print("-" * 60)
    print("GENE COUNT PER GO TERM")
    print("-" * 60)
    for term in target_terms:
        count = len(term_gene_map[term])
        print(f"{term}: {count} genes")

    return sorted_genes, term_gene_map


if __name__ == "__main__":
    csv_files =["/miridan-data/annaludmir/ndd_gene_modules/results/GO_terms/microcephaly_GO_enrichment/CellCyclePhase_-_S/significant_hits_fdr_0.05.csv",
"/miridan-data/annaludmir/ndd_gene_modules/results/GO_terms/microcephaly_GO_enrichment/CellCyclePhase_-_G2M/significant_hits_fdr_0.05.csv",
    "/miridan-data/annaludmir/ndd_gene_modules/results/GO_terms/microcephaly_GO_enrichment/CellCyclePhase_-_PostM/significant_hits_fdr_0.05.csv"]

    union_genes, term_map = extract_genes_union(
        csv_paths=csv_files,
        target_terms=TARGET_TERMS,
        term_col="Term",
        genes_col="Genes",
        output_txt="S_G2M_PostM_union_mitotic_go_terms_genes.txt",
    )
