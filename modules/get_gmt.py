import os
import pandas as pd
from pathlib import Path

def save_to_gmt(gene_list_path: str, output_path: str) -> str:
    """
    Create a GMT file from a gene list file.

    Parameters
    ----------
    gene_list_path : str
        Path to a file containing one gene per line (txt/csv).
    output_dir : str
        Folder where the GMT file will be created.

    Returns
    -------
    gmt_file : str
        Path to the created GMT file.
    """
    gene_list_path = Path(gene_list_path)

    if not gene_list_path.exists():
      raise FileNotFoundError(f"Gene list file not found: {gene_list_path}")

    # Load gene list
    df = pd.read_csv(gene_list_path)
    genes = df['gene'].dropna().astype(str).unique().tolist()

    if len(genes) == 0:
        raise ValueError(f"❌ Gene list is empty: {gene_list_path}")

    gene_set_name = gene_list_path.stem

    # Write GMT file
    print(f"📄 Creating GMT file: {output_path}")

    with open(output_path, "w") as f:
        genes_joined = "\t".join(genes)
        line = f"{gene_set_name}\tDescription\t{genes_joined}\n"
        f.write(line)

    print(f"✔ GMT created successfully ({len(genes)} genes)")
    return str(output_path)
