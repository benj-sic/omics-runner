
#!/usr/bin/env python3

import argparse
import glob
import logging
import os
import sys

import GEOparse
import mygene
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from pathlib import Path

# Logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler("pipeline.log"), logging.StreamHandler(sys.stdout)],
)

# Pipeline Functions
def fetch_geo_data(geo_id: str, output_dir: str = "./raw_data") -> tuple[GEOparse.GSE, str]:
    dataset_dir = os.path.join(output_dir, geo_id)
    supp_files_dir = os.path.join(dataset_dir, "supp_files")

    # Get GEO SOFT file
    gse = GEOparse.get_GEO(geo=geo_id, destdir=dataset_dir)

    # Get supplementary files
    gse.download_supplementary_files(directory=supp_files_dir, download_sra=False)

    return gse, supp_files_dir

def build_expression_matrix(supp_dir: str, target_metric: str = "expected_count", min_total_reads: int = 10,) -> pd.DataFrame:
    dataframes = []

    file_pattern = os.path.join(supp_dir, "**", "*.rsem.genes.results.txt.gz")
    rsem_files = glob.glob(file_pattern, recursive = True)

    if not rsem_files:
        raise FileNotFoundError(f"No RSEM files found matching {file_pattern}")

    for filepath in rsem_files:
        filename = os.path.basename(filepath)

        # Get sample id from filename
        sample_id = filename.split("_")[0]

        # Read only gene id + target metric cols
        df = pd.read_csv(filepath, sep="\t", usecols = ["gene_id", target_metric])

        # Set gene id as index and rename metric col to sample id
        df = df.set_index("gene_id").rename(columns={target_metric: sample_id})
        dataframes.append(df)

    # Concat horizontally
    expression_matrix = pd.concat(dataframes, axis=1)

    # Transpose and convert to rounded ints for PyDeseq2
    counts_df = expression_matrix.T.round().astype(int)

    filtered_counts_df = counts_df.loc[:, counts_df.sum(axis=0) >= min_total_reads]

    return filtered_counts_df

def run_deseq2_analysis(counts_df: pd.DataFrame, phenotype_df: pd.DataFrame, condition_col: str, contrast: list[str], n_cpus: int = 4) -> pd.DataFrame:
    # Align metadata rows to match count matrix rows
    aligned_metadata = phenotype_df.loc[counts_df.index]

    # Configure parallel inference handler
    inference = DefaultInference(n_cpus=n_cpus)

    # Initialize Deseq2 dataset object with design
    dds = DeseqDataSet(
        counts = counts_df,
        metadata = aligned_metadata,
        design = f"~`{condition_col}`",
        refit_cooks=True,
        inference=inference,
    )

    # Run dispersion estimation and log2FC fitting
    dds.deseq2()

    # Calculate hypothesis tests for the design contrast
    stat_res = DeseqStats(dds, contrast=contrast)
    stat_res.summary()

    # Extract raw results dataframe
    results_df = pd.DataFrame(stat_res.results_df)

    return results_df

def filter_deg_results(results_df: pd.DataFrame, padj_thresh: float = 0.05, lfc_thresh: float = 1.0,) -> tuple[pd.DataFrame, pd.DataFrame]:
    # Filter out ns and missing rows
    sig_genes = results_df.dropna(subset=["padj"]).query(f"padj < {padj_thresh}")

    # Extract and sort upregulated genes
    upregulated = sig_genes[sig_genes["log2FoldChange"] > lfc_thresh].sort_values(by="padj", ascending=True)

    # Extract and sort downregulated genes
    downregulated = sig_genes[sig_genes["log2FoldChange"] < -lfc_thresh].sort_values(by="padj", ascending=True)

    return upregulated, downregulated

def annotate_ensembl_ids(results_df: pd.DataFrame, species: str = "human") -> pd.DataFrame:
    df = results_df.copy()

    # Strip Ensembl version tags
    clean_ids = [ensg.split(".")[0] for ensg in df.index]

    # Bulk REST API query
    mg = mygene.MyGeneInfo()
    query_results = mg.querymany(
        clean_ids, scopes="ensembl.gene", fields="symbol", species=species
    )

    # Construct map with Ensembl ID fallback
    mapping = {
        hit["query"]: hit.get("symbol", hit["query"]) for hit in query_results
    }

    # Map symbols to df index
    df.index = [mapping.get(i, i) for i in clean_ids]

    return df

def analyze_and_save(
    geo_id,
    condition_col,
    test_group,
    reference_group,
    outdir,
    padj,
    lfc,
    n_cpus
):
    # Fetch
    gse, supp_dir = fetch_geo_data(geo_id=geo_id)

    # Build matrix
    counts_df = build_expression_matrix(supp_dir=supp_dir)

    # Deseq2
    raw_results = run_deseq2_analysis(
        counts_df=counts_df,
        phenotype_df=gse.phenotype_data,
        condition_col=condition_col,
        contrast=[condition_col, test_group, reference_group],
        n_cpus=n_cpus,
    )

    # Filter
    raw_up, raw_down = filter_deg_results(
        raw_results, padj_thresh=padj, lfc_thresh=lfc
    )

    # Annotate
    upregulated = annotate_ensembl_ids(raw_up)
    downregulated = annotate_ensembl_ids(raw_down)

    # Save results
    output_dir = Path(outdir)
    output_dir.mkdir(parents=True, exist_ok=True)

    raw_path = output_dir / f"{geo_id}_all_results_raw.csv"
    up_path = output_dir / f"{geo_id}_upregulated.csv"
    down_path = output_dir / f"{geo_id}_downregulated.csv"

    raw_results.to_csv(raw_path)
    upregulated.to_csv(up_path)
    downregulated.to_csv(down_path)

    logging.info(f"Pipeline finished successfully! Output written to: {outdir}")

    return {
        "status": "success",
        "geo_id": geo_id,
        "test_group": test_group,
        "reference_group": reference_group,
        "total_upregulated": len(upregulated),
        "total_downregulated": len(downregulated),
        "padj_threshold": padj,
        "lfc_threshold": lfc,
        "files": {
            "raw": str(raw_path),
            "upregulated": str(up_path),
            "downregulated": str(down_path),
        },
    }

# Main execution
def main ():
    parser = argparse.ArgumentParser(description="Run DESeq2 pipeline on GEO datasets.")
    parser.add_argument("--geo_id", type=str, required=True, help="GEO series ID (e.g., GSE159034)")
    parser.add_argument("--condition", type=str, required=True, help="Metadata colname for comparison condition")
    parser.add_argument("--contrast", nargs=2, required=True, help= "Test group and reference group (e.g., Responder Non-responder)")
    parser.add_argument("--outdir", type=str, default="./results", help="Directory to save output files")
    parser.add_argument("--cpus", type=int, default=4, help="Number of CPU cores for PyDESeq2")
    parser.add_argument("--padj", type=float, default=0.05, help="FDR threshold")
    parser.add_argument("--lfc", type=float, default=1.0, help="Log2 Fold Change threshold")

    args = parser.parse_args()

    try:
        result = analyze_and_save(
            geo_id=args.geo_id,
            condition_col=args.condition,
            test_group=args.contrast[0],
            reference_group=args.contrast[1],
            outdir=args.outdir,
            n_cpus=args.cpus,
            padj=args.padj,
            lfc=args.lfc)

    except Exception as e:
        logging.critical(f"Pipeline execution crashed: {e}", exc_info=True)
        sys.exit(1)

    print("Files written:")
    for label, path in result["files"].items():
        print(f"{label}: {path}")

if __name__ == "__main__":
    main()
