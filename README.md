
# omics-runner

A CLI agent that answers your scientific questions by executing e2e differential expression analysis on NCBI GEO transcriptomic data.

`omics-runner` couples a local LLM via Ollama with `GEOParse` and `PyDESeq2` to calculate differential gene expression between user-defined groups of interest and synthesize statistical findings into a natural language summary.

## Prerequisites

- **Python:** >= 3.10
- **Ollama:** Running locally (`ollama serve`)
- **System Memory:** 16 GB+ RAM recommended for expression matrix assembly

## Setup

```bash
pip install -r requirements.txt
ollama pull qwen2.5:32b
```

## Usage

Pass a specific hypothesis, biological question, or experimental comparison directly to `cli.py`:

```bash
python cli.py "Is TNF expression significantly elevated in anti-PD-1 non-responders within GSE159034?"
```

## Workflow

1. Agent extracts target GEO series accessions (GSE), metadata, and contrasts from your natural language prompt.
2. Agent queries NCBI GEO via `GEOparse` to parse sample annotations (GSM) and assemble gene-level counts or transcript abundances into an aligned matrix.
3. Agent executes differential expression analysis to derive log<sub>2</sub> fold changes, standard errors, and Benjamini-Hochberg adjusted *p*-values (*q*-values).
4. Agent evaluates the resulting DEG statistics against the original prompt and outputs a concise summary. For a request about a named human gene, the agent reolves its symbol to an Ensembl ID and looks up its statistics in the complete results. For a general request, it summarizes the run.

## Outputs

All artifacts are written to `results/` unless another output directory is specified:
- `<geo_id>_all_results_raw.csv` = Differential expression matrix across all evaluated genes using gene IDS.
- `<geo_id>_upregulated.csv` = Statistically significant upregulated genes using gene symbols (padj < 0.05, log2FC > 1.0).
- `<geo_id>_downregulated.csv` = Statistically significant downregulated genes using gene symbols (padj < 0.05, log2FC < -1.0).
