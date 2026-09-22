
DESEQ2_TOOL_SCHEMA = {
    "type": "function",
    "function": {
        "name": "run_full_deseq2_pipeline",
        "description": "Downloads GEO RNA-seq dataset, processes RSEM counts, runs differential expression using PyDESeq2, and maps gene symbols.",
        "parameters": {
            "type": "object",
            "properties": {
                "geo_id": {
                    "type": "string",
                    "description": "The GEO accession number, e.g., GSE159034",
                },
                "condition_col": {
                    "type": "string",
                    "description": "The metadata column name specifying the experimental condition or treatment.",
                },
                "test_group": {
                    "type": "string",
                    "description": "The test/experimental group name in the condition column (e.g., Responder, Treated, Mutant).",
                },
                "reference_group": {
                    "type": "string",
                    "description": "The control/reference group name in the condition column (e.g., Non-responder, Control, Wildtype).",
                },
                "padj": {
                    "type": "number",
                    "description": "FDR adjusted p-value cutoff for significance. Default is 0.05.",
                },
                "lfc": {
                    "type": "number",
                    "description": "Absolute log2 fold change cutoff. Default is 1.0.",
                },
            },
            "required": [
                "geo_id",
                "condition_col",
                "test_group",
                "reference_group",
            ],
        },
    },
}

INSPECT_METADATA_SCHEMA = {
    "type": "function",
    "function": {
        "name": "inspect_geo_metadata",
        "description": "ALWAYS run this first when the user does not provide the exact phenotype metadata column name. Returns available sample columns and conditions for a GEO accession.",
        "parameters": {
            "type": "object",
            "properties": {
                "geo_id": {
                    "type": "string",
                    "description": "The GEO accession number (e.g., 'GSE159034').",
                }
            },
            "required": ["geo_id"],
        },
    },
}

GENE_LOOKUP_SCHEMA = {
    "type": "function",
    "function": {
        "name": "lookup_gene",
        "description": (
            "Look up one human gene by its gene symbol in the current run's "
            "complete DESeq2 results. Use this tool before answering a request "
            "about a named gene. Returns its matching Ensembl IDs, log2 fold "
            "changes, p-values, adjusted p-values, and whether each row passes "
            "the run's DEG cutoffs. It can also report that the symbol could "
            "not be resolved or that the gene was absent from the result file."
        ),
        "parameters": {
            "type": "object",
            "properties": {
                "gene_symbol": {
                    "type": "string",
                    "description": (
                        "A single human gene symbol to look up, such as TNF. "
                        "Do not include a file path or a comparison group."
                    ),
                }
            },
            "required": ["gene_symbol"],
            "additionalProperties": False,
        },
    },
}
