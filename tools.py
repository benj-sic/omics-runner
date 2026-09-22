
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
