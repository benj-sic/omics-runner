
from argparse import Action

from agent import execute_deseq2_pipeline, inspect_geo_metadata, ActionDecision
import ollama

from tools import INSPECT_METADATA_SCHEMA

import json

def prompt_user_for_contrast(metadata: dict):
    candidate_columns = metadata.get("potential_condition_columns", {})

    col_keys = list(candidate_columns.keys())

    print("Candidate columns:")
    for col_index, col in enumerate(candidate_columns, 1):
        print(f"{col_index}. {col} ({', '.join(candidate_columns[col])})")

    col_choice = input("Select the number corresponding to the condition column: ")
    while not col_choice.isdigit() or not (1 <= int(col_choice) <= len(col_keys)):
        col_choice = input(f"Invalid! Enter a number between 1 and {len(col_keys)}: ")

    selected_col = col_keys[int(col_choice) - 1]

    values_list = candidate_columns[selected_col]

    print("Groups:")
    for value_index, value in enumerate(values_list, 1):
        print(f"{value_index}. {value}")

    test_choice = input("Select the number corresponding to the test group: ")
    while not test_choice.isdigit() or not (1 <= int(test_choice) <= len(values_list)):
        test_choice = input(f"Invalid! Enter a number between 1 and {len(values_list)}: ")

    selected_test = values_list[int(test_choice) - 1]

    reference_choice = input("Select the number corresponding to the reference: ")
    while (
        not reference_choice.isdigit()
        or not (1 <= int(reference_choice) <= len(values_list))
        or int(reference_choice) == int(test_choice)
    ):
        if reference_choice.isdigit() and int(reference_choice) == int(test_choice):
            reference_choice = input("Reference group cannot be the same as the test group! Enter a different number: ")
        else:
            reference_choice = input(f"Invalid! Enter a number between 1 and {len(values_list)}: ")

    selected_reference = values_list[int(reference_choice) - 1]

    user_choices = {
        "condition_col": selected_col,
        "test_group": selected_test,
        "reference_group": selected_reference
    }

    return user_choices

def main():
    print("Initiating omics-runner...")
    user_prompt = str(input("What analysis would you like to perform? "))

    print("Thinking...")
    response = ollama.chat(
        model="qwen2.5:32b",
        messages=[
            {"role": "system", "content": (
                "Extract GEO ID and contrast settings. If the exact metadata column name or "
                "contrast groups are missing or ambiguous, set tool_name to 'inspect_geo_metadata'."
            )
            },
            {"role": "user", "content": user_prompt},
        ],
        format=ActionDecision.model_json_schema()
    )

    action = ActionDecision.model_validate_json(response["message"]["content"])

    if (
        action.tool_name == "inspect_geo_metadata"
        or not action.condition_col
        or not action.test_group
        or not action.reference_group
    ):
        print("Fetching metadata...")

        fetched_metadata = inspect_geo_metadata(action.geo_id)
        user_contrast = prompt_user_for_contrast(fetched_metadata)
        user_condition = user_contrast["condition_col"]
        user_test = user_contrast["test_group"]
        user_ref = user_contrast["reference_group"]
    else:
        user_condition = action.condition_col
        user_test = action.test_group
        user_ref = action.reference_group

    results = execute_deseq2_pipeline(
        geo_id=action.geo_id,
        condition_col=user_condition,
        test_group = user_test,
        reference_group = user_ref,
        padj = 0.05,
        lfc = 1.0
    )

    print("DESeq2 Pipeline Results:")
    print(json.dumps(results, indent=2))

    print("LLM Summary:")
    summary_response = ollama.chat(
        model="qwen2.5:32b",
        messages=[
            {
                "role": "system",
                "content": "You are a bioinformatics assistant. Summarize DESeq2 results clearly for scientists.",
            },
            {
                "role": "user",
                "content": f"Summarize these differential expression results:\n{json.dumps(results)}"
            }
        ]
    )
    print(summary_response["message"]["content"])

if __name__ == "__main__":
    main()
