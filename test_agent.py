
import json

import ollama

from agent import (
    execute_deseq2_pipeline,
    inspect_geo_metadata,
    run_agent_query,
    run_agent_with_forced_json,
)
from tools import DESEQ2_TOOL_SCHEMA, INSPECT_METADATA_SCHEMA

print("Testing metadata inspection:")
inspect_res = inspect_geo_metadata("GSE159034")
print("Metadata inspection output:\n", inspect_res)

assert "potential_condition_columns" in inspect_res
print("Metadata inspection passed!")

print("Testing pipeline execution:")
pipeline_res = execute_deseq2_pipeline(
    geo_id="GSE159034",
    condition_col="characteristics_ch1.1.primary response",
    test_group="Responder",
    reference_group="Non-responder",
)
print("Pipeline execution output:\n", pipeline_res)

assert pipeline_res["status"] == "success"
print("Pipeline execution passes!")

vague_query = (
    "Can you run a DESeq2 differential expression analysis on study GSE159034? "
    "I want to compare responders versus non-responders."
)

print("Sending vague query to agent...\n")
result = run_agent_query(vague_query, model_name="qwen2.5:32b")

if not result:
    result = run_agent_with_forced_json(vague_query, model_name="qwen2.5:32b")

print("Final agent execution summary:\n", result)
assert result.get("status") == "success", "Agent failed to execute full pipeline."

def test_llm_tool_selection():
    print("Testing LLM function calling:")

    query = (
        "Run PyDESeq2 on dataset GSE159034. Use 'characteristics_ch1.1.primary response' "
        "to compare 'Responder' (test) against 'Non-responder' (reference)."
    )

    response = ollama.chat(
        model="qwen2.5:32b",
        messages=[{"role": "user", "content": query}],
        tools=[INSPECT_METADATA_SCHEMA, DESEQ2_TOOL_SCHEMA],
    )

    print("Raw Ollama response:", response)

    if hasattr(response, "message"):
        message = response.message
        tool_calls = message.tool_calls or []
    else:
        message = response.get("message", {})
        tool_calls = message.get("tools_calls", [])

    # Convert Pydantic ToolCall objects to dicts for json.dumps compatibility
    tool_calls_dict = [
        tc.model_dump() if hasattr(tc, "model_dump") else tc for tc in tool_calls
    ]

    print("Raw LLM tool call output:")
    print(json.dumps(tool_calls_dict, indent=2))

    assert len(tool_calls) > 0, "LLM failed to invoke any tool."

    func_name = (
        tool_calls[0].function.name
        if hasattr(tool_calls[0], "function")
        else tool_calls[0]["function"]["name"]
    )
    assert func_name == "run_full_deseq2_pipeline", "LLM picked the wrong tool."

    args = tool_calls[0]["function"]["arguments"]
    assert args["geo_id"] == "GSE159034"
    assert args["test_group"] == "Responder"

    print("\n LLM function calling passed!")

if __name__ == "__main__":
    test_llm_tool_selection()
