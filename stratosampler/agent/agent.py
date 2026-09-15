"""
StratoAgent — Groq-powered agent wrapping stratosampler tools.
"""

from __future__ import annotations

import json
import os
from typing import Iterator

from stratosampler.agent.prompts import SYSTEM_PROMPT
from stratosampler.agent.tools import OPENAI_TOOL_SCHEMAS, handle_tool_call

_SYSTEM = [{"role": "system", "content": SYSTEM_PROMPT}]

_BACKENDS = {
    "ollama": {
        "base_url": "http://localhost:11434/v1",
        "default_api_key": "ollama",
        "default_model": "llama3.2:3b",
    },
    "groq": {
        "base_url": "https://api.groq.com/openai/v1",
        "api_key_env": "GROQ_API_KEY",
        "default_model": "llama-3.3-70b-versatile",
    },
}


def _resolve_backend(
    backend: str, model: str | None, api_key: str | None
) -> tuple[str, str, str]:
    """Resolve (base_url, api_key, model) for a given backend name."""
    if backend not in _BACKENDS:
        raise ValueError(f"Unknown backend {backend!r}. Choose from {list(_BACKENDS)}.")
    cfg = _BACKENDS[backend]
    resolved_model = model or cfg["default_model"]

    if backend == "groq":
        resolved_key = api_key or os.environ.get(cfg["api_key_env"])
        if not resolved_key:
            raise ValueError(
                "GROQ_API_KEY not set and no api_key provided for backend='groq'."
            )
    else:
        resolved_key = api_key or cfg["default_api_key"]

    return cfg["base_url"], resolved_key, resolved_model


class StratoAgent:
    """Stateful agent with conversation history.

    Parameters
    ----------
    backend : str
        "ollama" (default, local, free) or "groq" (hosted, needs GROQ_API_KEY).
    model : str, optional
        Model ID. Defaults to a per-backend default (llama3.2:3b for Ollama,
        llama-3.3-70b-versatile for Groq).
    api_key : str, optional
        API key. Ignored for Ollama. Falls back to GROQ_API_KEY env var for Groq.
    max_tokens : int
        Max tokens per response turn.
    """

    def __init__(
        self,
        backend: str = "ollama",
        model: str | None = None,
        api_key: str | None = None,
        max_tokens: int = 4096,
    ):
        from openai import OpenAI

        base_url, resolved_key, resolved_model = _resolve_backend(backend, model, api_key)
        self.client = OpenAI(base_url=base_url, api_key=resolved_key)
        self.backend = backend
        self.model = resolved_model
        self.max_tokens = max_tokens
        self.history: list[dict] = []

    def reset(self) -> None:
        self.history.clear()

    def stream(self, message: str) -> Iterator[dict]:
        """Streaming generator. Yields dicts with keys:

        {"type": "text",        "content": str}
        {"type": "tool_call",   "name": str, "input": dict}
        {"type": "tool_result", "name": str}
        {"type": "pdb_viewer",  "pdb_ids": list, "structures": list}
        {"type": "done"}
        {"type": "error",       "content": str}
        """
        self.history.append({"role": "user", "content": message})

        while True:
            try:
                completion = self.client.chat.completions.create(
                    model=self.model,
                    max_tokens=self.max_tokens,
                    messages=_SYSTEM + self.history,
                    tools=OPENAI_TOOL_SCHEMAS,
                    tool_choice="auto",
                    stream=True,
                )

                accumulated_text = ""
                tool_call_chunks: dict[int, dict] = {}
                finish_reason = None

                for chunk in completion:
                    choice = chunk.choices[0]
                    finish_reason = choice.finish_reason or finish_reason
                    delta = choice.delta

                    if delta.content:
                        accumulated_text += delta.content
                        yield {"type": "text", "content": delta.content}

                    if delta.tool_calls:
                        for tc in delta.tool_calls:
                            idx = tc.index
                            if idx not in tool_call_chunks:
                                tool_call_chunks[idx] = {
                                    "id": "",
                                    "type": "function",
                                    "function": {"name": "", "arguments": ""},
                                }
                            if tc.id:
                                tool_call_chunks[idx]["id"] = tc.id
                            if tc.function:
                                if tc.function.name:
                                    tool_call_chunks[idx]["function"]["name"] += tc.function.name
                                if tc.function.arguments:
                                    tool_call_chunks[idx]["function"]["arguments"] += tc.function.arguments

            except Exception as exc:
                yield {"type": "error", "content": str(exc)}
                return

            if not tool_call_chunks or finish_reason == "stop":
                self.history.append({"role": "assistant", "content": accumulated_text or ""})
                yield {"type": "done"}
                return

            # Tool calls — add assistant turn then execute each tool
            tool_calls = [tool_call_chunks[i] for i in sorted(tool_call_chunks)]
            self.history.append({
                "role": "assistant",
                "content": accumulated_text or None,
                "tool_calls": tool_calls,
            })

            for tc in tool_calls:
                name = tc["function"]["name"]
                try:
                    args = json.loads(tc["function"]["arguments"])
                except json.JSONDecodeError:
                    args = {}

                yield {"type": "tool_call", "name": name, "input": args}
                result = handle_tool_call(name, args)
                yield {"type": "tool_result", "name": name}

                if name == "fetch_pdb_structures":
                    try:
                        rd = json.loads(result)
                        if rd.get("viewer_ready") and rd.get("pdb_ids"):
                            yield {
                                "type": "pdb_viewer",
                                "pdb_ids": rd["pdb_ids"],
                                "structures": rd.get("structures", []),
                            }
                    except Exception:
                        pass

                self.history.append({
                    "role": "tool",
                    "tool_call_id": tc["id"],
                    "content": result,
                })
