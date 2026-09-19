"""
Tests for StratoAgent.stream(): the streamed tool-call loop.

The OpenAI-compatible client is replaced with a scripted fake, so no model,
network, or API key is needed. Tool execution is patched to a stub.
"""

import json
from types import SimpleNamespace as NS

import pytest

pytest.importorskip("openai")

from stratosampler.agent import agent as agent_module
from stratosampler.agent.agent import StratoAgent


def chunk(content=None, tool_calls=None, finish=None):
    return NS(choices=[NS(finish_reason=finish, delta=NS(content=content, tool_calls=tool_calls))])


def tool_delta(index, id=None, name=None, args=None):
    return NS(index=index, id=id, function=NS(name=name, arguments=args))


class FakeClient:
    """Serves one scripted response (list of chunks, or an Exception) per create()."""

    def __init__(self, responses):
        self._responses = list(responses)
        self.calls = []
        self.chat = NS(completions=NS(create=self._create))

    def _create(self, **kwargs):
        self.calls.append(kwargs)
        response = self._responses.pop(0)
        if isinstance(response, Exception):
            raise response
        return iter(response)


@pytest.fixture
def tool_calls(monkeypatch):
    """Record handle_tool_call invocations; each returns {"ok": true}."""
    calls = []

    def fake_handle(name, args):
        calls.append((name, args))
        return json.dumps({"ok": True})

    monkeypatch.setattr(agent_module, "handle_tool_call", fake_handle)
    return calls


def make_agent(*responses):
    agent = StratoAgent(backend="ollama")
    agent.client = FakeClient(responses)
    return agent


def test_plain_text_response():
    agent = make_agent([chunk("Hel"), chunk("lo"), chunk(finish="stop")])

    events = list(agent.stream("hi"))

    assert events == [
        {"type": "text", "content": "Hel"},
        {"type": "text", "content": "lo"},
        {"type": "done"},
    ]
    assert agent.history == [
        {"role": "user", "content": "hi"},
        {"role": "assistant", "content": "Hello"},
    ]


def test_tool_call_round_trip(tool_calls):
    agent = make_agent(
        [
            chunk(tool_calls=[tool_delta(0, id="call_1", name="load_data", args='{"path": "a.')]),
            chunk(tool_calls=[tool_delta(0, args='csv"}')]),
            chunk(finish="tool_calls"),
        ],
        [chunk("Done"), chunk(finish="stop")],
    )

    events = list(agent.stream("load a.csv"))

    assert events == [
        {"type": "tool_call", "name": "load_data", "input": {"path": "a.csv"}},
        {"type": "tool_result", "name": "load_data"},
        {"type": "text", "content": "Done"},
        {"type": "done"},
    ]
    assert tool_calls == [("load_data", {"path": "a.csv"})]

    assistant_turn, tool_turn = agent.history[1], agent.history[2]
    assert assistant_turn["tool_calls"] == [
        {
            "id": "call_1",
            "type": "function",
            "function": {"name": "load_data", "arguments": '{"path": "a.csv"}'},
        }
    ]
    assert tool_turn == {"role": "tool", "tool_call_id": "call_1", "content": '{"ok": true}'}
    assert agent.client.calls[1]["messages"][-1] == tool_turn


def test_tool_calls_run_even_when_finish_reason_is_stop(tool_calls):
    """Some OpenAI-compatible servers (older Ollama) end tool-call turns with
    finish_reason="stop"; the accumulated tool calls must still be executed."""
    agent = make_agent(
        [
            chunk(tool_calls=[tool_delta(0, id="call_1", name="load_data", args='{"path": "a.csv"}')]),
            chunk(finish="stop"),
        ],
        [chunk("Done"), chunk(finish="stop")],
    )

    events = list(agent.stream("load a.csv"))

    assert tool_calls == [("load_data", {"path": "a.csv"})]
    assert events[-1] == {"type": "done"}
    assert {"type": "text", "content": "Done"} in events


def test_multiple_interleaved_tool_calls_keep_order(tool_calls):
    agent = make_agent(
        [
            chunk(tool_calls=[tool_delta(0, id="c0", name="load_data", args='{"path":')]),
            chunk(tool_calls=[tool_delta(1, id="c1", name="compute_properties", args='{"path": "b.csv"}')]),
            chunk(tool_calls=[tool_delta(0, args=' "a.csv"}')]),
            chunk(finish="tool_calls"),
        ],
        [chunk(finish="stop")],
    )

    list(agent.stream("go"))

    assert tool_calls == [
        ("load_data", {"path": "a.csv"}),
        ("compute_properties", {"path": "b.csv"}),
    ]
    tool_ids = [m["tool_call_id"] for m in agent.history if m["role"] == "tool"]
    assert tool_ids == ["c0", "c1"]


def test_pdb_viewer_event_follows_fetch_pdb_structures(monkeypatch):
    payload = {"viewer_ready": True, "pdb_ids": ["1M17"], "structures": [{"pdb_id": "1M17"}]}
    monkeypatch.setattr(agent_module, "handle_tool_call", lambda name, args: json.dumps(payload))
    agent = make_agent(
        [
            chunk(tool_calls=[tool_delta(0, id="c0", name="fetch_pdb_structures", args='{"target_name": "EGFR"}')]),
            chunk(finish="tool_calls"),
        ],
        [chunk(finish="stop")],
    )

    events = list(agent.stream("show EGFR"))

    types = [e["type"] for e in events]
    assert types[:3] == ["tool_call", "tool_result", "pdb_viewer"]
    assert events[2] == {
        "type": "pdb_viewer",
        "pdb_ids": ["1M17"],
        "structures": [{"pdb_id": "1M17"}],
    }


def test_malformed_tool_arguments_become_empty_dict(tool_calls):
    agent = make_agent(
        [
            chunk(tool_calls=[tool_delta(0, id="c0", name="load_data", args="{not json")]),
            chunk(finish="tool_calls"),
        ],
        [chunk(finish="stop")],
    )

    events = list(agent.stream("go"))

    assert events[0] == {"type": "tool_call", "name": "load_data", "input": {}}
    assert tool_calls == [("load_data", {})]


def test_api_error_becomes_error_event():
    agent = make_agent(RuntimeError("boom"))

    assert list(agent.stream("hi")) == [{"type": "error", "content": "boom"}]


def test_reset_clears_history():
    agent = make_agent([chunk("ok"), chunk(finish="stop")])
    list(agent.stream("hi"))
    assert agent.history

    agent.reset()

    assert agent.history == []
