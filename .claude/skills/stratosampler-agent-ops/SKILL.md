---
name: stratosampler-agent-ops
description: Operate StratoAgent, the natural-language layer for stratosampler (in stratosampler/agent/) — installing the `agent` extras, choosing a backend (local Ollama by default, free and unrate-limited, or hosted Groq opt-in via --backend groq), running `stratosampler chat` (interactive REPL or one-shot) and `stratosampler serve` (FastAPI + SSE), calling the streaming POST /chat endpoint with curl or Python, resetting a session via DELETE /sessions/{id}, and fixing "uvicorn not installed" / Ollama-connection / GROQ_API_KEY errors. Use this whenever the user wants to run, deploy, script against, or debug StratoAgent, `stratosampler chat`, or `stratosampler serve` — even if they just say "start the agent server" or "why won't chat work." This is about *operating* the agent CLI/server, not about the underlying molecular dataset splitting library (load_data, split_dataset, etc.) — for that, see the stratosampler-split skill.
---

# Operating StratoAgent

StratoAgent is the natural-language wrapper around stratosampler's tools
(`load_data`, `split_dataset`, `compare_strategies`, `visualize_split`,
`compute_properties`, `fetch_pdb_structures`, `find_mcs`). It ships as a CLI
(`stratosampler chat` / `stratosampler serve`) plus a FastAPI SSE server, all
under `stratosampler/agent/`. This skill is about running and debugging that
layer — not about the splitting library itself.

## 1. Install the agent extras

The base `stratosampler` install does **not** pull in the `openai` client
(used for both backends), FastAPI, uvicorn, or click — they live behind the
`agent` extras group in `pyproject.toml`:

```bash
pip install -e '.[agent,rdkit]'
```

Quote the extras (`'.[agent]'`) — an unquoted `.[agent]` gets glob-expanded by
zsh and fails. Verify the install succeeded with:

```bash
python -c "import uvicorn, fastapi, openai, click; print('agent extras OK')"
```

If any of those four imports fail, the extras weren't installed (or were
installed into a different environment/interpreter than the one running
`stratosampler`).

## 2. Pick a backend

Two backends, chosen with `--backend` (default: `ollama`):

- **`ollama`** (default) — fully local and free, no rate limits. Needs
  [Ollama](https://ollama.com) installed and running (`ollama serve`, or
  the app's background service) plus a model pulled:
  ```bash
  ollama pull llama3.2:3b
  ```
  No API key needed — StratoAgent talks to `http://localhost:11434/v1`.
- **`groq`** (opt-in) — hosted, stronger default model
  (`llama-3.3-70b-versatile`), but Groq's free tier is request-rate-limited
  and StratoAgent makes **one API call per tool call**, not per user
  message — a multi-step task (load → split → visualize) burns through
  that quota fast. Needs a key from console.groq.com:
  ```bash
  export GROQ_API_KEY="gsk_..."
  ```

Both `chat` and `serve` read `GROQ_API_KEY` automatically when
`--backend groq` is used. You can override it per-invocation of `chat` with
`--api-key`, but `serve` has no such flag — it always relies on the
environment variable being set in the process that runs `stratosampler serve`.
`--api-key` and `GROQ_API_KEY` are ignored entirely under `--backend ollama`.

## 3. `stratosampler chat`

Two modes, same command:

**One-shot** — pass the message as an argument, get one streamed reply, exit:

```bash
stratosampler chat "split my data.csv on MolLogP and TPSA, 20% test"
```

**Interactive REPL** — omit the message:

```bash
stratosampler chat
```

Inside the REPL:
- Type `exit` or `quit` to leave.
- Type `reset` to clear conversation history (calls `agent.reset()`) without
  restarting the process.
- Empty input is ignored (re-prompts).

Useful flags on `chat` (also apply conceptually to `serve`):
- `--backend` — `ollama` (default) or `groq`.
- `--model` — model ID, defaults to `llama3.2:3b` on Ollama or
  `llama-3.3-70b-versatile` on Groq.
- `--api-key` — overrides `GROQ_API_KEY` for this invocation (Groq only).

While streaming, tool calls print dimmed as `[tool_name...]` before their
result is folded back into the conversation — that's expected, not an error.

## 4. `stratosampler serve`

Starts the FastAPI app under uvicorn:

```bash
stratosampler serve --host 127.0.0.1 --port 8000                       # local, Ollama
stratosampler serve --backend groq --model llama-3.3-70b-versatile     # hosted, Groq
```

- `--reload` auto-restarts on code changes — dev only, don't use it in
  anything resembling production.
- The server also serves a minimal built-in web UI at `GET /` — open
  `http://127.0.0.1:8000/` in a browser for a point-and-click chat window
  instead of curl/Python.
- The server has **no authentication** and its tools read/write arbitrary host
  paths. Keep `--host 127.0.0.1`; `serve` prints a warning for any non-loopback
  host (e.g. `0.0.0.0`) — don't suggest exposing it without a trusted network
  or a reverse proxy providing auth.
- `/docs` and `/redoc` are intentionally disabled (`docs_url=None`,
  `redoc_url=None` in `create_app`), so don't expect Swagger UI there.
- Conversations live in an **in-memory** `dict` inside the server process —
  they are lost on restart and are not shared across multiple worker
  processes. Idle sessions are evicted automatically after 30 minutes
  (`_SESSION_TTL_SECONDS = 1800`), checked lazily on each `/chat` call.

## 5. Talk to the SSE `/chat` endpoint directly

`POST /chat` takes JSON `{"message": str, "session_id": str | None}` and
responds with a `text/event-stream`. The **first** event is always
`{"type": "session_id", "session_id": "<uuid>"}` — capture it and send it
back on subsequent requests to keep the conversation's history. Omitting
`session_id` (or sending an id the server doesn't recognize, e.g. after TTL
eviction) just starts a fresh agent under that id.

Remaining event `type`s mirror `StratoAgent.stream()`: `text`, `tool_call`,
`tool_result`, `pdb_viewer`, `done`, `error`.

### curl

```bash
curl -N -X POST http://127.0.0.1:8000/chat \
  -H "Content-Type: application/json" \
  -d '{"message": "split my data.csv on MolLogP and TPSA, 20% test"}'
```

`-N` disables curl's output buffering so SSE chunks print as they arrive
instead of all at once at the end.

### Python

```python
import json
import requests

resp = requests.post(
    "http://127.0.0.1:8000/chat",
    json={"message": "split my data.csv on MolLogP and TPSA, 20% test"},
    stream=True,
)

session_id = None
for line in resp.iter_lines(decode_unicode=True):
    if not line or not line.startswith("data: "):
        continue
    event = json.loads(line[len("data: "):])
    if event["type"] == "session_id":
        session_id = event["session_id"]
    elif event["type"] == "text":
        print(event["content"], end="", flush=True)
    elif event["type"] == "tool_call":
        print(f"\n[{event['name']}...]")
    elif event["type"] == "error":
        print(f"\nError: {event['content']}")

# Reuse session_id to keep history across turns:
resp2 = requests.post(
    "http://127.0.0.1:8000/chat",
    json={"message": "now visualize that split", "session_id": session_id},
    stream=True,
)
```

## 6. Reset a session

```bash
curl -X DELETE http://127.0.0.1:8000/sessions/<session_id>
```

This calls `agent.reset()` on that session's history — the `session_id`
stays valid and reusable afterward (the endpoint does not delete the entry
from the sessions dict). It always returns `{"ok": true}`, even for an
unknown or already-expired `session_id`, so a 200 response doesn't guarantee
the id was actually tracked.

## 7. Troubleshooting

### "uvicorn not installed" ImportError

Running `stratosampler serve` and seeing:

```
Error: uvicorn not installed. Run: pip install 'stratosampler[agent]'
```

means the `agent` extras weren't installed into the environment/interpreter
that `stratosampler` is running from. Fix:

```bash
pip install -e '.[agent]'
```

Then re-verify with the import check from step 1. If `pip install` reports
success but the error persists, you likely have more than one Python
environment in play (e.g. installed into a venv but running a
globally-installed `stratosampler` script, or vice versa) — check
`which stratosampler` and `python -c "import stratosampler; print(stratosampler.__file__)"`
resolve to the environment you just installed into.

### Related but different: missing `openai` package

`cli.py` only special-cases uvicorn's `ImportError` inside `serve`. The
`openai` import happens lazily inside `StratoAgent.__init__`, used by both
`chat` and `serve`, so a missing `openai` package surfaces as a raw
`ModuleNotFoundError` traceback rather than a friendly click error. Same
fix applies: `pip install -e '.[agent]'`.

### Connection error on `--backend ollama` (the default)

Means Ollama isn't reachable at `http://localhost:11434`. Start it
(`ollama serve`, or launch the Ollama app) and make sure the requested
model has been pulled — `ollama pull llama3.2:3b` (or whatever `--model`
you passed). The CLI appends an "Is Ollama running?" hint to this error
when it detects a connection failure.

### No / bad GROQ_API_KEY (only relevant to `--backend groq`)

`StratoAgent.__init__` raises `ValueError: GROQ_API_KEY not set and no
api_key provided for backend='groq'` immediately at construction if
neither `GROQ_API_KEY` nor `--api-key` is set — before any chat turn is
attempted. Confirm the key is visible in the same shell/process that runs
`stratosampler`: `echo $GROQ_API_KEY`. This check does not apply to
`--backend ollama` at all.

### A long-lived session "forgets" everything

If a `session_id` sits idle for more than 30 minutes, the server evicts it
on the next `/chat` call and silently starts a brand-new agent under that
same id — there's no error, just a fresh, empty history. If a client needs
guaranteed persistence beyond 30 minutes of inactivity, it must send a
keep-alive message periodically or treat `session_id` as best-effort.
