# StratoAgent

StratoAgent is a natural-language interface over stratosampler. It wraps the
library's core operations — loading data, splitting, comparing strategies,
visualizing, computing properties, MCS analysis, and fetching PDB structures
— as tools an LLM can call, exposed via a CLI and a small web app.

By default it runs entirely **local and free**, against
[Ollama](https://ollama.com). A hosted [Groq](https://console.groq.com)
backend is available as an opt-in for when you want a bigger model's
better tool-calling — but Groq's free tier is request-rate-limited, and
each StratoAgent tool call is its own API round-trip, so a multi-step task
burns through that quota fast. Ollama has no such limit.

---

## Installation

```bash
pip install "stratosampler[agent,rdkit]"
```

This pulls in the extra dependencies the agent needs on top of the base
library: `openai` (used as a generic client for both backends — Ollama and
Groq each expose an OpenAI-compatible endpoint), `fastapi`, `uvicorn`, `click`.

For the default local backend, also install
[Ollama](https://ollama.com/download) itself and pull a model:

```bash
ollama pull llama3.2:3b
```

For development (editable install from a checkout):

```bash
pip install -e ".[agent,rdkit]"
```

`rdkit` is included alongside `agent` because most tools (splitting, MCS,
property computation) need it to parse SMILES.

---

## Configuration

There is no stratosampler-side config file — everything is passed as CLI
flags or environment variables.

| Setting | How to set it | Default |
|---|---|---|
| Backend | `--backend ollama\|groq` | `ollama` |
| Model | `--model` | `llama3.2:3b` (ollama) / `llama-3.3-70b-versatile` (groq) |
| Groq API key | `GROQ_API_KEY` env var, or `--api-key` (chat only) | — (required only for `--backend groq`) |
| Max response tokens | `StratoAgent(max_tokens=...)` (Python API only) | `4096` |
| Server host/port | `--host` / `--port` (serve only) | `127.0.0.1` / `8000` |
| Session idle TTL | hardcoded `_SESSION_TTL_SECONDS` in `agent/server.py` | `1800` (30 min) |

`--backend ollama` needs nothing beyond `ollama serve` running locally (the
default install starts this automatically). `--backend groq` needs a key
from [console.groq.com](https://console.groq.com):

```bash
export GROQ_API_KEY="gsk_..."
```

Any model your backend hosts with tool-calling support can be passed via
`--model` — e.g. `qwen2.5:7b-instruct` for Ollama, or
`llama-3.1-8b-instant` for Groq.

---

## Running it

### Interactive chat (CLI)

```bash
stratosampler chat
```

```
StratoAgent ready.  Type 'exit' to quit, 'reset' to clear history.

> split my data.csv on MolLogP and TPSA, scaffold-aware, 20% test
```

Or one-shot, non-interactive:

```bash
stratosampler chat "split molecules.csv on MolLogP and TPSA, 20% test"
```

Useful flags:

```bash
stratosampler chat --model qwen2.5:7b-instruct "quick property summary of data.csv"
stratosampler chat --backend groq "split molecules.csv on MolWt and TPSA"  # hosted, needs GROQ_API_KEY
stratosampler chat --backend groq --api-key gsk_... "..."   # overrides GROQ_API_KEY
```

Inside the REPL, `reset` clears conversation history, `exit`/`quit` ends the
session.

### Web server

```bash
stratosampler serve --port 8080
stratosampler serve --backend groq --port 8080  # hosted, needs GROQ_API_KEY
```

Serves a small web UI at `http://127.0.0.1:8080/` and a streaming chat API
at `POST /chat` (Server-Sent Events). `--reload` auto-restarts on code
changes during development.

**Security:** the server has no authentication, and the agent's tools can read and
write arbitrary files on the host (`load_data`, `split_dataset`, `visualize_split`).
Keep the default `--host 127.0.0.1`. Binding to `0.0.0.0` or any non-loopback
address prints a warning, because anyone who can reach the port can use those tools.

Talking to the API directly:

```bash
curl -N -X POST http://127.0.0.1:8080/chat \
  -H "Content-Type: application/json" \
  -d '{"message": "split molecules.csv on MolWt and TPSA"}'
```

The first SSE event returns a `session_id`; pass it back on subsequent
requests (`{"message": "...", "session_id": "..."}`) to keep conversation
history. Sessions idle for more than 30 minutes are evicted automatically.

To clear a session's history without losing the id:

```bash
curl -X DELETE http://127.0.0.1:8080/sessions/<session_id>
```

---

## Testing

### Automated

```bash
pip install -e ".[agent,dev,rdkit]"
pytest tests/agent/ -v
```

These tests cover `stratosampler/agent/tools.py`'s handlers and the backend
resolution logic in `agent.py`. They mock the Groq/Ollama HTTP calls and
RCSB network calls — no running Ollama, no `GROQ_API_KEY`, and no internet
access are required to run them.

### Manual smoke test

To confirm your actual local setup (Ollama + a pulled model) works
end-to-end:

```bash
ollama list                    # confirm llama3.2:3b (or your --model) is pulled
stratosampler chat "load examples/egfr_stratified_sample.csv, smiles column is 'smiles'"
```

Expect a dimmed `[load_data...]` line (the tool call) followed by the
model's summary of the dataset's shape and columns. If instead you get a
connection error, see [Troubleshooting](#troubleshooting) below.

For the hosted backend:

```bash
export GROQ_API_KEY="gsk_..."
stratosampler chat --backend groq "load examples/egfr_stratified_sample.csv, smiles column is 'smiles'"
```

---

## Available tools

| Tool | Purpose |
|---|---|
| `load_data` | Load a CSV/SDF dataset and report shape, columns, sample rows |
| `split_dataset` | Run `PropertyStratifiedSplitter`, save train/test(/val) CSVs + metrics |
| `compare_strategies` | Benchmark random vs. stratified vs. scaffold-aware splitting |
| `visualize_split` | Save a property-distribution plot comparing train vs. test |
| `compute_properties` | Summary stats (mean/std/min/max) for molecular properties |
| `find_mcs` | Maximum Common Substructure across a list of SMILES |
| `fetch_pdb_structures` | Look up 3D structures for a protein target on RCSB PDB |

See [`stratosampler-split`](https://github.com/rubbs14/stratosampler/blob/main/.claude/skills/stratosampler-split/SKILL.md),
[`stratosampler-mcs`](https://github.com/rubbs14/stratosampler/blob/main/.claude/skills/stratosampler-mcs/SKILL.md), and
[`stratosampler-agent-ops`](https://github.com/rubbs14/stratosampler/blob/main/.claude/skills/stratosampler-agent-ops/SKILL.md)
for detailed workflows if you're using Claude Code against this repo.

---

## Claude Code skills

This repo ships three [Claude Code](https://claude.com/claude-code) skills
under `.claude/skills/` that teach an agent working in this codebase how to
run splits, do MCS analysis, and operate StratoAgent itself — using the
real function signatures, not guesses.

### Project-scoped (default, no install step)

Skills under a repo's `.claude/skills/<name>/SKILL.md` are auto-discovered
by Claude Code whenever it's run inside that repo (or a subdirectory of
it). There's nothing to register — clone the repo, open Claude Code there,
and `stratosampler-split`, `stratosampler-mcs`, and `stratosampler-agent-ops`
are available automatically. Claude picks the right one based on each
skill's `description` frontmatter, or you can invoke one explicitly:

```
Skill(skill="stratosampler-split")
```

### Making them available in other projects

To use these skills outside this repo, copy the skill directories into your
user-level skills folder:

```bash
cp -r .claude/skills/stratosampler-split      ~/.claude/skills/
cp -r .claude/skills/stratosampler-mcs        ~/.claude/skills/
cp -r .claude/skills/stratosampler-agent-ops  ~/.claude/skills/
```

Skills in `~/.claude/skills/` load for every project, not just this one.

### Verifying they're loaded

Start a Claude Code session in the repo and ask what skills are available,
or open a `SKILL.md` directly to confirm its `name`/`description`
frontmatter — that description is what Claude matches against your request
to decide whether to trigger it.

---

## Troubleshooting

**`uvicorn not installed`** — `stratosampler serve` needs the `agent` extra:
`pip install "stratosampler[agent]"`.

**`ModuleNotFoundError: No module named 'openai'`** — same fix; this one
isn't caught by a friendly `ClickException`, so it surfaces as a raw
traceback if `agent` extras weren't installed.

**Connection error on `--backend ollama` (default)** — Ollama isn't
running, or the model hasn't been pulled. Run `ollama serve` in another
terminal and `ollama pull llama3.2:3b` (or your `--model` choice) first.

**Agent errors immediately with an auth-related message on `--backend
groq`** — `GROQ_API_KEY` is missing or invalid. Check `echo $GROQ_API_KEY`,
or pass `--api-key` explicitly to `chat`.

**Split/property tools fail on valid-looking SMILES** — make sure `rdkit`
is installed (`pip install "stratosampler[rdkit]"`); most tools import it
lazily and fail with an import error otherwise.
