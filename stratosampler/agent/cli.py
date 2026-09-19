"""
CLI entry point: `stratosampler chat` and `stratosampler serve`.
"""

from __future__ import annotations

import ipaddress

import click


def _is_loopback(host: str) -> bool:
    if host == "localhost":
        return True
    try:
        return ipaddress.ip_address(host).is_loopback
    except ValueError:
        return False


@click.group()
def main() -> None:
    """StratoAgent — stratosampler natural language interface."""


@main.command()
@click.argument("message", required=False)
@click.option(
    "--backend",
    type=click.Choice(["ollama", "groq"]),
    default="ollama",
    show_default=True,
    help="Ollama (local, free) or Groq (hosted, needs GROQ_API_KEY).",
)
@click.option("--model", default=None, help="Model ID. Defaults to a per-backend default.")
@click.option(
    "--api-key",
    envvar="GROQ_API_KEY",
    help="Groq API key (or set GROQ_API_KEY). Ignored for --backend ollama.",
)
def chat(message: str | None, backend: str, model: str | None, api_key: str | None) -> None:
    """Chat with StratoAgent.

    Pass MESSAGE for a single-shot query, or omit to enter an interactive REPL.

    Examples:

    \b
        stratosampler chat "split my data.csv on MolLogP and TPSA, 20% test"
        stratosampler chat  # REPL mode, local Ollama backend
        stratosampler chat --backend groq "..."  # hosted, needs GROQ_API_KEY
    """
    from stratosampler.agent.agent import StratoAgent

    try:
        agent = StratoAgent(backend=backend, model=model, api_key=api_key)
    except ValueError as exc:
        raise click.ClickException(str(exc))

    def run(msg: str) -> None:
        for event in agent.stream(msg):
            if event["type"] == "text":
                print(event["content"], end="", flush=True)
            elif event["type"] == "tool_call":
                print(f"\n\033[2m[{event['name']}...]\033[0m", flush=True)
            elif event["type"] == "done":
                print()
            elif event["type"] == "error":
                content = event["content"]
                if backend == "ollama" and "Connection" in content:
                    content += "\nIs Ollama running? Try: ollama serve"
                print(f"\n\033[31mError: {content}\033[0m")

    if message:
        run(message)
        return

    click.echo("StratoAgent ready.  Type 'exit' to quit, 'reset' to clear history.\n")
    while True:
        try:
            user_input = click.prompt("\033[1m>\033[0m", prompt_suffix=" ").strip()
        except (EOFError, KeyboardInterrupt):
            click.echo()
            break
        if user_input.lower() in ("exit", "quit"):
            break
        if user_input.lower() == "reset":
            agent.reset()
            click.echo("History cleared.\n")
            continue
        if not user_input:
            continue
        run(user_input)


@main.command()
@click.option("--host", default="127.0.0.1", show_default=True, help="Bind host.")
@click.option("--port", default=8000, show_default=True, type=int, help="Bind port.")
@click.option(
    "--backend",
    type=click.Choice(["ollama", "groq"]),
    default="ollama",
    show_default=True,
    help="Ollama (local, free) or Groq (hosted, needs GROQ_API_KEY).",
)
@click.option("--model", default=None, help="Model ID. Defaults to a per-backend default.")
@click.option("--reload", is_flag=True, help="Auto-reload on code changes (dev only).")
def serve(host: str, port: int, backend: str, model: str | None, reload: bool) -> None:
    """Start the StratoAgent web server.

    Example:

    \b
        stratosampler serve --port 8080
        stratosampler serve --backend groq --port 8080
    """
    try:
        import uvicorn
    except ImportError:
        raise click.ClickException("uvicorn not installed. Run: pip install 'stratosampler[agent]'")

    from stratosampler.agent.server import create_app

    if not _is_loopback(host):
        click.secho(
            f"WARNING: binding to {host} exposes StratoAgent beyond this machine. "
            "It has no authentication, and its tools can read and write arbitrary "
            "files on the host. Use 127.0.0.1 unless the network is trusted.",
            fg="yellow",
            err=True,
        )

    app = create_app(model=model, backend=backend)
    click.echo(f"StratoAgent running at http://{host}:{port} (backend={backend})")
    uvicorn.run(app, host=host, port=port, reload=reload)
