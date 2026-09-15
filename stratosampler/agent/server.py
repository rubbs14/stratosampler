"""
FastAPI server for StratoAgent.  Exposes SSE streaming at POST /chat.
"""

from __future__ import annotations

import asyncio
import json
import threading
import time
import uuid
from pathlib import Path

from fastapi import FastAPI
from fastapi.responses import HTMLResponse, StreamingResponse
from pydantic import BaseModel

_STATIC = Path(__file__).parent / "static"
_SESSION_TTL_SECONDS = 1800  # evict sessions idle longer than this


class ChatRequest(BaseModel):
    message: str
    session_id: str | None = None


def create_app(model: str | None = None, backend: str = "ollama") -> FastAPI:
    from stratosampler.agent.agent import StratoAgent

    app = FastAPI(title="StratoAgent", docs_url=None, redoc_url=None)
    sessions: dict[str, tuple[StratoAgent, float]] = {}

    def _purge_expired() -> None:
        cutoff = time.time() - _SESSION_TTL_SECONDS
        for sid in [s for s, (_, last_used) in sessions.items() if last_used < cutoff]:
            del sessions[sid]

    @app.get("/", response_class=HTMLResponse)
    async def index() -> str:
        return (_STATIC / "index.html").read_text()

    @app.post("/chat")
    async def chat(request: ChatRequest) -> StreamingResponse:
        _purge_expired()
        session_id = request.session_id or str(uuid.uuid4())
        agent, _ = sessions.get(session_id, (None, 0.0))
        if agent is None:
            agent = StratoAgent(backend=backend, model=model)
        sessions[session_id] = (agent, time.time())

        async def event_stream():
            q: asyncio.Queue = asyncio.Queue()
            loop = asyncio.get_running_loop()
            sentinel = object()

            def produce() -> None:
                try:
                    for event in agent.stream(request.message):
                        asyncio.run_coroutine_threadsafe(q.put(event), loop).result()
                except Exception as exc:
                    asyncio.run_coroutine_threadsafe(
                        q.put({"type": "error", "content": str(exc)}), loop
                    ).result()
                finally:
                    asyncio.run_coroutine_threadsafe(q.put(sentinel), loop).result()

            thread = threading.Thread(target=produce, daemon=True)
            thread.start()

            # First event: send session_id so client can persist it
            yield f"data: {json.dumps({'type': 'session_id', 'session_id': session_id})}\n\n"

            while True:
                event = await q.get()
                if event is sentinel:
                    break
                yield f"data: {json.dumps(event)}\n\n"

        return StreamingResponse(
            event_stream(),
            media_type="text/event-stream",
            headers={
                "Cache-Control": "no-cache",
                "X-Accel-Buffering": "no",
            },
        )

    @app.delete("/sessions/{session_id}")
    async def reset_session(session_id: str) -> dict:
        if session_id in sessions:
            agent, _ = sessions[session_id]
            agent.reset()
        return {"ok": True}

    return app
