"""Tests for the server's env-driven app factory (used by `serve --reload`)."""

import pytest

pytest.importorskip("fastapi")

from stratosampler.agent import server


@pytest.fixture
def captured(monkeypatch):
    seen = {}

    def fake_create_app(model=None, backend="ollama"):
        seen.update(model=model, backend=backend)
        return "APP"

    monkeypatch.setattr(server, "create_app", fake_create_app)
    return seen


def test_factory_reads_backend_and_model_from_env(monkeypatch, captured):
    monkeypatch.setenv("STRATOSAMPLER_BACKEND", "groq")
    monkeypatch.setenv("STRATOSAMPLER_MODEL", "some-model")

    assert server.create_app_from_env() == "APP"
    assert captured == {"model": "some-model", "backend": "groq"}


def test_factory_defaults_when_env_unset(monkeypatch, captured):
    monkeypatch.delenv("STRATOSAMPLER_BACKEND", raising=False)
    monkeypatch.delenv("STRATOSAMPLER_MODEL", raising=False)

    server.create_app_from_env()

    assert captured == {"model": None, "backend": "ollama"}
