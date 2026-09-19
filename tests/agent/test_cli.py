"""
Tests for `stratosampler serve` host handling.

uvicorn.run is patched out, so no server is started and no network is used.
"""

import os

import pytest

pytest.importorskip("uvicorn")
pytest.importorskip("fastapi")

from click.testing import CliRunner

from stratosampler.agent import cli


@pytest.fixture
def uvicorn_calls(monkeypatch):
    import uvicorn

    calls = []
    monkeypatch.setattr(uvicorn, "run", lambda app, **kwargs: calls.append({"app": app, **kwargs}))
    return calls


@pytest.fixture
def clean_env(monkeypatch):
    # setenv first so monkeypatch restores/removes them after the CLI mutates os.environ
    monkeypatch.setenv("STRATOSAMPLER_BACKEND", "placeholder")
    monkeypatch.setenv("STRATOSAMPLER_MODEL", "placeholder")


def _serve(*args):
    return CliRunner().invoke(cli.main, ["serve", *args])


@pytest.mark.parametrize("host", ["0.0.0.0", "192.168.1.20", "::", "myserver.local"])
def test_warns_when_binding_beyond_loopback(uvicorn_calls, host):
    result = _serve("--host", host)

    assert result.exit_code == 0
    assert "no authentication" in result.output.lower()
    assert uvicorn_calls[0]["host"] == host


@pytest.mark.parametrize("host", ["127.0.0.1", "localhost", "::1"])
def test_no_warning_on_loopback(uvicorn_calls, host):
    result = _serve("--host", host)

    assert result.exit_code == 0
    assert "authentication" not in result.output.lower()
    assert uvicorn_calls[0]["host"] == host


def test_reload_passes_an_import_string_factory(uvicorn_calls, clean_env, monkeypatch):
    """uvicorn can only reload an import string, not an app object."""
    monkeypatch.setenv("GROQ_API_KEY", "gsk_test")
    result = _serve("--reload", "--backend", "groq", "--model", "some-model")

    assert result.exit_code == 0
    call = uvicorn_calls[0]
    assert call["app"] == "stratosampler.agent.server:create_app_from_env"
    assert call["factory"] is True
    assert call["reload"] is True
    assert os.environ["STRATOSAMPLER_BACKEND"] == "groq"
    assert os.environ["STRATOSAMPLER_MODEL"] == "some-model"


def test_reload_without_model_clears_model_env(uvicorn_calls, clean_env):
    _serve("--reload")

    assert os.environ["STRATOSAMPLER_BACKEND"] == "ollama"
    assert "STRATOSAMPLER_MODEL" not in os.environ


def test_without_reload_passes_the_app_object(uvicorn_calls):
    _serve()

    call = uvicorn_calls[0]
    assert not isinstance(call["app"], str)
    assert call["reload"] is False
    assert "factory" not in call


def test_groq_without_key_fails_at_startup(uvicorn_calls, monkeypatch):
    """Fail fast instead of starting a server whose every /chat returns a 500."""
    monkeypatch.delenv("GROQ_API_KEY", raising=False)

    result = _serve("--backend", "groq")

    assert result.exit_code != 0
    assert "GROQ_API_KEY" in result.output
    assert uvicorn_calls == []


def test_groq_without_key_fails_at_startup_with_reload(uvicorn_calls, clean_env, monkeypatch):
    monkeypatch.delenv("GROQ_API_KEY", raising=False)

    result = _serve("--backend", "groq", "--reload")

    assert result.exit_code != 0
    assert "GROQ_API_KEY" in result.output
    assert uvicorn_calls == []


def test_groq_with_key_starts(uvicorn_calls, monkeypatch):
    monkeypatch.setenv("GROQ_API_KEY", "gsk_test")

    result = _serve("--backend", "groq")

    assert result.exit_code == 0
    assert len(uvicorn_calls) == 1


def test_ollama_needs_no_key(uvicorn_calls, monkeypatch):
    monkeypatch.delenv("GROQ_API_KEY", raising=False)

    result = _serve()

    assert result.exit_code == 0
    assert len(uvicorn_calls) == 1
