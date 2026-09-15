"""
Tests for stratosampler.agent.agent's backend resolution (Ollama/Groq).

Pure-function tests only — no network calls, no real openai.OpenAI client
construction, no Groq/Ollama server required.
"""

import pytest

from stratosampler.agent.agent import _resolve_backend


class TestResolveBackendOllama:
    def test_defaults(self):
        base_url, api_key, model = _resolve_backend("ollama", model=None, api_key=None)
        assert base_url == "http://localhost:11434/v1"
        assert api_key == "ollama"
        assert model == "llama3.2:3b"

    def test_custom_model_overrides_default(self):
        _, _, model = _resolve_backend("ollama", model="qwen2.5:7b-instruct", api_key=None)
        assert model == "qwen2.5:7b-instruct"

    def test_base_url_unaffected_by_api_key_arg(self):
        base_url, api_key, _ = _resolve_backend("ollama", model=None, api_key="whatever")
        assert base_url == "http://localhost:11434/v1"
        assert api_key == "whatever"


class TestResolveBackendGroq:
    def test_uses_env_var_when_no_explicit_key(self, monkeypatch):
        monkeypatch.setenv("GROQ_API_KEY", "gsk_from_env")
        base_url, api_key, model = _resolve_backend("groq", model=None, api_key=None)
        assert base_url == "https://api.groq.com/openai/v1"
        assert api_key == "gsk_from_env"
        assert model == "llama-3.3-70b-versatile"

    def test_explicit_api_key_overrides_env(self, monkeypatch):
        monkeypatch.setenv("GROQ_API_KEY", "gsk_from_env")
        _, api_key, _ = _resolve_backend("groq", model=None, api_key="gsk_explicit")
        assert api_key == "gsk_explicit"

    def test_custom_model_overrides_default(self, monkeypatch):
        monkeypatch.setenv("GROQ_API_KEY", "gsk_from_env")
        _, _, model = _resolve_backend("groq", model="llama-3.1-8b-instant", api_key=None)
        assert model == "llama-3.1-8b-instant"

    def test_missing_key_raises(self, monkeypatch):
        monkeypatch.delenv("GROQ_API_KEY", raising=False)
        with pytest.raises(ValueError, match="GROQ_API_KEY"):
            _resolve_backend("groq", model=None, api_key=None)


class TestResolveBackendUnknown:
    def test_unknown_backend_raises(self):
        with pytest.raises(ValueError, match="ollama"):
            _resolve_backend("bogus", model=None, api_key=None)
