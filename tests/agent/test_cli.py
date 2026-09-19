"""
Tests for `stratosampler serve` host handling.

uvicorn.run is patched out, so no server is started and no network is used.
"""

import pytest

pytest.importorskip("uvicorn")
pytest.importorskip("fastapi")

from click.testing import CliRunner

from stratosampler.agent import cli


@pytest.fixture
def uvicorn_calls(monkeypatch):
    import uvicorn

    calls = []
    monkeypatch.setattr(uvicorn, "run", lambda app, **kwargs: calls.append(kwargs))
    return calls


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
