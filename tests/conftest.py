"""Test-suite guards that must apply before any test imports run.

**The suite must never touch the user's real protein cache.** `proteins.build()`
defaults to `~/.teaid/proteins`, and any test that reaches the CLI's protein
path on a machine with BATH installed would otherwise fetch 130 Pfam models
from InterPro and write ~240 MB into it — turning a seconds-long suite into a
multi-minute one and, worse, leaving a cache the developer did not ask for.

So `TEAID_CACHE` is pointed at a throwaway directory unless the developer has
set it deliberately. Setting it to the prebuilt cache is the normal way to
exercise the real library::

    TEAID_CACHE=dev-data/protein-cache pytest
"""

from __future__ import annotations

import os

import pytest


@pytest.fixture(scope="session", autouse=True)
def _isolate_protein_cache(tmp_path_factory):
    """Point the protein cache somewhere disposable for the whole session."""
    if os.environ.get("TEAID_CACHE"):
        yield  # the developer chose one; respect it
        return
    home = tmp_path_factory.mktemp("teaid-cache")
    previous = os.environ.get("TEAID_CACHE")
    os.environ["TEAID_CACHE"] = str(home)
    try:
        yield
    finally:
        if previous is None:
            os.environ.pop("TEAID_CACHE", None)
        else:
            os.environ["TEAID_CACHE"] = previous
