"""
Tiny results registry: every analysis script writes the numbers the
manuscript quotes to results/<script>.json, and make_numbers.py turns
all of them into LaTeX macros (numbers.tex) that article.tex \\input's.
The paper therefore cannot quote a number the code did not produce.

Value convention: a plain number, a string, or a (number, fmt) pair
where fmt is a Python format spec (e.g. ".4f") controlling how the
macro renders. Plain numbers get a sensible default in make_numbers.py.
"""

import json
import os

RESULTS_DIR = "results"


def write_results(name, values):
    """Write `values` (dict of key -> number | str | (number, fmt)) to
    results/<name>.json, creating the directory if needed."""
    os.makedirs(RESULTS_DIR, exist_ok=True)
    payload = {}
    for key, val in values.items():
        if isinstance(val, tuple):
            number, fmt = val
            payload[key] = {"value": _to_builtin(number), "fmt": fmt}
        else:
            payload[key] = {"value": _to_builtin(val)}
    path = os.path.join(RESULTS_DIR, f"{name}.json")
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
    return path


def read_results(name):
    """Read results/<name>.json back as {key: value} (formats dropped).
    Lets one script use another's output (e.g. the fitted omega_dot
    precision) instead of hardcoding a number that can go stale."""
    path = os.path.join(RESULTS_DIR, f"{name}.json")
    with open(path, encoding="utf-8") as handle:
        payload = json.load(handle)
    return {key: entry["value"] for key, entry in payload.items()}


def _to_builtin(val):
    """Coerce numpy scalars to JSON-serializable Python builtins."""
    if hasattr(val, "item"):
        return val.item()
    return val
