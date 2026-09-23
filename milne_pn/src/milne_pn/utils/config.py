"""
config.py
---------
Thin YAML config loader for config/config.yaml, with sane defaults so the
package works even if no config file is supplied.
"""

import os

DEFAULTS = {
    "solver": {
        "n_max": 79,
        "n_plot": 39,
        "zero_tol": 1e-9,
        "quad_multiplier": 4,
    },
    "runge_analysis": {
        "orders": [3, 9, 21, 41],
        "n_mu": 2000,
    },
    "paths": {
        "data_output": "data/output",
        "data_input": "data/input",
        "plots": "data/output/plots",
    },
    "reference": {
        "z0_exact": 0.710446,
    },
}


def _deep_merge(base, override):
    out = dict(base)
    for k, v in override.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = v
    return out


def load_config(path=None):
    if path is None:
        here = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(here, "..", "..", "..", "config", "config.yaml")
    path = os.path.abspath(path)

    if not os.path.exists(path):
        return dict(DEFAULTS)

    try:
        import yaml
    except ImportError:
        return dict(DEFAULTS)

    with open(path) as f:
        user_cfg = yaml.safe_load(f) or {}
    return _deep_merge(DEFAULTS, user_cfg)
