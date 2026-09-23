#!/usr/bin/env python3
"""
scripts/run_gui.py
--------------------
Launch the Milne PN Tkinter desktop GUI.

Usage:
    python scripts/run_gui.py

Equivalent to `python -m milne_pn.gui` or, after `pip install -e .`, the
`milne-gui` console script.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

from milne_pn.gui.app import main

if __name__ == "__main__":
    main()
