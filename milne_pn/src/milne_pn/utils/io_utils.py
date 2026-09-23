"""
io_utils.py
-----------
Minimal CSV read/write helpers (dict-of-columns <-> csv), kept dependency-
free (no pandas) since numpy + the stdlib csv module cover everything this
project needs.
"""

import csv
import os
import numpy as np


def write_csv(path, columns, fmt="{:.10g}"):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    headers = list(columns.keys())
    n = len(next(iter(columns.values())))
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        for i in range(n):
            writer.writerow([fmt.format(columns[h][i]) for h in headers])


def read_csv(path):
    with open(path, newline="") as f:
        reader = csv.reader(f)
        header = [h.strip() for h in next(reader)]
        rows = [row for row in reader if row]
    return {name: np.array([float(r[i]) for r in rows]) for i, name in enumerate(header)}
