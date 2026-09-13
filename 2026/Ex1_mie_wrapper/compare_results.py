import json
import csv
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt


BASE_DIR = Path(__file__).resolve().parent

OLD_RESULT_FILE = BASE_DIR / "mie_result_old.json"
NEW_RESULT_FILE = BASE_DIR / "mie_result_new.json"

OLD_TIMING_FILE = BASE_DIR / "timing_results_old.json"
NEW_TIMING_FILE = BASE_DIR / "timing_results_new.json"

NUMERICAL_CSV = BASE_DIR / "numerical_comparison.csv"
TIMING_CSV = BASE_DIR / "timing_comparison.csv"
TIMING_PNG = BASE_DIR / "timing_comparison.png"


# ============================================================
# UTILITY
# ============================================================

def load_json(path):
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    with open(path, "r") as f:
        return json.load(f)


# ============================================================
# LOAD DATA
# ============================================================

old_result = load_json(OLD_RESULT_FILE)
new_result = load_json(NEW_RESULT_FILE)

old_timing = load_json(OLD_TIMING_FILE)
new_timing = load_json(NEW_TIMING_FILE)


# ============================================================
# NUMERICAL COMPARISON
# ============================================================

print("=" * 80)
print("NUMERICAL COMPARISON")
print("=" * 80)

numerical_rows = []


for key in old_result:

    if key not in new_result:
        print(f"{key}: missing in NEW results")
        continue

    old_value = np.asarray(
        old_result[key],
        dtype=float
    )

    new_value = np.asarray(
        new_result[key],
        dtype=float
    )

    if old_value.shape != new_value.shape:

        print(
            f"{key}: shape mismatch "
            f"{old_value.shape} != "
            f"{new_value.shape}"
        )

        continue


    old_flat = old_value.ravel()
    new_flat = new_value.ravel()


    # ========================================================
    # ABSOLUTE ERROR
    # ========================================================

    abs_error = np.abs(
        new_flat - old_flat
    )


    # ========================================================
    # RELATIVE ERROR
    # ========================================================

    rel_error = np.full(
        abs_error.shape,
        np.nan,
        dtype=float
    )

    nonzero = np.abs(old_flat) > 0

    rel_error[nonzero] = (
        abs_error[nonzero]
        / np.abs(old_flat[nonzero])
    )


    max_abs_error = np.max(
        abs_error
    )

    finite_rel = rel_error[
        np.isfinite(rel_error)
    ]

    if len(finite_rel) > 0:

        max_rel_error = np.max(
            finite_rel
        )

    else:

        max_rel_error = np.nan


    # ========================================================
    # LEGENDRE MOMENTS
    # ========================================================

    if key == "legendre_moments":

        mean_abs_error = np.mean(
            abs_error
        )

        if len(finite_rel) > 0:

            mean_rel_error = np.mean(
                finite_rel
            )

        else:

            mean_rel_error = np.nan


        numerical_rows.append(
            {
                "quantity": key,
                "absolute_error":
                    "",
                "relative_error":
                    "",
                "max_absolute_error":
                    max_abs_error,
                "mean_absolute_error":
                    mean_abs_error,
                "max_relative_error":
                    max_rel_error,
                "mean_relative_error":
                    mean_rel_error,
            }
        )


        print(f"\n{key}")

        print(
            f"  max absolute error  : "
            f"{max_abs_error:.6e}"
        )

        print(
            f"  mean absolute error : "
            f"{mean_abs_error:.6e}"
        )

        print(
            f"  max relative error  : "
            f"{max_rel_error:.6e}"
        )

        print(
            f"  mean relative error : "
            f"{mean_rel_error:.6e}"
        )


    # ========================================================
    # SCALAR QUANTITIES
    # ========================================================

    else:

        abs_error_scalar = float(
            max_abs_error
        )

        rel_error_scalar = float(
            max_rel_error
        )


        numerical_rows.append(
            {
                "quantity": key,
                "absolute_error":
                    abs_error_scalar,
                "relative_error":
                    rel_error_scalar,
                "max_absolute_error":
                    "",
                "mean_absolute_error":
                    "",
                "max_relative_error":
                    "",
                "mean_relative_error":
                    "",
            }
        )


        print(f"\n{key}")

        print(
            f"  absolute error : "
            f"{abs_error_scalar:.6e}"
        )

        print(
            f"  relative error : "
            f"{rel_error_scalar:.6e}"
        )


# ============================================================
# SAVE NUMERICAL COMPARISON
# ============================================================

with open(
    NUMERICAL_CSV,
    "w",
    newline=""
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=[
            "quantity",
            "absolute_error",
            "relative_error",
            "max_absolute_error",
            "mean_absolute_error",
            "max_relative_error",
            "mean_relative_error",
        ],
    )

    writer.writeheader()

    writer.writerows(
        numerical_rows
    )


# ============================================================
# TIMING COMPARISON
# ============================================================

print("\n")
print("=" * 80)
print("TIMING COMPARISON")
print("=" * 80)


old_by_n = {
    int(row["n_calls"]): row
    for row in old_timing
}


new_by_n = {
    int(row["n_calls"]): row
    for row in new_timing
}


common_n = sorted(
    set(old_by_n.keys())
    &
    set(new_by_n.keys())
)


if not common_n:

    raise RuntimeError(
        "No common n_calls values found "
        "between old and new timing files."
    )


timing_rows = []


print(
    f"\n{'N':>5} "
    f"{'OLD total [s]':>16} "
    f"{'NEW total [s]':>16} "
    f"{'OLD s/call':>16} "
    f"{'NEW s/call':>16} "
    f"{'Speedup':>10}"
)

print("-" * 90)


for n in common_n:

    old_row = old_by_n[n]
    new_row = new_by_n[n]


    old_total = float(
        old_row["total_time_s"]
    )

    new_total = float(
        new_row["total_time_s"]
    )


    old_per_call = float(
        old_row["time_per_call_s"]
    )

    new_per_call = float(
        new_row["time_per_call_s"]
    )


    old_calls_per_s = float(
        old_row["calls_per_second"]
    )

    new_calls_per_s = float(
        new_row["calls_per_second"]
    )


    speedup = (
        old_total
        /
        new_total
    )


    timing_rows.append(
        {
            "n_calls":
                n,

            "old_total_time_s":
                old_total,

            "new_total_time_s":
                new_total,

            "old_time_per_call_s":
                old_per_call,

            "new_time_per_call_s":
                new_per_call,

            "old_calls_per_second":
                old_calls_per_s,

            "new_calls_per_second":
                new_calls_per_s,

            "speedup":
                speedup,
        }
    )


    print(
        f"{n:5d} "
        f"{old_total:16.6e} "
        f"{new_total:16.6e} "
        f"{old_per_call:16.6e} "
        f"{new_per_call:16.6e} "
        f"{speedup:10.2f}"
    )


# ============================================================
# SAVE TIMING TABLE
# ============================================================

with open(
    TIMING_CSV,
    "w",
    newline=""
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=[
            "n_calls",
            "old_total_time_s",
            "new_total_time_s",
            "old_time_per_call_s",
            "new_time_per_call_s",
            "old_calls_per_second",
            "new_calls_per_second",
            "speedup",
        ],
    )

    writer.writeheader()

    writer.writerows(
        timing_rows
    )


# ============================================================
# TIMING PLOT
# ============================================================

n_values = np.array(
    [
        row["n_calls"]
        for row in timing_rows
    ]
)

old_times = np.array(
    [
        row["old_total_time_s"]
        for row in timing_rows
    ]
)

new_times = np.array(
    [
        row["new_total_time_s"]
        for row in timing_rows
    ]
)


plt.figure(
    figsize=(8, 5)
)


plt.plot(
    n_values,
    old_times,
    marker="o",
    label="Old wrapper"
)


plt.plot(
    n_values,
    new_times,
    marker="o",
    label="New wrapper"
)


plt.xlabel(
    "Number of successive calls"
)

plt.ylabel(
    "Total execution time [s]"
)

plt.title(
    "Mie wrapper performance comparison"
)

plt.grid(
    True,
    alpha=0.3
)

plt.legend()

plt.tight_layout()


plt.savefig(
    TIMING_PNG,
    dpi=200
)


plt.close()


# ============================================================
# FINAL SUMMARY
# ============================================================

print("\n")
print("=" * 80)
print("FILES GENERATED")
print("=" * 80)

print(NUMERICAL_CSV)
print(TIMING_CSV)
print(TIMING_PNG)