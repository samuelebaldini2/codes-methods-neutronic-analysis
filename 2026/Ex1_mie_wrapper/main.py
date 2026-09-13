import time
import json
import miescat

# ============================================================
# REFERENCE DATA
# ============================================================

m_real = 1.02239
m_img = -0.0119
radius = 2.0
wavelength = 0.5

n_calls_list = [1, 5, 10, 20]


# ============================================================
# SINGLE RUN
# ============================================================

print("==================================================================")
print("compute_mie_scattering")
print("==================================================================")

result = miescat.compute_mie_scattering(
    m_real,
    m_img,
    radius,
    wavelength,
)

print("\nRESULTS")
print("------------------------------------------------------------------")

for key, value in result.items():
    print(f"{key:<35} : {value}")


# ============================================================
# SAVE NUMERICAL RESULT
# ============================================================

with open("mie_result.json", "w") as f:
    json.dump(result, f, indent=4)


# ============================================================
# PERFORMANCE TEST
# ============================================================

print("\n")
print("==================================================================")
print("PERFORMANCE TEST")
print("==================================================================")

print(
    f"{'N calls':>10} "
    f"{'Total time [s]':>20} "
    f"{'Time/call [s]':>20} "
    f"{'Calls/s':>15}"
)

print("-" * 70)


timing_results = []


for n_calls in n_calls_list:

    start = time.perf_counter()

    for _ in range(n_calls):

        miescat.compute_mie_scattering(
            m_real,
            m_img,
            radius,
            wavelength,
        )

    end = time.perf_counter()

    total_time = end - start

    time_per_call = total_time / n_calls

    calls_per_second = n_calls / total_time


    print(
        f"{n_calls:10d} "
        f"{total_time:20.8e} "
        f"{time_per_call:20.8e} "
        f"{calls_per_second:15.3f}"
    )


    timing_results.append(
        {
            "n_calls": n_calls,
            "total_time_s": total_time,
            "time_per_call_s": time_per_call,
            "calls_per_second": calls_per_second,
        }
    )


# ============================================================
# SAVE TIMING RESULTS
# ============================================================

with open("timing_results.json", "w") as f:
    json.dump(
        timing_results,
        f,
        indent=4,
    )


print("\n")
print("==================================================================")
print("OUTPUT FILES")
print("==================================================================")

print("Numerical results : mie_result.json")
print("Timing results    : timing_results.json")
