#!/usr/bin/env python3

# SPDX-FileCopyrightText: 2026 Alberto P
#
# SPDX-License-Identifier: MPL-2.0

import subprocess
from pathlib import Path
import argparse
import json


def run_console(program, workdir, *args):
    console_result = subprocess.run(
        [program] + list(args),
        capture_output=True,
        text=True,
        check=True,
        cwd=workdir
    )

    return console_result.stdout, console_result.stderr, console_result.returncode


def compute_mie_scattering(executable, m_real, m_img, radius, wavelength):
    output, _, _ = run_console(executable,
                               Path.cwd(),
                               str(m_real),
                               str(m_img),
                               str(radius),
                               str(wavelength)
                            )
    mie_scatt_params = json.loads(output)
    return mie_scatt_params


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", default=Path(__file__).parent / "miescat")
    parser.add_argument("--m_real", type=float, required=True)
    parser.add_argument("--m_img", type=float, required=True)
    parser.add_argument("--radius", type=float, required=True)
    parser.add_argument("--wavelength", type=float, required=True)

    args = parser.parse_args()

    mie_scatt_params = compute_mie_scattering(
        args.executable,
        args.m_real,
        args.m_img,
        args.radius,
        args.wavelength
    )
    print(json.dumps(mie_scatt_params))


if __name__ == "__main__":
    main()
