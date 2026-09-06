# SPDX-FileCopyrightText: 2026 Alberto P
#
# SPDX-License-Identifier: MPL-2.0

import argparse
import json
from .mie import compute_mie_scattering, DEFAULT_EXECUTABLE


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Compute Mie scattering using the MIEV0 routine.")
    parser.add_argument("--executable", default=DEFAULT_EXECUTABLE)
    parser.add_argument("--m_real", type=float, required=True)
    parser.add_argument("--m_img", type=float, required=True)
    parser.add_argument("--radius", type=float, required=True)
    parser.add_argument("--wavelength", type=float, required=True)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    mie_scatt_params = compute_mie_scattering(
        args.m_real,
        args.m_img,
        args.radius,
        args.wavelength,
        args.executable
    )
    print(json.dumps(mie_scatt_params))
