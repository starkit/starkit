# -*- coding: utf-8 -*-
from __future__ import print_function
import os
import re
import json

def parse_model_filename(filename):
    """
    Parses a stellar model filename into its 16 component parameters.

    Supports filenames with variations in case (e.g., Vinfty/Vinfty, finfty1.0, etc.)
    Ignores file extensions at the end.
    """
    filename = os.path.basename(filename)
    # Remove file extension if present
    filename = os.path.splitext(filename)[0]

    pattern = re.compile(
        r"^(?P<model_family>[A-Z])_"
        r"(?P<model_type>[a-z])"
        r"(?P<t_eff>\d+)"
        r"g(?P<gravity>[\d.]+)"
        r"z(?P<metallicity>[\d.-]+)"
        r"(t|f)(?P<microturbulent_velocity>[\d.-]+)_"
        r"a(?P<alpha_fe>[\d.-]+)"
        r"c(?P<c_fe>[\d.-]+)"
        r"n(?P<n_fe>[\d.-]+)"
        r"o(?P<o_fe>[\d.-]+)_"
        r"Mdot(?P<mass_loss>[\d.-]+)"
        r"(V|v)infty(?P<terminal_velocity>[\d]+)"
        r"beta(?P<beta>[\d.]+)"
        r"finfty(?P<first_clumping_parameter>[\d.]+)"
        r"vcl(?P<second_clumping_parameter>[\d.]+)_"
        r"(?P<spectral_domain>[A-Z]+)$"
    )

    match = pattern.match(filename)

    if not match:
        return None

    params = match.groupdict()

    type_conversions = {
        't_eff': int,
        'gravity': float,
        'metallicity': float,
        'microturbulent_velocity': float,
        'alpha_fe': float,
        'c_fe': float,
        'n_fe': float,
        'o_fe': float,
        'mass_loss': float,
        'terminal_velocity': int,
        'beta': float,
        'first_clumping_parameter': float,
        'second_clumping_parameter': float,
    }

    for key, type_func in type_conversions.items():
        if key in params:
            params[key] = type_func(params[key])

    return params

# --- Example Usage ---
if __name__ == "__main__":
    example_filename = "C_s35000g3.40z0.0f5.0_a0.00c0.00n0.00o0.00_Mdot-5.41vinfty1760beta0.8finfty1vcl0_VIS"
    
    parsed_data = parse_model_filename(example_filename)

    if parsed_data:
        # Replaced f-strings with the .format() method for Python 2.7 compatibility
        print(u"✅ Successfully parsed filename: '{}'".format(example_filename))
        print(json.dumps(parsed_data, indent=4))
        print(u"\nTotal variables parsed: {}".format(len(parsed_data)))
    else:
        print(u"❌ Failed to parse filename: '{}'".format(example_filename))