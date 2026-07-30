"""Command-line interface for Starkit analysis workflows."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="starkit-analysis",
        description="Run and validate Starkit stellar-analysis workflows.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run", help="Run an analysis from a JSON configuration.")
    run_parser.add_argument("config", help="Path to the analysis JSON configuration file.")
    run_parser.set_defaults(handler=_run)

    validate_parser = subparsers.add_parser(
        "validate",
        help="Validate an analysis configuration before running it.",
    )
    validate_parser.add_argument("config", help="Path to the analysis JSON configuration file.")
    validate_parser.add_argument(
        "--shallow",
        action="store_true",
        help="Skip the deep preflight that loads the grid, spectrum, and fit model.",
    )
    validate_parser.add_argument(
        "--json",
        action="store_true",
        help="Emit a machine-readable JSON report.",
    )
    validate_parser.set_defaults(handler=_validate)

    inspect_parser = subparsers.add_parser(
        "inspect-grid",
        help="Inspect a Starkit grid and optionally generate a starter configuration.",
    )
    inspect_parser.add_argument("grid", help="Path to the Starkit grid HDF5 file.")
    inspect_parser.add_argument(
        "--template-out",
        help="Write a starter analysis JSON configuration to this path.",
    )
    inspect_parser.add_argument(
        "--summary-json-out",
        help="Write machine-readable grid inspection JSON to this path.",
    )
    inspect_parser.add_argument("--name", help="Override run.name in the generated configuration.")
    inspect_parser.add_argument(
        "--output-directory",
        help="Override outputs.directory in the generated configuration.",
    )
    inspect_parser.set_defaults(handler=_inspect_grid)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.handler(args)


def _run(args: argparse.Namespace) -> int:
    from .runner import run_analysis

    return run_analysis(args.config)


def _validate(args: argparse.Namespace) -> int:
    from .validation import format_report, verify_config

    report = verify_config(args.config, deep=not args.shallow)
    if args.json:
        print(json.dumps(report.to_dict(), indent=2))
    else:
        print(format_report(report))
    return 0 if report.ok else 1


def _inspect_grid(args: argparse.Namespace) -> int:
    from starkit.gridkit.inspector import (
        format_inspection_summary,
        inspect_grid,
        inspection_to_dict,
    )

    from .templates import build_template_config

    inspection = inspect_grid(args.grid)
    print(format_inspection_summary(inspection))

    if args.summary_json_out:
        summary_path = Path(args.summary_json_out).expanduser()
        summary_path.write_text(
            json.dumps(inspection_to_dict(inspection), indent=2, sort_keys=True) + "\n",
            encoding="ascii",
        )
        print(f"\nWrote inspection JSON to {summary_path}")

    if args.template_out:
        template = build_template_config(
            inspection,
            run_name=args.name,
            output_directory=args.output_directory,
        )
        template_path = Path(args.template_out).expanduser()
        template_path.write_text(
            json.dumps(template, indent=2, sort_keys=False) + "\n",
            encoding="ascii",
        )
        print(f"Wrote starter config to {template_path}")

    return 0
