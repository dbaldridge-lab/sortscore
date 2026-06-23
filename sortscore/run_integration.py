"""
Main entry point for integration-specific export helpers.
"""

import argparse
import logging
import sys
from pathlib import Path

from sortscore.integrations.lilace import score_table_to_lilace_csv


def _build_lilace_input(args: argparse.Namespace) -> Path:
    metadata_columns = None
    metadata_constants = None
    if args.score_column:
        metadata_columns = {"sortscore_score": args.score_column}
        metadata_constants = {"sortscore_score_column": args.score_column}

    return score_table_to_lilace_csv(
        args.input,
        args.output,
        mutagenesis_type=args.mutagenesis_type,
        batch=args.batch,
        metadata_columns=metadata_columns,
        metadata_constants=metadata_constants,
    )


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s:%(name)s:%(message)s")

    parser = argparse.ArgumentParser(description="Export integration-specific files from sortscore outputs.")
    subparsers = parser.add_subparsers(dest="integration", required=True)

    lilace_parser = subparsers.add_parser(
        "lilace",
        help="Convert a sortscore score table into a Lilace input CSV.",
    )
    lilace_parser.add_argument("-i", "--input", required=True, help="Path to a sortscore score table CSV")
    lilace_parser.add_argument("-o", "--output", required=True, help="Path for the Lilace input CSV")
    lilace_parser.add_argument(
        "--mutagenesis-type",
        choices=["aa", "codon"],
        default="aa",
        help="Mutagenesis type to export",
    )
    lilace_parser.add_argument(
        "-b",
        "--batch",
        help="Optional batch/tile label to export, e.g. tile4",
    )
    lilace_parser.add_argument(
        "--score-column",
        default="avgscore",
        help="Optional sortscore score column to preserve as metadata",
    )
    lilace_parser.set_defaults(handler=_build_lilace_input)

    args = parser.parse_args()

    try:
        output_path = args.handler(args)
    except Exception as exc:
        logging.error("Failed to export integration file: %s", exc)
        sys.exit(1)

    print(f"Integration file written to {output_path}")


if __name__ == "__main__":
    main()
