"""Shared argument parser for example plot scripts."""

import argparse
import pathlib


def create_argument_parser(
    description: str, caller_file: str
) -> argparse.ArgumentParser:
    """Return an argument parser with common ``--dir`` and ``--format`` options.

    Parameters
    ----------
    description:
        Short description shown in ``--help`` (typically ``__doc__``).
    caller_file:
        ``__file__`` of the calling script, used to resolve the default
        ``--dir`` relative to the script location.
    """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "--dir",
        type=pathlib.Path,
        default=pathlib.Path(caller_file).resolve().parent / "build",
        help="Directory containing CSV output files (default: <script_dir>/build)",
    )
    parser.add_argument(
        "--format",
        choices=["pdf", "png"],
        default="png",
        help="Output image format (default: png)",
    )
    return parser
