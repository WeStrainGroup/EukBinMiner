#!/usr/bin/env python3
"""Launch the exact VAMB console script after seeding process RNGs.

VAMB 5.0.4 records ``--seed`` but constructs the Taxometer ``VAMB2Label``
network before applying that seed.  Seeding Python, NumPy, and PyTorch before
the console script is imported makes the advertised seed cover Taxometer too.
"""

from __future__ import annotations

import os
import random
import runpy
import sys
from pathlib import Path
from typing import Optional, Sequence


REEXEC_MARKER = "EUKBINMINER_VAMB_SEED_REEXEC"


def parse_args(argv: Sequence[str]) -> tuple[int, list[str]]:
    if len(argv) < 4 or argv[0] != "--seed" or argv[2] != "--":
        raise ValueError(
            "usage: run_seeded_vamb.py --seed INT -- VAMB_EXE [VAMB_ARGS...]"
        )
    try:
        seed = int(argv[1])
    except ValueError as error:
        raise ValueError("VAMB seed must be an integer") from error
    command = list(argv[3:])
    if not command:
        raise ValueError("VAMB command is empty")
    return seed, command


def console_interpreter(executable: Path) -> Optional[Path]:
    try:
        with executable.open("rb") as handle:
            first_line = handle.readline(4096)
    except OSError as error:
        raise ValueError(f"cannot read VAMB executable: {executable}") from error
    if not first_line.startswith(b"#!"):
        return None
    shebang = first_line[2:].decode("utf-8", errors="strict").strip()
    if not shebang.startswith("/") or any(character.isspace() for character in shebang):
        return None
    interpreter = Path(shebang)
    if not interpreter.is_file() or not os.access(interpreter, os.X_OK):
        raise ValueError(f"VAMB interpreter is not executable: {interpreter}")
    return interpreter


def reexec_with_console_python(
    seed: int,
    command: Sequence[str],
) -> None:
    interpreter = console_interpreter(Path(command[0]))
    if interpreter is None or os.environ.get(REEXEC_MARKER) == "1":
        return
    if interpreter.resolve() == Path(sys.executable).resolve():
        return
    environment = dict(os.environ)
    environment[REEXEC_MARKER] = "1"
    os.execve(
        interpreter,
        [
            str(interpreter),
            str(Path(__file__).resolve()),
            "--seed",
            str(seed),
            "--",
            *command,
        ],
        environment,
    )


def seed_process(seed: int) -> None:
    import numpy
    import torch

    random.seed(seed)
    numpy.random.seed(seed % (2**32))
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def main(argv: Optional[Sequence[str]] = None) -> int:
    arguments = list(sys.argv[1:] if argv is None else argv)
    seed, command = parse_args(arguments)
    executable = Path(command[0]).expanduser().resolve()
    if not executable.is_file():
        raise ValueError(f"VAMB executable does not exist: {executable}")
    command[0] = str(executable)
    reexec_with_console_python(seed, command)
    seed_process(seed)
    os.environ.pop(REEXEC_MARKER, None)
    sys.argv = command
    runpy.run_path(str(executable), run_name="__main__")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        print(f"run_seeded_vamb.py: error: {error}", file=sys.stderr)
        raise SystemExit(2)
