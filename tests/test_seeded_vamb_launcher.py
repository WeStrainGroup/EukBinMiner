from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
LAUNCHER = PROJECT / "scripts" / "run_seeded_vamb.py"
SPEC = importlib.util.spec_from_file_location("run_seeded_vamb", LAUNCHER)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class SeededVambLauncherTests(unittest.TestCase):
    def test_argument_contract(self) -> None:
        seed, command = MODULE.parse_args(
            ["--seed", "6311", "--", "/example/vamb", "bin", "default"]
        )
        self.assertEqual(seed, 6311)
        self.assertEqual(command, ["/example/vamb", "bin", "default"])
        with self.assertRaisesRegex(ValueError, "usage"):
            MODULE.parse_args(["--seed", "6311", "/example/vamb"])

    def test_absolute_console_interpreter_is_discovered(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            executable = Path(temporary) / "vamb"
            executable.write_text(
                f"#!{Path(sys.executable)}\n",
                encoding="utf-8",
            )
            executable.chmod(0o755)
            interpreter = MODULE.console_interpreter(executable)
            self.assertEqual(interpreter, Path(sys.executable))


if __name__ == "__main__":
    unittest.main()
