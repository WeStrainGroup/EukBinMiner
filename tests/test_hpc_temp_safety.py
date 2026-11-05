from __future__ import annotations

import ast
import re
import unittest
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
PRODUCTION_ROOTS = (PROJECT / "bin", PROJECT / "scripts")
PYTHON_SOURCES = (
    PROJECT / "bin" / "eukbinminer",
    PROJECT / "bin" / "eukmagminer",
    *sorted((PROJECT / "scripts").glob("*.py")),
)


class HpcTemporaryStorageSafetyTests(unittest.TestCase):
    def test_python_tempfile_calls_have_explicit_parent(self) -> None:
        unsafe: list[str] = []
        tempfile_calls = {"mkdtemp", "mkstemp", "TemporaryDirectory"}
        for path in PYTHON_SOURCES:
            tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
            for node in ast.walk(tree):
                if not isinstance(node, ast.Call):
                    continue
                function = node.func
                if not (
                    isinstance(function, ast.Attribute)
                    and isinstance(function.value, ast.Name)
                    and function.value.id == "tempfile"
                    and function.attr in tempfile_calls
                ):
                    continue
                directory = next(
                    (
                        keyword.value
                        for keyword in node.keywords
                        if keyword.arg == "dir"
                    ),
                    None,
                )
                if directory is None or (
                    isinstance(directory, ast.Constant)
                    and directory.value is None
                ):
                    unsafe.append(f"{path}:{node.lineno}")
        self.assertEqual(unsafe, [])

    def test_shell_mktemp_calls_use_explicit_templates(self) -> None:
        unsafe: list[str] = []
        for path in sorted((PROJECT / "scripts").glob("*.sh")):
            lines = path.read_text(encoding="utf-8").splitlines()
            for index, line in enumerate(lines):
                if not re.search(r"\bmktemp\b", line):
                    continue
                context = " ".join(lines[index : index + 3])
                if "XXXXXX" not in context or "${" not in context:
                    unsafe.append(f"{path}:{index + 1}")
        self.assertEqual(unsafe, [])

    def test_production_code_has_no_node_root_temp_literal(self) -> None:
        unsafe: list[str] = []
        for root in PRODUCTION_ROOTS:
            for path in sorted(root.glob("*")):
                if not path.is_file():
                    continue
                for line_number, line in enumerate(
                    path.read_text(encoding="utf-8").splitlines(),
                    start=1,
                ):
                    if '"/tmp' in line or "'/tmp" in line:
                        unsafe.append(f"{path}:{line_number}")
                    if '"/var/tmp' in line or "'/var/tmp" in line:
                        unsafe.append(f"{path}:{line_number}")
        self.assertEqual(unsafe, [])


    def test_taxvamb_environment_uses_kunpeng_not_kraken(self) -> None:
        environment = (
            PROJECT / "envs" / "vamb_taxvamb.yml"
        ).read_text(
            encoding="utf-8"
        )
        self.assertIn("kun-peng=0.7.12", environment)
        self.assertNotIn("kraken2=", environment)


if __name__ == "__main__":
    unittest.main()
