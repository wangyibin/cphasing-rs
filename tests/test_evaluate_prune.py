import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "examples" / "kprune" / "evaluate_prune.py"


class EvaluatePruneTests(unittest.TestCase):
    def test_canonicalizes_pairs_and_ignores_unrelated_predictions(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            contigs = directory / "contigsizes"
            contacts = directory / "contacts"
            prune_table = directory / "prune.table"
            contigs.write_text(
                "1A.a\t100\n"
                "1A.b\t100\n"
                "1B.c\t100\n"
                "1B.d\t100\n"
                "2A.e\t100\n",
                encoding="utf-8",
            )
            contacts.write_text(
                "1A.a\t1B.c\t1\n"
                "1B.d\t1A.a\t0.5\n"
                "1A.a\t1A.b\t1\n"
                "1A.a\t2A.e\t1\n",
                encoding="utf-8",
            )
            prune_table.write_text(
                "1B.c\t1A.a\t0\t0\t0\t0\t1\n"
                "1A.a\t1A.b\t0\t0\t0\t0\t1\n"
                "1A.a\t2A.e\t0\t0\t0\t0\t1\n",
                encoding="utf-8",
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    str(contigs),
                    str(prune_table),
                    "--contacts",
                    str(contacts),
                ],
                check=True,
                capture_output=True,
                text=True,
            )

        self.assertEqual(result.stdout, "1A.a\t1B.d\t0.5\n")
        self.assertIn("Precision:0.5\n", result.stderr)
        self.assertIn("Recall:0.5\n", result.stderr)
        self.assertIn("F1 score:0.5\n", result.stderr)


if __name__ == "__main__":
    unittest.main()
