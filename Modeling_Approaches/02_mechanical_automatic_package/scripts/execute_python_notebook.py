"""Execute Python code cells from a notebook without requiring Jupyter."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path


def execute_notebook(path: Path) -> None:
    notebook = json.loads(path.read_text(encoding="utf-8"))
    source = "\n\n".join(
        cell["source"] if isinstance(cell["source"], str) else "".join(cell["source"])
        for cell in notebook["cells"]
        if cell["cell_type"] == "code"
    )
    os.chdir(path.parent)
    exec(compile(source, str(path), "exec"), {"__name__": "__main__"})


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("notebook", type=Path)
    args = parser.parse_args()
    execute_notebook(args.notebook.resolve())
