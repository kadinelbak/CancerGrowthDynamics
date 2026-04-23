#!/usr/bin/env python3
from __future__ import annotations

import ast
import re
from pathlib import Path
from typing import Iterable

import pandas as pd


METRIC_COLUMNS = {"sse", "aic", "bic", "loss", "ssr", "converged"}
META_COLUMNS = {
    "file",
    "density",
    "mix",
    "mix_label",
    "cellline",
    "cell_line",
    "model",
    "title",
    "params",
}


def find_repo_root(start: Path) -> Path:
    current = start.resolve()
    for _ in range(8):
        if (current / "Modeling_Approaches/01_mechanical_manual").exists() and (current / "Processed_Datasets").exists():
            return current
        if current.parent == current:
            break
        current = current.parent
    raise FileNotFoundError("Could not find repository root containing Modeling_Approaches/01_mechanical_manual and Processed_Datasets")


def to_float(value: object) -> float | None:
    try:
        if value is None:
            return None
        text = str(value).strip()
        if text == "":
            return None
        return float(text)
    except Exception:
        return None


def parse_params_text(params_text: str) -> list[tuple[str, float | str]]:
    text = (params_text or "").strip()
    if not text:
        return []

    if text.startswith("[") and text.endswith("]"):
        try:
            raw = ast.literal_eval(text)
            if isinstance(raw, (list, tuple)):
                parsed: list[tuple[str, float | str]] = []
                for idx, value in enumerate(raw, start=1):
                    numeric = to_float(value)
                    parsed.append((f"p{idx}", numeric if numeric is not None else str(value)))
                return parsed
        except Exception:
            pass

    if "=" in text:
        pieces = [chunk.strip() for chunk in text.split(",") if chunk.strip()]
        parsed = []
        for piece in pieces:
            if "=" not in piece:
                continue
            key, val = piece.split("=", 1)
            key = key.strip()
            val = val.strip()
            numeric = to_float(val)
            parsed.append((key, numeric if numeric is not None else val))
        if parsed:
            return parsed

    return [("params_raw", text)]


def normalize_columns(frame: pd.DataFrame) -> pd.DataFrame:
    renamed = {}
    for col in frame.columns:
        cleaned = re.sub(r"\s+", "_", str(col).strip().lower())
        renamed[col] = cleaned
    return frame.rename(columns=renamed)


def extract_parameter_rows(source_group: str, file_path: Path) -> list[dict[str, object]]:
    df = pd.read_csv(file_path)
    df = normalize_columns(df)
    rows: list[dict[str, object]] = []

    for _, row in df.iterrows():
        density = row.get("density")
        mix_label = row.get("mix_label", row.get("mix"))
        cell_line = row.get("cell_line", row.get("cellline"))
        model = row.get("model", row.get("title"))

        metrics = {
            "sse": row.get("sse", row.get("ssr")),
            "aic": row.get("aic"),
            "bic": row.get("bic"),
            "loss": row.get("loss"),
        }

        param_entries: list[tuple[str, float | str]] = []

        if "params" in row and pd.notna(row.get("params")):
            param_entries.extend(parse_params_text(str(row.get("params"))))

        for col in df.columns:
            if col in META_COLUMNS or col in METRIC_COLUMNS:
                continue
            value = row.get(col)
            if pd.isna(value):
                continue
            numeric = to_float(value)
            if numeric is None:
                continue
            param_entries.append((col, numeric))

        for parameter_name, parameter_value in param_entries:
            rows.append(
                {
                    "source_group": source_group,
                    "source_file": str(file_path).replace("\\", "/"),
                    "density": density,
                    "mix_label": mix_label,
                    "cell_line": cell_line,
                    "model": model,
                    "parameter_name": parameter_name,
                    "parameter_value": parameter_value,
                    "sse": metrics["sse"],
                    "aic": metrics["aic"],
                    "bic": metrics["bic"],
                    "loss": metrics["loss"],
                }
            )

    return rows


def iter_parameter_sources(root: Path) -> Iterable[tuple[str, Path]]:
    base = root / "Modeling_Approaches/01_mechanical_manual"
    expected = [
        ("simplemodels", base / "SimpleModels" / "coculture_option1_params.csv"),
        ("simplemodels", base / "SimpleModels" / "coculture_option2_params.csv"),
        ("simplemodels", base / "SimpleModels" / "coculture_fit_params.csv"),
        ("untreated_monoculture", base / "Untreated MonoCulture" / "untreated_monoculture_logistic_parameters.csv"),
        ("untreated_monoculture", base / "Untreated MonoCulture" / "untreated_monoculture_all_models_parameters.csv"),
        ("treated_monoculture", base / "Treated MonoCulture" / "outputs" / "treated_monoculture_best_models.csv"),
        ("treated_monoculture", base / "Treated MonoCulture" / "outputs" / "treated_monoculture_model_bic_results.csv"),
        ("untreated_coculture", base / "Untreated CoCulture" / "untreated_coculture_simplified_fits.csv"),
    ]
    for source_group, path in expected:
        if path.exists():
            yield source_group, path


def main() -> int:
    root = find_repo_root(Path(__file__).resolve().parent)
    out_file = root / "Processed_Datasets" / "standardized_parameters.csv"
    out_file.parent.mkdir(parents=True, exist_ok=True)

    all_rows: list[dict[str, object]] = []
    for source_group, csv_path in iter_parameter_sources(root):
        extracted = extract_parameter_rows(source_group, csv_path)
        all_rows.extend(extracted)
        print(f"Loaded {len(extracted):4d} rows from {csv_path}")

    if not all_rows:
        print("No parameter rows extracted. Nothing written.")
        return 1

    out_df = pd.DataFrame(all_rows)
    for col in ("density", "mix_label", "cell_line", "model", "parameter_name"):
        out_df[col] = out_df[col].fillna("")

    out_df = out_df.sort_values(
        by=["source_group", "density", "mix_label", "cell_line", "model", "parameter_name"],
        kind="stable",
    ).reset_index(drop=True)

    out_df.to_csv(out_file, index=False, float_format="%.8g")
    print(f"Wrote {len(out_df)} standardized parameter rows -> {out_file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())