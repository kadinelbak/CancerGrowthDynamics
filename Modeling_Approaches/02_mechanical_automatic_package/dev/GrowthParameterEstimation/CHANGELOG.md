# Changelog

All notable changes to this project will be documented in this file.

## v0.3.0 - 2026-04-09

### Breaking changes

- Breaking: this is a pre-1.0 minor release, so `v0.3.0` is treated as a breaking update relative to `v0.2.x` for Julia package registration and downstream compatibility expectations.
- Breaking: staged and workflow-oriented usage now relies on the stricter schema-validation and metadata helpers introduced in this release. Existing workflow inputs may need to be normalized to the canonical columns and metadata expected by `validate_strict_schema`, `build_conditions`, and `run_staged_pipeline`.
- Breaking: workflow exports are now organized into structured output directories such as `tables/`, `params/`, `diagnostics/`, and `figures/`. Downstream scripts that assumed the older flat export layout should be updated.

### Added

- Added staged pipeline execution with checkpoint/manual modes via `run_staged_pipeline`.
- Added workflow configuration, manifests, QC reporting, and resume support.
- Added population and cell-line stage templates with parameter inheritance helpers.
- Added simulation sweep utilities and joint fitting helpers for multi-state or multi-dataset workflows.
- Added bootstrap stage uncertainty summaries and stricter data/schema validation helpers.

### Notes

- For registration release notes, reference this changelog and call out the breaking items above.