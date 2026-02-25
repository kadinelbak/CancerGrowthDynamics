# Julia Plotting Guidelines (Publication Standard)

These are the default plotting rules for this project unless explicitly overridden.

## Core defaults
- Use **CairoMakie** for all plots unless told otherwise.
- Provide **complete, runnable Julia scripts** (no partial snippets).
- Use a clean scientific style appropriate for journals.

## Style rules (mandatory)
- Export vector-safe outputs: PDF + SVG.
- Also export high-resolution PNG.
- Typical figure size: `size = (600, 420)` (about 500–700 px width equivalent).
- No chartjunk: keep grid off unless explicitly requested.
- Keep axis spines visible (`xspinevisible = true`, `yspinevisible = true`) and axis line width around 1–1.5.
- Keep ticks clear with consistent lengths.
- Font sizes:
  - axis labels: 12–14
  - ticks: 10–12
  - panel titles: 13–16
- Use neutral colors by default; use distinct mapping when categories require it.
- Never rely on color alone: also use markers and/or line styles.
- Show individual datapoints whenever possible.
- Error bars must explicitly state **SEM** or **CI**.
- Avoid bar plots when raw data is available.

## Preferred data presentation
- Scatter + mean ± SEM/CI.
- Box/violin + points.
- Regression line + confidence band.
- For grouped data:
  - maintain consistent ordering across figures.
  - maintain consistent color mapping across thesis/project.
- Apply slight jitter to reduce overlap.

## Layout rules
- Use `Figure` + `Axis` explicitly.
- Multi-panel figures must have aligned axes.
- Use shared axis limits when comparing conditions.

## Export block (always include)
```julia
save("figure_name.pdf", fig)
save("figure_name.svg", fig)
save("figure_name.png", fig; px_per_unit=3)
```

## Reproducibility (always include)
```julia
using Random
Random.seed!(1)
```

## Code quality
- Define and use a reusable publication theme.
- Use descriptive variable names.
- Add brief comments for each major section.
- Avoid interactive-only features.

## Output format expectation
- Return one complete script.
- End with a short **"How to adapt this"** section.
