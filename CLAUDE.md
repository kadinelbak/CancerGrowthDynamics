# CLAUDE.md

Research repo for cancer cell growth modelling (A2780 naive / cisplatin-resistant
lines, mono- and co-culture, treated and untreated). Mostly Julia notebooks and
scripts plus some Python data processing.

## Layout
- `Processing Files/`: Python/Julia scripts that turn raw plate data into
  `Processed_Datasets/`. Many expect to be run from the dataset's directory.
- `Modeling_Approaches/01_mechanical_manual/`: hand-built ODE fits in notebooks.
- `Modeling_Approaches/02_mechanical_automatic_package/`: package-driven fits;
  the only part covered by CI (`.github/workflows/validate-strobl-models.yml`,
  `julia --project=. -e 'using Pkg; Pkg.test()'` from that folder).
- `Modeling_Approaches/03_automatic_neural_network/`: neural network approaches.
- `docs/`, `index.html`, `build_*_report.jl`: published reports.
- `plotting_guidelines.md`: required style for every figure (CairoMakie, PDF+SVG+PNG).
- Each modelling folder has its own `Project.toml`; activate the one next to the
  file you are running.

## Model delegation

Keep the main session's model for work that needs judgment. Hand mechanical,
well-scoped work to cheaper subagents with the Agent tool's `model` parameter.

**Use `model: "haiku"` for** lookups and mechanical edits where the answer is
easy to check:
- Finding files, functions, or where a symbol is used (use the `Explore` agent).
- Listing or summarizing directory contents, CSV headers, or output files.
- Renaming, reformatting, fixing typos, updating paths or docstrings.
- Running a command and reporting what it printed (tests, `Pkg.status()`, git log).

**Use `model: "sonnet"` for** contained implementation work with a clear spec:
- Writing or updating a single function, script, or test from an explicit description.
- Plotting scripts that follow an existing pattern or style guide.
- Summarizing a long notebook, log, or report into the facts you need.
- Drafting docs, READMEs, changelog entries, or PR descriptions.

**Keep on the main model** (do not delegate):
- Choosing or changing model structure: ODE right-hand sides, priors, bounds,
  identifiability, likelihoods, BIC/AIC comparisons.
- Debugging wrong numbers, failed fits, or solver instability.
- Interpreting results and deciding what they mean for the science.
- Cross-file design decisions and anything the user will act on directly.
- Reviewing and integrating whatever a subagent returns.

**How to delegate well**
- Give the subagent everything it needs in the prompt: file paths, the exact
  change, and what "done" looks like. It does not see this conversation.
- Ask for a short result (paths, line numbers, a diff summary), not a dump.
- Run independent subagents in parallel in one message.
- Verify a subagent's edits before building on them: read the diff and run the
  relevant check. If a cheap model gets it wrong twice, do it yourself.
- Don't delegate a task smaller than writing its prompt. A single grep or a
  one-line edit is faster to do directly.

**Repo-specific examples**
- Haiku: "Which notebooks read `Processed_Datasets/Treated Coculture`?",
  "List the columns in every CSV under `outputs/csv`", regenerating a figure
  after a path change.
- Sonnet: a new plotting script that follows `plotting_guidelines.md`, a
  processing script for a new well layout that mirrors an existing one,
  summarizing a notebook's results section.
- Main model: choosing between logistic / Gompertz / theta-logistic / death or
  interaction terms, fixing a fit whose parameters hit bounds, deciding whether
  a BIC difference is meaningful.

Large binary outputs (PNG, PDF, Manifest.toml) don't need to be read in full;
have a cheap subagent report what is there instead.
