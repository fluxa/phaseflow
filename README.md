# PhaseFlow

Exploratory study of whether interference in a locally propagated complex
field is a useful routing signal on classical networks. Start with
`docs/phaseflow_kickoff.md` and `docs/phaseflow_research_charter.md`.

## Reproduce Experiment 0 (smoke test)

```bash
uv sync
uv run pytest -q
uv run python -m phaseflow.experiment0   # the predeclared sweep, about 8 seconds
uv run python -m phaseflow.analysis      # post-hoc tables, written to posthoc.json
```

Outputs go to `results/experiment0-smoke/`. The run is deterministic. Three
large CSVs are git-ignored and are recreated by the first command.
Design: `docs/superpowers/specs/2026-10-02-experiment0-smoke-design.md`.
Findings and the open decision: `docs/experiment-0.md`.
