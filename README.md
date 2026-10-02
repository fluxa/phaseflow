# PhaseFlow

Exploratory study of whether interference in a locally propagated complex
field is a useful routing signal on classical networks. Start with
`docs/phaseflow_kickoff.md` and `docs/phaseflow_research_charter.md`.

## Reproduce Experiment 0 (smoke test)

```bash
uv sync
uv run pytest -q
uv run python -m phaseflow.experiment0
```

Outputs go to `results/experiment0-smoke/`. The run is deterministic.
Design: `docs/superpowers/specs/2026-10-02-experiment0-smoke-design.md`.
Findings and the open decision: `docs/experiment-0.md`.
