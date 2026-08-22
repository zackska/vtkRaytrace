# Job search agent

A specification, not an implementation. `SPEC.md` describes an agent that polls
Stockholm-area career pages, scores each opening against a fixed candidate
profile, drafts a tailored CV and cover letter, and queues the package for a
human to review and submit.

It lives in this repository because the candidate profile it scores against is
largely the same body of work the rest of the repository documents — CFD, optical
diagnostics, instrumented test hardware. Nothing here is imported by the tracers
or the analysis scripts, and nothing here affects the build. Note that a GitHub
release archives the whole tree, so this directory does travel into the Zenodo
record — which is the reason §3 keeps the contact and compensation values out of
`SPEC.md` and in a gitignored file instead.

## Files

| File | Purpose |
|---|---|
| `SPEC.md` | The specification. §2 lists the operating principles that override everything else; §9 explains why the agent stops before the submit button. |
| `profile.local.md.example` | Template for the values `SPEC.md` §3 deliberately does not commit. |

## Local, uncommitted state

Three paths are gitignored, because §12 of the spec requires the candidate's
personal data to stay local and this repository is public:

| Path | Holds |
|---|---|
| `profile.local.md` | Phone, email, contract rate, salary reference |
| `queue/` | Drafted application packages (§8) |
| `tracker.csv` | One row per role, incl. agency and consent columns (§10) |

Copy `profile.local.md.example` to `profile.local.md` and fill it in before
running anything against the spec.

## Status

Specification only — no code yet. The scoring model (§5) and the CV variant
library (§6) are the two pieces that need to exist before any of the polling in
§4 is worth automating.
