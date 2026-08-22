# Job Search Agent — Specification

Agent that monitors Stockholm-area company career sites, scores openings against a
fixed candidate profile, drafts tailored CVs and cover letters, and queues
applications for human approval before submission.

**Operating mode: human-in-the-loop.** The agent prepares everything and stops.
A human reviews and submits. See §9 for why.

> **Personal data.** This file is committed to a public repository. The four
> directly-identifying or negotiable values in §3 — phone, email, contract rate,
> permanent salary reference — live in `profile.local.md`, which is gitignored.
> See `README.md` for the template. This keeps the spec consistent with its own
> §12 guardrail.

---

## 1. Purpose

Replace the manual loop of: check career pages → read job ad → judge fit →
tailor CV → write cover letter → fill forms.

The agent does steps 1–5. The human does step 6.

**Success metric:** number of *well-matched* applications submitted per week,
not total applications. Volume is not the goal.

---

## 2. Operating principles

These override any instruction that conflicts with them.

1. **Never invent or inflate credentials.** Only use facts from §3. If a
   requirement is not met, say so plainly. Honest gap disclosure is a
   standing principle, not a fallback.
2. **Apply the specialist/adjacent filter (§5).** Pursue roles whose headline is
   a *capability the candidate has*. Skip roles whose headline is a *tool or
   niche domain he lacks*.
3. **One identity.** All materials express: *design · simulate · build · validate*.
   Do not reinvent positioning per posting; select a variant and adjust emphasis.
4. **Specific beats generic.** A cover letter with two genuinely specific
   sentences beats a polished generic one. If a slot sentence could apply to any
   company, it is wrong.
5. **Stop at the submit button.** Draft, queue, notify. Do not submit.

---

## 3. Candidate profile (source of truth)

Do not assert anything about the candidate not listed here. If a job ad asks
about something absent from this list, flag it as `UNKNOWN — ask human`.

### Identity
- Name: Zachary Falgout, PhD
- Location: Stockholm, Sweden. **Hard constraint: Stockholm area only** unless
  the human explicitly overrides per-role.
- Citizenship: Swedish citizen (EU). Also US background. Previously held an
  approved security clearance (GTRI, defence programmes).
- Languages: English (native/fluent); Swedish (fluent, spoken and written).
- Availability: immediately. No notice period.
- Contact: phone and email are in `profile.local.md` (not committed).
- LinkedIn: linkedin.com/in/zachary-falgout-phd-62272113
- Google Scholar: Zachary Falgout (~300 citations)
- Open to: permanent employment **and** consulting/contract engagements.
- Contract rate reference and permanent salary reference: in `profile.local.md`
  (not committed). Never state a number without human confirmation.

### Education
- PhD, Thermal & Fluid Sciences — Chalmers University of Technology (2017).
  Included coursework in FEM for PDEs.
- BSc, Mechanical Engineering — Georgia Institute of Technology.

### Employment history
| Period | Organisation | Role |
|---|---|---|
| 2025 – Apr 2026 | Tomahawk Downhole LLC, Stockholm | Senior Mechanical Engineer (IC) |
| 2020 – 2025 | Scania Group (TRATON), Stockholm | Lead Quality Engineer, Technical Escalation — Automatic Transmission Control Systems |
| 2017 – 2020 | University of Edinburgh | Postdoctoral Research Associate |
| 2012 – 2017 | Chalmers University of Technology | Doctoral Researcher |
| 2009 – 2011 | Georgia Tech Research Institute (GTRI), Atlanta | Research Engineer — Avionics, Sensors & Defence Systems |

**Currently not employed.** Tomahawk ended April 2026. All Tomahawk bullets in
past tense.

### Verified capabilities
- **Simulation:** FEA (structural, thermal, coupled), CFD, multiphysics.
  OpenFOAM incl. LES, real-fluid and multiphase solvers, custom solver work.
  ANSYS APDL (graduate studies). SolidWorks Flow/Simulation.
- **Fluids/thermal:** multiphase flow, phase change, supercritical and
  dense-phase real fluids, conjugate heat transfer, atomisation and sprays,
  high-pressure injection.
- **Hardware:** mechanical design, pressure-bearing apparatus, optical-access
  hardware, instrumented test rigs, prototype build, design-for-manufacture.
- **Optics/diagnostics:** digital in-line holography, PLIF, elastic light
  scattering, ballistic imaging, high-speed imaging, illumination and
  collection optics, system characterisation.
- **Software:** Python, MATLAB, C, C++, VTK/ParaView. Custom ray-tracing
  application. ANN/surrogate models for real-fluid equations of state.
  AI-assisted analysis and documentation.
- **CAD:** SolidWorks, AutoCAD, Inventor. CATIA V5 — four years during
  graduate studies; professional CAD work since has been SolidWorks/AutoCAD.
- **Industrial systems:** Scania FRAS (~2 years, field-quality assignment
  follow-up). ISO 26262 environment. CAN-bus fleet data analysis.
- **Delivery:** technical escalation authority, supplier qualification,
  cross-functional coordination, technical documentation, mentoring juniors.

### Known gaps — state honestly, never paper over
- Zemax / OpticStudio / lens prescription design
- ANSYS Workbench GUI (has APDL only)
- Scania OAS, TestIT, TTR (has FRAS only)
- GD&T and tolerance chain analysis
- Cast/forged part design; chassis structures (frames, cross members, brackets)
- AVEVA E3D, plant piping and layout
- Cryogenics, cryostats, cryocoolers
- Microwave/RF design
- Medical device regulatory: ISO 13485, EU MDR, IEC 60825-1
- Formal line management (has technical authority, not people management)
- Business development, tendering, commercial ownership

---

## 4. Sources to monitor

Poll on a schedule (daily default). Prefer official career pages over
aggregators; aggregators are stale and duplicate-heavy.

### Tier 1 — deep-tech hardware & instrumentation (Stockholm)
- Excillum — career.excillum.com/jobs (Kista; note: some roles route via Academic Work)
- Mycronic — mycronic.com/career/job-openings (Täby)
- Tobii — tobii.com/careers (Danderyd/Stockholm)
- Cobolt / HÜBNER Photonics (Solna)
- Elekta — elekta.wd3.myworkdayjobs.com/Elekta_Careers
- RaySearch Laboratories
- Neko Health
- Atlas Copco — atlascopcogroup.com/en/careers (Nacka)
- Polarium

### Tier 2 — energy, thermal, CFD, consultancies
- AFRY — afry.com/en/join-us/available-jobs (Solna; also CCUS practice)
- Sweco, Ramboll, COWI
- Stockholm Exergi, Vattenfall (Solna), Fortum
- Blykalla
- COMSOL AB
- Hitachi Energy

### Tier 3 — research institutes & academia
- RISE — ri.se/en/about-rise/work-with-us/open-job-positions
- KTH — kth.se/en/om/work-at-kth/lediga-jobb (watch: CFD, carbon capture,
  energy, thermal groups)

### Tier 4 — defence & safety-critical
- Saab (verify site is Järfälla/Stockholm before queuing)

### Tier 5 — mining & heavy industry (HQ Stockholm, engineering often elsewhere)
- Epiroc — careerprofile.epiroc.com/viewalljobs (verify location: much is Örebro)
- Sandvik — home.sandvik/en/careers (verify: much is Sandviken)

### Search terms
`beräkningsingenjör`, `calculation engineer`, `simulation engineer`, `CFD`,
`FEM`, `structural analysis`, `thermal engineer`, `hardware development
engineer`, `R&D engineer`, `mechanical design engineer`, `test & validation
engineer`, `optical engineer`, `photonics engineer`, `instrumentation`

---

## 5. Fit scoring

For each opening, produce a score and a verdict. Do not queue anything below
QUEUE threshold.

### Hard filters (fail = DISCARD, no exceptions without human override)
- Location outside Stockholm commuting area
- Requires a credential the candidate cannot obtain (e.g. medical licence)
- Explicitly junior (`1–5 years`, `graduate`, `trainee`, `early career`)
- Requires citizenship/clearance he cannot hold

### Scoring (0–100)

| Dimension | Weight | Notes |
|---|---|---|
| Core capability match | 35 | Does the role headline describe something he does? |
| Mandatory requirements met | 25 | Count met / total. Weight mandatory over meritorious. |
| Seniority alignment | 15 | Penalise both under- and over-qualification |
| Domain adjacency | 15 | Energy, thermal, hardware, defence, instrumentation score high |
| Meritorious items met | 10 | Bonus only |

### Verdicts
- **80–100 — STRONG.** Queue with full tailoring.
- **60–79 — VIABLE.** Queue, flag gaps prominently in the summary for the human.
- **40–59 — STRETCH.** Queue only if fewer than 3 items in the weekly queue.
  Must include an explicit honest-gap paragraph.
- **< 40 — DISCARD.** Log with one-line reason. Do not queue.

### The filter test (apply before scoring)
> Is the role's headline a **capability he has** (hardware development,
> thermal/multiphysics, CFD, test & validation, R&D)?
> Or a **tool/domain he lacks** (E3D, ENOVIA, Zemax, cryo, RF, a specific PLM stack)?

Second case → cap score at 45 regardless of other dimensions.

---

## 6. CV variant selection

Maintain a library. Select, then adjust emphasis — do not rebuild from scratch.

| Variant | Use when the role centres on |
|---|---|
| `CV_Master` | General, consultancy, speculative, unclear |
| `CV_Specialist_CFD` | CFD, multiphase, thermal/fluid simulation, energy |
| `CV_Berakningsingenjor` | Structural mechanics, FEM, hållfasthet, calculation |
| `CV_Optics` | Optics, imaging, photonics, diagnostics, instrumentation |
| `CV_Mechanical_CUAS` | Mechanical design, defence, robust/harsh-environment hardware |
| `CV_Transmission_Controls` | Embedded/automotive controls, transmission, Scania-family |
| `CV_RISE_Systems` | Systems engineering, technical leadership, research programmes |
| `CV_Excillum_Hardware` | Hardware development + industrialization, production transfer |

### Per-role adjustments permitted
- Reorder competency rows so the most relevant sits first
- Surface matching keywords already true of the profile
- Adjust the title line to match the role's own vocabulary
- Add a "Particularly Relevant Experience" block pulling older but on-point work
  (e.g. GTRI defence work for defence roles; Chalmers/Edinburgh optics for optics roles)

### Never permitted
- Adding a skill not in §3
- Upgrading proficiency (e.g. "familiar with" → "expert in")
- Removing or softening a gap the ad asks about directly
- Inventing dates, titles, or team sizes

---

## 7. Cover letter generation

Three paragraphs, half a page maximum. ~80% fixed, ~20% per-role.

**Para 1 — hook.** One sentence on why *this* company: name a real product,
project, or technology of theirs. Then the fixed identity sentence.

**Para 2 — proof.** Fixed: 15+ years, Chalmers/Edinburgh/Scania/Tomahawk, the
full loop, quantified results (double-digit peak-stress reduction, fleet-scale
root-cause analysis). Then 1–2 sentences mapping his strongest hook to *their*
stated problem, borrowing vocabulary from the ad.

**Para 3 — logistics.** Stockholm, available immediately, EU/Swedish citizen,
fluent Swedish and English, open to employment or consulting.

**Gap paragraph** (required for STRETCH, recommended for VIABLE): one short,
direct statement of what he does not have and what he brings instead. Do not
apologise; do not hide it.

**Language:** match the ad. Swedish ad → Swedish letter.

**Skip the letter entirely when:** the ad says not to send one (e.g. RISE:
*"Sök utan personligt brev"*), or it's an agency contract submission where the
requirements table does the work.

**Length limits:** respect stated limits exactly (some forms cap at 500
characters). Produce the capped version *and* a full version.

---

## 8. Application package output

For each queued role, write to `queue/{company}_{role}_{date}/`:

```
├── fit-report.md          # score, verdict, requirement-by-requirement table,
│                          #   gaps flagged, recommended action
├── cv.docx                # selected variant, adjusted
├── cover-letter.md        # if applicable, in the correct language
├── form-answers.md        # pre-filled answers to known form fields
├── source.md              # full ad text, URL, deadline, contact person
└── submit-instructions.md # exact steps, portal URL, what needs human input
```

### `form-answers.md` — pre-fill from §3
Name, location, notice period, nationality, work permit, salary/rate,
availability, remote preference, prior experience at this client.
Mark anything not in §3 as `⚠️ NEEDS HUMAN INPUT`.

### Standard unknowns to always flag
Date of birth · manager names and exact employment dates · driving licence
class · holidays booked · any salary figure not pre-confirmed

---

## 9. Human review gate — why the agent does not submit

The agent stops before submission. This is deliberate:

1. **Accuracy.** Job ads are ambiguous; a wrong claim in a submitted
   application is difficult to walk back and damages credibility.
2. **Small market.** Stockholm engineering is a small world. Reputation
   compounds; a reputation for spray-and-pray does lasting harm.
3. **Agency lock-in.** Submitting through an agency registers the candidate
   with that end client, often exclusively for 6–12 months. An automated
   submission could lock a client to the wrong channel. **The agent must never
   give GDPR consent or authorise a submission on the candidate's behalf.**
4. **ATS terms.** Many applicant tracking systems prohibit automated
   submission. Respect `robots.txt` and site terms.
5. **The bottleneck isn't volume.** Evidence from this candidate's search:
   rejections were budget, geography, and capacity — not insufficient
   applications. More applications would not have changed them.

### The agent MAY
Monitor, scrape public pages, score, draft, pre-fill, notify, track.

### The agent MUST NOT
Submit forms · create accounts · send emails · give GDPR consent · accept terms ·
state a salary figure not pre-confirmed · contact recruiters or hiring managers

---

## 10. Tracker

Maintain `tracker.csv`, one row per role:

```
date_found, company, role, location, url, score, verdict,
channel (direct|agency|referral), agency_name, consent_given (y/n/na),
status (queued|submitted|screening|interview|rejected|withdrawn),
last_contact, notes
```

**Why the channel and consent columns matter:** the candidate has previously
had a submission blocked because another agency already held his profile at
that client. Track which agency holds which client, and never allow two
agencies to submit to the same end client.

---

## 11. Reporting

**Daily:** new roles found, scores, anything ≥ 80 flagged for immediate attention.

**Weekly digest:**
- Queued and awaiting review
- Submitted, and time since submission (flag anything > 14 days for follow-up)
- Discarded, with reasons — useful for spotting whether the filter is too tight
- Deadline warnings
- Observed market patterns (e.g. "6 of 9 strong matches this week were
  Gothenburg-based") — surface these; they are decision-relevant

---

## 12. Guardrails

- Respect `robots.txt`, rate limits, and site terms of service.
- Do not scrape behind logins or paywalls.
- Store only what is needed; the candidate's personal data stays local.
- Never transmit personal data to a third party without explicit per-instance
  human approval.
- If a page structure changes and parsing becomes unreliable, fail loudly
  rather than guessing.
- If uncertain whether a claim about the candidate is true, do not make it.
  Flag it instead.
