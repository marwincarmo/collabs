# Deviations from Gross et al. and from the preregistration

Every departure the analysis code embodies, why it exists, and what it costs.
Each is switched or asserted in the code rather than left implicit; the
`` names refer to `R/00_setup.R`.

Status legend: **instrument** = forced by how the Brazilian survey was built;
**analytic** = a choice made by this team; **prereg** = differs from the
registered plan.

---

## 1. No no-inner-speech branch  *(instrument)*

Gross et al. asked the thought-content items twice — once of participants
reporting inner speech, once of those reporting none — and coalesced the pairs
into `clarity_all`, `valence_all`, `reaction_all`, `interesting_all`. The
Brazilian instrument asks them only on the inner-speech branch.

**Consequence.** The `_all` coalescing has no counterpart here and is not
implemented. More importantly it removes the base Gross used for person-mean
centering (see §2).

**Verified by.** `R/02_recode.R` prints a raw-label x value table for every
thought item, over inner-speech prompts only.

## 2. Person-mean centering base  *(analytic, follows from §1)*

Gross computed person means over **all** moments using the coalesced items.
Here they are computed over **inner-speech moments only** —
`centering_base = "inner_speech"`.

**Consequence.** Person means rest on roughly 13 moments rather than roughly 42,
so the `_cw2` predictors carry more sampling noise than in the original, and
`wb_cw2` is not comparable in construction to Gross's.

**Cost, measured.** Mood is asked at every prompt, so its person mean *could*
be taken over all moments. `centering_base = "all_moments"` does exactly
that and is reported in `output/sensitivity.csv`. It leaves Models A and B
untouched (their outcome is raw `wb`) but moves the total-effect model
materially, so the sensitivity table includes that model specifically to make
the difference visible.

## 3. No `reaction` item  *(instrument)*

Gross measured thought reaction; the Brazilian ESM survey does not. No analysis
in the prereg depends on it.

## 4. Branch integrity is not enforced by the platform  *(instrument)*

A small number of moments report **no** inner speech yet carry thought-content
answers — 3 of 303 (1.0%) in the pilot — most plausibly participants who
backtracked and changed their inner-speech answer.

**Resolution.** The inner-speech item is authoritative. Such moments are
flagged as `branch_odd` in `R/02_recode.R`, which prints every offending row,
and excluded by the `inner_speech == "Yes"` filter in `R/04_derive.R`.

## 5. FFMQ has 14 items, not 15  *(instrument, differs from prereg)*

The prereg names the FFMQ-15. The instrument administers 14: three each for
Observe, Describe, Acting with Awareness and Non-judging, but only **two** for
Non-reactivity.

**Consequence.** Facet scores are means, so they stay on the 1–5 metric and
remain interpretable. The **total is not comparable to a published FFMQ-15
total** and should be described as a 14-item mean.

**Open.** Confirm with the team whether an item was dropped deliberately. If it
was dropped in error and the final instrument restores it, the `ffmq_facet` and `ffmq_rev`
vectors in `R/02_recode.R` need a fifteenth entry, and the item table printed
just below them will show the change.

## 6. Attention-state contrast sign  *(analytic)*

`contr.sum(2)` assigns +1 to the **first** factor level. These scripts set
`attention_levels = c("Present", "MW")`, so `attention1` is positive when
mood is higher while present — the direction H1 predicts.

`R/power/power_gross2024.R` let the levels fall out alphabetically
(`"MW"`, `"Present"`), which gives `attention1` the opposite sign. Power depends
only on the magnitude, so the power analysis is unaffected, but **the sign of
`attention1` is not comparable between `R/power/` and `R/0*.R`**.

`attention1` is half the Present − MW difference, and every script that prints
it doubles it explicitly. `R/09_check_recovery.R` recovers a known effect and
prints the sign check. **Superseded as the primary coding by §13**, which drops
the contr.sum factor in favour of a 0/1 within/between split; `attention1`
survives only in the `Blended` row of the sensitivity table.

## 7. Estimator: REML  *(prereg)*

The prereg specifies Restricted Maximum Likelihood, and `reml = TRUE`
follows it. `R/power/power_gross2024.R` used `REML = FALSE` for simulation
speed. `R/07_sensitivity.R` refits under ML so the difference is visible.

## 8. Exclusion criteria are operationalised, not transcribed  *(prereg)*

The prereg names three criteria: not taking the study seriously, failing
attention checks, and having no inner-speech prompt. Gross had separate validity
checks (`VAL1`/`VAL2`) and three attention checks (excluding at more than one
failure); the Brazilian instrument has **one** seriousness/effort item per
survey, so the first two criteria collapse into a single effort screen.

- **Threshold.** `effort_pass_max = 2`: options 1–2 pass, 3–4 fail. Option 2
  ("went through quickly but still answered correctly") explicitly denies the
  thing the prereg excludes on; options 3 and 4 assert it. The strict reading
  (option 1 only) is reported in the sensitivity table.
- **Level.** Participant-level only. The per-prompt effort item is asked *only*
  on the inner-speech branch, so excluding moments on it would thin the analysis
  sample non-randomly with respect to the outcome. The moment-level screen is a
  sensitivity analysis (`drop_prompts_on_effort`).
- **Age.** Not a prereg criterion, but consent restricts the study to adults and
  Gross applied the same screen. Applied as its own logged step
  (`min_age`), and reported without it in the sensitivity table.

## 9. Prompt order is derived, not exported  *(instrument)*

`Occasion within Day` is constant (`1`) for every row of the pilot export, so
within-day prompt order is reconstructed from `Start Date`. Which source was
used is printed by `R/02_recode.R` and carried on the data as the
`occasion_source` attribute rather than assumed. If the final
export populates the column, the pipeline uses it automatically.

## 10. Gender is collapsed for moderation  *(analytic)*

The item has seven options. Non-cis and undisclosed responses are very small
cells (4 of 204 in the pilot), so the preregistered demographic moderation uses
three levels — Woman, Man, Other / undisclosed (the `gender3` recode in
`R/02_recode.R`). The
ungrouped counts are reported alongside. The collapsed category is
heterogeneous and its interaction estimate should not be interpreted as an
effect for any specific group.

## 11. Demographics: no children item  *(instrument)*

The consent form describes collecting age, presence of children and sex. The
day-1 survey collects age, gender identity and the MacArthur subjective social
status ladder. No analysis in the prereg depends on a children variable.

## 12. Personal data in the raw exports  *(handling)*

The exports carry `ip address`, `Location - Lat` and `Location - Long`.
`R/02_recode.R` keeps only the columns it needs and drops all three, so nothing downstream — including everything written to `output/` —
contains them. The raw files still do: `.gitignore` excludes `data/raw/`, and
those columns must be stripped before any public data deposit.

## 13. Attention state is split within/between, not entered raw  *(analytic, differs from Gross et al. and from the prereg)*

**What Gross did.** Two different codings, in two different places.

- Table 3 (`model1a_fit`/`model1b_fit`) enters `attention` raw alongside
  `(1 | PID)`, with raw mood as the outcome.
- The mediation chain (`fit.totaleffect1`, `fit.mediator1`, `fit.dv1`) uses
  **person-mean-centred** mood and valence as outcomes.

**Why neither is used here.**

A random intercept does not absorb the between-person component of a level-1
predictor. With `attention` entered raw, the coefficient is a
precision-weighted blend of the within-person effect (the same person, Present
versus Mind-Wandering) and the between-person association (people who are
Present more often report higher mood overall). Only the first is what H1 and H2
claim. See Neuhaus & Kalbfleisch (1998), Enders & Tofighi (2007), Curran & Bauer
(2011).

Centring the *outcome* does isolate the within-person effect, but at two costs.
It forces every person's mean to exactly zero, leaving `(1 | PID)` with no
variance to estimate — every such model is singular by construction, not because
of anything about the data. And a participant with a single prompt then
contributes an outcome of exactly zero against a predictor that still varies,
which drags the slope toward zero. That second cost is severe here and was
negligible for Gross: their participants averaged ~13 inner-speech prompts,
where half of ours currently have one.

**What is done instead.** `R/04_derive.R` builds `att_num` (1 = Present,
0 = MW), its person mean `att_mean2`, and the split `att_cw2` (within) /
`att_cb2` (between), on the same `_cw2`/`_c2`/`_cb2` pattern as every other
variable. Both terms enter every H1 and H2 model. `att_cw2` is the **full**
Present − MW within-person difference — `att_num` is 0/1, so unlike
`attention1` (§6) there is nothing to double.

A single-prompt participant has `att_cw2` exactly 0: they contribute nothing to
the within slope while still informing the intercept and the variance
components, which is the correct treatment rather than a harmful one.

**Consequence, measured** (pilot: 183 prompts, 95 participants, Model A
covariates):

| coding | Present − MW | p | singular |
|---|---|---|---|
| `att_cw2` + `att_cb2` (used here) | +0.526 (SE 0.348) | .134 | no |
| `attention` raw, contr.sum (Gross Table 3) | +0.894 (SE 0.237) | .0002 | no |
| `wb_cw2` outcome (Gross mediation chain) | +0.222 (SE 0.154) | .150 | **yes** |

The between-person coefficient is +1.215 (p = .0003), which is what inflates
the raw-coding row. The preregistered direction is significant only under the
coding that confounds it with a between-person difference.

**Confirmed by simulation.** `R/09_check_recovery.R` generates data with a known
within-person effect of +0.500 and a deliberate between-person confound, then
fits all three codings over 25 replications. Averaged bias:

| coding | 300 x 13 prompts | half the sample at one prompt |
|---|---|---|
| `att_cw2` + `att_cb2` | +0.008 | +0.017 |
| `attention` raw | +0.075 | +0.109 |
| `wb_cw2` outcome | −0.077 (singular 100%) | −0.100 (singular 100%) |

The two rejected codings miss in opposite directions, and both get worse as the
data thin — which is the regime this study is currently in.

**Cost.** This is a departure from the registered analysis plan, which follows
Gross et al. All three codings are reported for every data-construction variant
in `R/07_sensitivity.R` and `output/sensitivity.csv`, so the comparison above is
reproduced on the final data rather than asserted from the pilot.

**One compromise inside H2.** `mediation::mediate()` requires the mediator to be
the same variable in the a-path and b-path models, and the a-path's outcome
cannot be `valence_cw2` without reintroducing the centred-outcome problem. The
mediator is therefore raw `valence`, making the b path a within/between blend.
`R/05_models.R` prints it both ways; in the pilot they differ by 0.004
(+0.518 blended, +0.522 within). Re-check on the final data.

**Not affected.** `clarity_cw2`, `interesting_cw2` and `valence_cw2` are already
person-mean centred, and a person-mean-centred predictor is orthogonal to its
own person mean, so omitting their between-person counterparts does not bias
their within-person coefficients. Only `attention` needed the split.
