# Deviations from Gross et al. and from the preregistration

Every departure the analysis code embodies, why it exists, and what it costs.
Each is switched or asserted in the code rather than left implicit; the
`CFG$` names refer to `R/00_config.R`.

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

**Verified by.** `R/tests/test_ingest.R` asserts the thought items are absent
for more than 95% of non-inner-speech moments.

## 2. Person-mean centering base  *(analytic, follows from §1)*

Gross computed person means over **all** moments using the coalesced items.
Here they are computed over **inner-speech moments only** —
`CFG$centering_base = "inner_speech"`.

**Consequence.** Person means rest on roughly 13 moments rather than roughly 42,
so the `_cw2` predictors carry more sampling noise than in the original, and
`wb_cw2` is not comparable in construction to Gross's.

**Cost, measured.** Mood is asked at every prompt, so its person mean *could*
be taken over all moments. `CFG$centering_base = "all_moments"` does exactly
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
flagged (`branch_inconsistent`), counted in the run log, and excluded by the
`inner_speech == "Yes"` filter. `R/tests/test_ingest.R` asserts none reaches the
analysis sample and that the rate stays under 5%.

## 5. FFMQ has 14 items, not 15  *(instrument, differs from prereg)*

The prereg names the FFMQ-15. The instrument administers 14: three each for
Observe, Describe, Acting with Awareness and Non-judging, but only **two** for
Non-reactivity.

**Consequence.** Facet scores are means, so they stay on the 1–5 metric and
remain interpretable. The **total is not comparable to a published FFMQ-15
total** and should be described as a 14-item mean.

**Open.** Confirm with the team whether an item was dropped deliberately. If it
was dropped in error and the final instrument restores it, `FFMQ_KEY` in
`R/01_codebook.R` needs the fifteenth row and `01_codebook.R`'s `n_choices`
assertions will catch the change on ingest.

## 6. Attention-state contrast sign  *(analytic)*

`contr.sum(2)` assigns +1 to the **first** factor level. These scripts set
`CFG$attention_levels = c("Present", "MW")`, so `attention1` is positive when
mood is higher while present — the direction H1 predicts.

`R/power/power_gross2024.R` let the levels fall out alphabetically
(`"MW"`, `"Present"`), which gives `attention1` the opposite sign. Power depends
only on the magnitude, so the power analysis is unaffected, but **the sign of
`attention1` is not comparable between `R/power/` and `R/0*.R`**.

`attention1` is half the Present − MW difference. The reporting layer emits
`contrast_present_minus_mw` (twice the coefficient, with a matching interval)
so the two are never confused. `R/tests/test_recovery.R` asserts both the sign
and the doubling.

## 7. Estimator: REML  *(prereg)*

The prereg specifies Restricted Maximum Likelihood, and `CFG$reml = TRUE`
follows it. `R/power/power_gross2024.R` used `REML = FALSE` for simulation
speed. `08_sensitivity.R` refits under ML so the difference is visible.

## 8. Exclusion criteria are operationalised, not transcribed  *(prereg)*

The prereg names three criteria: not taking the study seriously, failing
attention checks, and having no inner-speech prompt. Gross had separate validity
checks (`VAL1`/`VAL2`) and three attention checks (excluding at more than one
failure); the Brazilian instrument has **one** seriousness/effort item per
survey, so the first two criteria collapse into a single effort screen.

- **Threshold.** `CFG$effort_pass_max = 2`: options 1–2 pass, 3–4 fail. Option 2
  ("went through quickly but still answered correctly") explicitly denies the
  thing the prereg excludes on; options 3 and 4 assert it. The strict reading
  (option 1 only) is reported in the sensitivity table.
- **Level.** Participant-level only. The per-prompt effort item is asked *only*
  on the inner-speech branch, so excluding moments on it would thin the analysis
  sample non-randomly with respect to the outcome. The moment-level screen is a
  sensitivity analysis (`CFG$drop_prompts_on_effort`).
- **Age.** Not a prereg criterion, but consent restricts the study to adults and
  Gross applied the same screen. Applied as its own logged step
  (`CFG$min_age`), and reported without it in the sensitivity table.

## 9. Prompt order is derived, not exported  *(instrument)*

`Occasion within Day` is constant (`1`) for every row of the pilot export, so
within-day prompt order is reconstructed from `Start Date`. Which source was
used is recorded in `output/run_manifest.txt` rather than assumed. If the final
export populates the column, the pipeline uses it automatically.

## 10. Gender is collapsed for moderation  *(analytic)*

The item has seven options. Non-cis and undisclosed responses are very small
cells (4 of 204 in the pilot), so the preregistered demographic moderation uses
three levels — Woman, Man, Other / undisclosed (`GENDER_COLLAPSE`). The
ungrouped counts are reported alongside. The collapsed category is
heterogeneous and its interaction estimate should not be interpreted as an
effect for any specific group.

## 11. Demographics: no children item  *(instrument)*

The consent form describes collecting age, presence of children and sex. The
day-1 survey collects age, gender identity and the MacArthur subjective social
status ladder. No analysis in the prereg depends on a children variable.

## 12. Personal data in the raw exports  *(handling)*

The exports carry `ip address`, `Location - Lat` and `Location - Long`.
`03_clean.R` keeps an explicit whitelist of metadata columns and drops all
three, so nothing downstream — including everything written to `output/` —
contains them. The raw files still do: `.gitignore` excludes `data/raw/`, and
those columns must be stripped before any public data deposit.
