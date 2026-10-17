# genieLOT algorithm specification

Version 0.1.0 · Applies to AACR Project GENIE BPC data (tested on NSCLC 2.0-public)

This document states every rule genieLOT uses, so results can be reviewed and
reproduced. Every setting below can be changed by the user; defaults are shown.

## 1. Inputs

| Table (BPC file) | Used for | Key columns (default names) |
|---|---|---|
| Cancer diagnosis (`cancer_level_dataset_index`) | stage, distant metastasis, follow-up | `record_id`, `ca_seq`, `stage_dx_iv`, `dx_to_dmets_*`, `tt_os_dx_*` |
| Regimens (`regimen_cancer_level_dataset`) | drugs and dates | `regimen_drugs`, `dx_reg_start_int`, `dx_reg_end_all_int` |
| Imaging (`imaging_level_dataset`) | progression | `dx_scan_days`, `image_overall` |
| Oncologist notes (`med_onc_note_level_dataset`) | progression | `dx_md_visit_days`, `md_ca_status` |
| Optional: pathology reports (`pathology_report_level_dataset`) | surgery dates (proxy) | `path_proc_type` = "Surgical pathology", `dx_path_proc_days` via `bpc_surgery_events()` |
| Optional: radiation or other events | initial therapy | any table via `therapy_events()` |

All times are days since diagnosis of the index cancer. A case is one
`record_id` + `ca_seq`. Column names are set in one place, `bpc_columns()`.

## 2. Data preparation

1. Each regimen is split into one record per drug, sharing the regimen's
   start and end dates.
2. **Missing end dates** are filled, in order, from: the last recorded dose of
   any drug in the regimen (`dx_drug_end_or_lastadm_int_*`); the first drug's
   end (`dx_reg_end_any_int`); else the start date. The source is recorded.
3. A **progression event** is any imaging or oncologist-note record whose
   status contains "progress" (case-insensitive).
4. **Early stage** = stage I–III; stage IV and unknown are kept separate.

## 3. Index date (where counting starts)

Chosen with `set_index(event = ...)`:

| Event | Index date |
|---|---|
| `recurrence` (default) | Stage I–III: earliest of distant metastasis or progression on or after day 90. Stage IV: diagnosis. |
| `progression_after_therapy` | First progression after initial therapy (drug, surgery or radiation) has started, on or after day 90. This follows the challenge definition. Stage IV: diagnosis. |
| `diagnosis` | Day 0 (the naive approach). |
| `metastasis` | First distant metastasis. |
| `progression` | First progression event. |
| `initial_therapy` | First treatment of any kind. |
| `custom` | User-supplied dates. |

Options: progression source (imaging, notes, either, or both agreeing within
30 days); the 90-day minimum (`min_days`); stage IV handling.

**Upstaged patients:** stage I–III patients whose distant metastasis is found
within `min_days` of diagnosis had spread at diagnosis; they are indexed at
diagnosis, like stage IV.

Patients with no index date (for example early-stage patients who never
recurred) receive no lines.

## 4. Line of therapy rules

Drugs starting before the index date (minus `lookback`, default 0) are
excluded. The remaining drugs, in date order, are assigned to lines:

1. **Regimen window (28 days):** drugs that start within 28 days of a line's
   first drug belong to that line.
2. **New drug:** after the window, a drug not already in the line starts a
   new line, unless it is an allowed add-on or a swap.
3. **Treatment gap (120 days):** a break longer than 120 days after all of the
   line's drugs have stopped starts a new line, even if the same drugs restart.
4. **Dropping a drug** (e.g. maintenance after induction) never starts a new line.
5. **Swaps:** listed drug groups (default carboplatin/cisplatin) are treated as
   the same drug.
6. **Add-ons:** listed drugs may be added later without a new line.
7. Optional: start a new line after every progression event (off by default;
   published comparisons found progression data did not improve accuracy).
8. Optional: drop listed drugs entirely (e.g. "Investigational Drug").

Drug names match by prefix and ignore case ("Goserelin" matches "Goserelin
Acetate").

**Disease presets** (`lot_preset()`): NSCLC, CRC (fluorouracil/capecitabine
swap; bevacizumab and leucovorin add-ons), BrCa (aromatase-inhibitor swaps;
ovarian suppression and HER2 antibody add-ons; 180-day gap). Presets are
starting points and should be reviewed by a clinician for each analysis.

## 5. Outputs

- `calculate_lot()`: one row per patient line: regimen, start, end, days from
  index, duration, time to next line.
- `summarize_lot()`, `compare_lot()` (naive vs chosen definition),
  `validate_index()` (vs BPC's `dmets_stage_i_iii` flag), `lot_attrition()`,
  `lot_sensitivity()`, `cohort_summary()`.

## 6. Privacy

Summary outputs mask counts below 10 (`suppress_small()`). The app's privacy
mode hides patient timelines. No patient-level data is written by the
reproduce script. GENIE data must not be committed to the repository.

## 7. Known limitations

- Surgery dates come from surgical pathology reports, a proxy: some are diagnostic
  procedures rather than treatment. Radiation is used only when supplied.
- Masked trial drugs ("Investigational Drug") cannot be told apart.
- Progression is taken from structured status fields; free-text notes are not read.
- Presets for CRC and BrCa have not yet been reviewed by clinicians.

## 8. Validation

Unit tests cover each rule (`tests/testthat`). On NSCLC 2.0-public, the
recurrence index handled all 583 patients flagged by BPC's
`dmets_stage_i_iii` (536 as recurrences, 47 upstaged), with 88.0% overall
agreement on early-stage patients. See RESULTS.md.

## References

Default window and new-drug rules follow the published line-of-therapy
algorithm in: An Automated Line-of-Therapy Algorithm for Adults With
Metastatic Non–Small Cell Lung Cancer: Validation Study Using Blinded Manual
Chart Review. JMIR Medical Informatics, 2021;9(10):e29017.
