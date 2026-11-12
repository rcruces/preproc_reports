#!/usr/bin/env python3
"""
2026_bids-eval_runs.py — minimal BIDS acquisition-family run/part checker.

The 'run' entity may ONLY be used to distinguish repeated acquisitions.
It must never be used to distinguish part-mag from part-phase.

Usage:
    python 2026_bids-eval_runs.py <bids_dataset_root> [--outdir OUTDIR]

Output:
    summary.tsv  - one row per (sub, ses, datatype, acquisition_family)
    details.tsv  - one row per finding
"""

import os
import re
import sys
import argparse
import csv
from collections import defaultdict, Counter

SUPPORTED_TYPES = {"anat", "func", "fmap", "dwi"}
VALID, WARNING, ERROR = "VALID", "WARNING", "ERROR"
SEVERITY_ORDER = {VALID: 0, WARNING: 1, ERROR: 2}

# family identity = every entity except these
MEMBER_ENTITIES = {"run", "part"}
ENTITY_RE = re.compile(r"^([A-Za-z0-9]+)-(.+)$")


###############################################################################
# SCAN + PARSE
###############################################################################

def scan_dataset(root):
    kept = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if not d.startswith(".") and d != "derivatives"]
        if os.path.basename(dirpath) not in SUPPORTED_TYPES:
            continue
        for fname in filenames:
            if fname.startswith(".") or not (fname.endswith(".nii.gz") or fname.endswith(".json")):
                continue
            kept.append(os.path.join(dirpath, fname))
    return kept


def parse_filename(path):
    filename = os.path.basename(path)
    if filename.endswith(".nii.gz"):
        stem = filename[:-len(".nii.gz")]
    elif filename.endswith(".json"):
        stem = filename[:-len(".json")]
    else:
        return None

    datatype = os.path.basename(os.path.dirname(path))
    if datatype not in SUPPORTED_TYPES:
        return None

    tokens = stem.split("_")
    if len(tokens) < 2:
        return None
    suffix = tokens[-1]

    entities = {}
    for tok in tokens[:-1]:
        m = ENTITY_RE.match(tok)
        if m:
            entities[m.group(1)] = m.group(2)

    if "sub" not in entities:
        return None

    family_parts = [f"{k}-{v}" for k, v in entities.items() if k not in MEMBER_ENTITIES and k not in ("sub", "ses")]
    family = "_".join(family_parts + [suffix])

    return {
        "path": path,
        "filename": filename,
        "sub": entities["sub"],
        "ses": entities.get("ses", "NA"),
        "datatype": datatype,
        "run": entities.get("run"),
        "part": entities.get("part"),
        "family": family,
    }


###############################################################################
# GROUP
###############################################################################

def build_families(parsed_files):
    groups = defaultdict(list)
    for pf in parsed_files:
        key = (pf["sub"], pf["ses"], pf["datatype"], pf["family"])
        groups[key].append(pf)
    return groups


###############################################################################
# RULES
###############################################################################

def duplicate_rule(files):
    findings = []
    for filename, count in Counter(f["filename"] for f in files).items():
        if count > 1:
            findings.append((ERROR, "DUPLICATE_FILENAME", f"duplicate filename: {filename} ({count} occurrences)"))
    return findings


def run_part_rule(files):
    """Everything about run/part legality lives here."""
    # runs present (as ints where possible), and the set of parts seen per run
    runs = set()
    run_parts = defaultdict(set)
    for f in files:
        run = f["run"]
        if run is not None:
            try:
                run = int(run)
            except ValueError:
                pass
            runs.add(run)
            if f["part"] is not None:
                run_parts[run].add(f["part"])

    if not runs:
        return []  # no run entity at all -> fine

    runs_sorted = sorted(runs, key=lambda r: (isinstance(r, str), r))

    # single run present -> should have been omitted
    if len(runs_sorted) == 1:
        return [(ERROR, "RUN_SHOULD_BE_OMITTED",
                 f"only run-{runs_sorted[0]} exists; a single run should be omitted from the filename")]

    # must be contiguous starting at 1 (e.g. run-1,2 or run-1,2,3 — not run-1,3)
    expected = list(range(1, len(runs_sorted) + 1))
    if runs_sorted != expected:
        return [(ERROR, "RUN_NONCONTIGUOUS",
                 f"non-contiguous or missing run(s): found run-{', run-'.join(str(r) for r in runs_sorted)}")]

    # contiguous, 2+ runs: is run being used to split mag/phase?
    mag_only = [r for r, parts in run_parts.items() if parts == {"mag"}]
    phase_only = [r for r, parts in run_parts.items() if parts == {"phase"}]
    if mag_only and phase_only:
        return [(ERROR, "RUN_SEPARATES_MAG_PHASE",
                 "run incorrectly separates magnitude and phase; "
                 "this is a single acquisition — remove the run entity "
                 "(e.g. part-mag / part-phase with no run-)")]

    # every run has both mag and phase -> genuine repeat, but verify
    complete_runs = [r for r, parts in run_parts.items() if {"mag", "phase"} <= parts]
    if run_parts and len(complete_runs) == len(runs_sorted):
        return [(WARNING, "REPEATED_ACQUISITION",
                 "repeated acquisition detected (each run contains both part-mag "
                 "and part-phase); verify run labels against the acquisition "
                 "protocol or DICOM metadata")]

    # contiguous runs, no part entity involved (or partial data) -> fine
    return []


def validate_family(files):
    findings = []
    findings += duplicate_rule(files)
    findings += run_part_rule(files)
    return findings


###############################################################################
# OUTPUT
###############################################################################

def write_tsv(path, fieldnames, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dataset_root")
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    files = scan_dataset(args.dataset_root)
    parsed = [p for p in (parse_filename(f) for f in files) if p]
    groups = build_families(parsed)

    summary_rows, detail_rows = [], []
    for (sub, ses, datatype, family), group_files in groups.items():
        findings = validate_family(group_files)
        if findings:
            severity = max(findings, key=lambda x: SEVERITY_ORDER[x[0]])[0]
            message = "; ".join(msg for _, _, msg in findings)
        else:
            severity, message = VALID, "complete"

        summary_rows.append({
            "sub": sub, "ses": ses, "type": datatype,
            "acquisition_family": family, "eval": severity, "message": message,
        })
        for sev, rule, msg in findings:
            detail_rows.append({
                "sub": sub, "ses": ses, "type": datatype,
                "acquisition_family": family, "severity": sev, "rule": rule, "message": msg,
            })

    summary_rows.sort(key=lambda r: (r["sub"], r["ses"], r["type"], r["acquisition_family"]))
    detail_rows.sort(key=lambda r: (r["sub"], r["ses"], r["type"], r["acquisition_family"]))

    os.makedirs(args.outdir, exist_ok=True)
    write_tsv(os.path.join(args.outdir, "summary.tsv"),
              ["sub", "ses", "type", "acquisition_family", "eval", "message"], summary_rows)
    write_tsv(os.path.join(args.outdir, "details.tsv"),
              ["sub", "ses", "type", "acquisition_family", "severity", "rule", "message"], detail_rows)

    counts = Counter(r["eval"] for r in summary_rows)
    print(f"Families scanned: {len(summary_rows)}  "
          f"VALID={counts.get(VALID,0)} WARNING={counts.get(WARNING,0)} ERROR={counts.get(ERROR,0)}")


if __name__ == "__main__":
    sys.exit(main())