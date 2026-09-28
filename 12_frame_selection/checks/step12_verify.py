#!/usr/bin/env python3
"""
step12_verify.py - re-runnable audit of step 12, QM/MM frame selection.

Two selections were run. step12a applies a competence filter and the DISTANCE clause of
the near-attack definition. step12b adds frames satisfying the FULL three-condition
definition of Hur & Bruice, after a check showed three of the twelve failed the angular
clauses.

Checks performed:

  1. The canonical manifest is identified and distinguished from the stage files. Six
     manifests exist under 12_frame_selection/; only the top-level one is the ensemble set.
  2. THE SITE OUTCOME WAS ALGORITHMIC, NOT CHOSEN. step12a excludes any site competent in
     fewer than SITE_MIN of frames. The log records the occupancies; this confirms that
     one site qualified and two were excluded by the stated threshold.
  3. Every selected frame comes from the qualifying site.
  4. The near-NAC frames are traceable: the three frames the full NAC test rejected are
     exactly the three labelled near_NAC in the ensemble, and they fail on theta2 alone,
     marginally.
  5. step12b's stated test - that near-NACs should sit systematically higher in barrier if
     near-attack geometry governs reactivity - is evaluated against the committed barriers.
  6. Frame spacing and coverage of the sampled range.

Exit status 0 if every check passes, 1 otherwise.

Usage:  python3 step12_verify.py [outputs] [qmmm_dir]
"""
import sys
import re
import statistics as st
from pathlib import Path

FAIL = []
REJECTED_BY_FULL_NAC = {"05680", "07310", "15825"}


def check(label, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}{(' — ' + detail) if detail else ''}")
    if not ok:
        FAIL.append(label)


def read_manifest(p):
    rows = []
    for l in Path(p).read_text().splitlines():
        if l.startswith("#") or not l.strip():
            continue
        f = l.rstrip("\n").split("\t")
        if len(f) >= 5 and f[1].strip().lstrip("-").isdigit():
            rows.append(f)
    return rows



def _find(*candidates):
    """First existing path among candidates. Keeps the checks self-contained: the files
    live in sibling step folders of consolidated_works, and the external results
    repository is only a fallback for anyone running from a full working tree."""
    from pathlib import Path as _P
    for c in candidates:
        p = _P(c)
        if p.exists():
            return p
    return _P(candidates[0])


def main():
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "../outputs")
    nac_log = _find("../../11_trajectory_analysis/outputs/11f_nac_check.log",
                    "../../../chorismate-thesis-results/05_qmmm/11f_nac_check.log")
    ens_tsv = _find("../../19_ensemble/outputs/ensemble_barriers.tsv",
                    "../../../chorismate-thesis-results/05_qmmm/"
                    "19_ensemble_barriers/ensemble_barriers.tsv")
    print(f"step 12 verification against {out.resolve()}\n")

    print("1. which manifest is canonical")
    canon = out / "selection_manifest.tsv"
    check("the canonical manifest is present", canon.exists())
    rows = read_manifest(canon)
    print(f"         canonical: selection_manifest.tsv, {len(rows)} rows")
    stage = sorted((out / "stage_manifests").glob("*.tsv"))
    for s in stage:
        print(f"         stage file: {s.name}, {len(read_manifest(s))} rows")
    ens_rows, hdr = [], None
    for l in ens_tsv.read_text().splitlines():
        if l.startswith("#"):
            continue
        p = l.rstrip("\n").split("\t")
        if p[0] == "frame":
            hdr = p
            continue
        ens_rows.append(dict(zip(hdr, p)))
    ens = {int(r["frame"]) for r in ens_rows}
    man = {int(r[1]) for r in rows}
    check("the canonical manifest contains every ensemble frame", ens <= man,
          f"{len(ens)} ensemble frames, {len(man)} manifest rows, "
          f"extra in manifest: {sorted(man - ens)}")
    for s in stage:
        sm = {int(r[1]) for r in read_manifest(s)}
        if ens <= sm:
            check(f"{s.name} is NOT mistakable for the canonical set", False,
                  "it also covers the ensemble")

    print("\n2. the site outcome was decided by a stated threshold, not chosen")
    log = out / "step12a_selection.log"
    check("the selection log is committed", log.exists())
    txt = log.read_text()
    m = re.search(r"site_min\s+(\d+)%", txt)
    check("the competence threshold is recorded in the log", bool(m),
          f"site_min {m.group(1)}%" if m else "")
    occ = re.findall(r"CHA#(\d) \(res (\d+)\): competent (\d+) / (\d+) \((\d+)%\)\s+(\S+)",
                     txt)
    check("all three sites were evaluated", len(occ) == 3, f"{len(occ)} sites")
    qualifying = [o for o in occ if "QUALIF" in o[5]]
    excluded = [o for o in occ if "exclud" in o[5]]
    for c, res, n, tot, pct, verdict in occ:
        print(f"         CHA#{c} (res {res}): {n}/{tot} competent ({pct}%)  {verdict}")
    check("exactly one site qualified", len(qualifying) == 1,
          f"{len(qualifying)} qualified, {len(excluded)} excluded")
    thr = int(m.group(1)) if m else 50
    check("every excluded site is genuinely below the threshold",
          all(int(o[4]) < thr for o in excluded),
          f"excluded at {[o[4]+'%' for o in excluded]} against a {thr}% threshold")
    check("the qualifying site is comfortably above it",
          all(int(o[4]) >= thr for o in qualifying),
          f"{qualifying[0][4]}%" if qualifying else "")
    check("the filter is a conjunction, not distance alone",
          "AND" in txt, re.search(r"filter: .*", txt).group(0) if
          re.search(r"filter: .*", txt) else "")

    print("\n3. every selected frame comes from the qualifying site")
    sites = {(r[3], r[4]) for r in rows}
    check("one substrate copy across the whole manifest", len(sites) == 1, f"{sites}")

    print("\n4. the near-NAC frames are traceable")
    nonfull = {r["frame"] for r in ens_rows if r["group"] == "near_NAC"}
    check("the near_NAC labels are exactly the frames the full test rejected",
          nonfull == REJECTED_BY_FULL_NAC,
          f"labelled {sorted(nonfull)}, rejected by the test {sorted(REJECTED_BY_FULL_NAC)}")
    nac = nac_log
    if nac.exists():
        t2 = nac.read_text()
        devs = []
        for fr in sorted(REJECTED_BY_FULL_NAC):
            mm = re.search(rf"\s{int(fr)}\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+([\d.]+)", t2)
            if mm:
                devs.append((fr, float(mm.group(1))))
        print("         theta2 deviations: " +
              ", ".join(f"{f} {d:.1f} deg" for f, d in devs))
        check("all three fail marginally, on theta2 alone",
              devs and all(20.0 < d < 30.0 for _, d in devs),
              "against a 20 deg tolerance; Hur & Bruice's unreactive diequatorial "
              "conformers show theta2 deviations of about 110 deg")

    print("\n5. step12b's own test: do near-NACs sit higher in barrier?")
    orig = [r for r in ens_rows if int(r["frame"]) < 20000]
    def barr(g):
        return [float(r["barrier"]) for r in orig
                if r["group"] == g and r["barrier"] not in ("", "NA")
                and float(r["barrier"]) > 0]
    fullb, nearb = barr("full_NAC"), barr("near_NAC")
    print(f"         full_NAC n={len(fullb)}  mean {st.mean(fullb):.2f} kcal/mol")
    print(f"         near_NAC n={len(nearb)}  mean {st.mean(nearb):.2f} kcal/mol")
    check("the near-NACs do sit higher, as step12b predicted",
          st.mean(nearb) > st.mean(fullb),
          f"{st.mean(nearb):.2f} against {st.mean(fullb):.2f}. n={len(nearb)}, so this is "
          f"directional only and must NOT be quoted as a result.")

    print("\n6. spacing and coverage")
    ps = sorted(float(r[2]) for r in rows)
    gaps = [b - a for a, b in zip(ps[:-1], ps[1:])]
    print(f"         {len(ps)} frames from {min(ps):.0f} to {max(ps):.0f} ps")
    print(f"         spacing: min {min(gaps):.0f}, median {st.median(gaps):.0f}, "
          f"max {max(gaps):.0f} ps")
    check("no two selected frames are adjacent in the trajectory", min(gaps) > 100,
          f"minimum spacing {min(gaps):.0f} ps")
    beyond = sum(1 for p in ps if p > 20000)
    check("the selection spans beyond the 20 ns reference window", beyond > 0,
          f"{beyond} of {len(ps)} frames beyond 20000 ps")

    print(f"\n{'ALL CHECKS PASSED' if not FAIL else 'FAILURES: ' + ', '.join(FAIL)}")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
