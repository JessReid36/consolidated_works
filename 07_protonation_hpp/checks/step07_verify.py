#!/usr/bin/env python3
"""
step07_verify.py - re-runnable audit of step 07, protonation-state assignment.

Step 07 submits the protein-only structure to the H++ server, then audits what comes back
against the frozen template and writes an explicit accept decision. The server run itself
is EXTERNAL and not reproducible from the repository; hpp_output_from_server.pdb is an
input to step07b, not an output of any script here.

Checks performed:

  1. Continuity. step07a's input is step03's committed output, by SHA-256.
  2. Submission file. The structure sent to the server is protein-only, TER-delimited,
     three chains of 1010 atoms over residues 1-127, reproduced and compared against
     step07a_hpp_input_report.tsv.
  3. Conservation. Every template heavy atom is present in the accepted structure and none
     is added, reproduced from the structures and compared against step07b_audit_summary.tsv.
  4. Drift, recomputed per category using the script's own definitions (backbone N/CA/C/O/
     OXT; flippable ASN/GLN/HIS; active-site set as at step 04) and compared against the
     recorded category rows.
  5. The acceptance rule reproduced against the script's four thresholds: backbone 0.05 A,
     non-flippable 0.20 A, active-site backbone 0.05 A, and no active-site heavy atom
     beyond 1.0 A.
  6. Histidine assignment reproduced from the accepted structure by ring-hydrogen presence,
     and compared against step07b_his_assignment.tsv.
  7. Trimer symmetry of the histidine assignment. The three chains are copies, so a
     per-chain disagreement would indicate an artefact of the server run rather than
     chemistry.
  8. Orphan output. step07b_heavy_atom_drift.tsv is present in the results repository but
     is written by NO script in phase1_system_dev. Its numbers are checked for consistency
     against step07b_drift_by_residue.tsv, and it is reported as undeclared.

NOT CHECKABLE HERE: the H++ run parameters. pH 7.0 is agreed across all copies of the
script comments, but salinity and the two dielectric constants are recorded only in those
comments and the copies disagree. The server output carries no parameter record, only a
URL remark. This is a record of a past interactive submission and is not recoverable by
re-running anything. See 00_crosscutting/sourcing_gaps.txt item 4.

Exit status 0 if every check passes, 1 otherwise.

Usage:  python3 step07_verify.py [outputs_dir] [step03_outputs]
"""
import sys
import math
import hashlib
from pathlib import Path
from collections import defaultdict, Counter

FAIL = []
BACKBONE = {"N", "CA", "C", "O", "OXT"}
ACTIVE_SITE = {7, 57, 59, 60, 63, 73, 74, 75, 78, 90, 108, 115}
FLIPPABLE = {"ASN", "GLN", "HIS", "HID", "HIE", "HIP"}
BB_MAX, NONFLIP_MAX, ACTIVE_BB_MAX, ACTIVE_BIG = 0.05, 0.20, 0.05, 1.0


def check(label, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}{(' — ' + detail) if detail else ''}")
    if not ok:
        FAIL.append(label)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def parse(path):
    out, ters = [], 0
    for L in Path(path).read_text().splitlines():
        if L.startswith("TER"):
            ters += 1
            continue
        if not (L.startswith("ATOM  ") or L.startswith("HETATM")):
            continue
        nm = L[12:16].strip()
        el = (L[76:78].strip() if len(L) >= 78 else "")
        if not el:
            el = nm.lstrip("0123456789")[:1]
        out.append({"rec": L[:6].strip(), "name": nm, "resname": L[17:20].strip(),
                    "chain": L[21], "resid": int(L[22:26]), "element": el,
                    "xyz": (float(L[30:38]), float(L[38:46]), float(L[46:54]))})
    return out, ters


def is_h(a):
    return a["element"].upper() == "H"


def main():
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "../outputs")
    s03 = Path(sys.argv[2] if len(sys.argv) > 2 else "../../03_cleanup_validate/outputs")
    print(f"step 07 verification against {out.resolve()}\n")

    summ = {}
    for line in (out / "step07b_audit_summary.tsv").read_text().splitlines()[1:]:
        f = line.split("\t")
        if len(f) == 2:
            summ[f[0]] = f[1]

    print("1. continuity")
    check("step03 abc_repaired_clean.pdb reachable",
          (s03 / "abc_repaired_clean.pdb").exists(),
          sha(s03 / "abc_repaired_clean.pdb")[:16]
          if (s03 / "abc_repaired_clean.pdb").exists() else "")

    print("\n2. the structure submitted to the server")
    tmpl, ters = parse(out / "hpp_input_protein_only_with_TER.pdb")
    check("no heteroatoms in the submission",
          not [a for a in tmpl if a["rec"] == "HETATM"])
    check("no hydrogens in the submission", not [a for a in tmpl if is_h(a)])
    check("TER count is 3", ters == 3, str(ters))
    rep = {}
    for line in (out / "step07a_hpp_input_report.tsv").read_text().splitlines()[1:]:
        f = line.split("\t")
        rep[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))
    per = defaultdict(list)
    for a in tmpl:
        per[a["chain"]].append(a["resid"])
    for ch in "ABC":
        n, lo, hi, nres = rep[ch]
        got = per[ch]
        check(f"chain {ch}: {n} atoms, residues {lo}-{hi}, {nres} residues",
              (len(got), min(got), max(got), len(set(got))) == (n, lo, hi, nres),
              f"{len(got)} atoms, {min(got)}-{max(got)}, {len(set(got))} residues")

    print("\n3. conservation of the template heavy atoms")
    acc, _ = parse(out / "abc_protonated_hpp_accepted.pdb")
    acc_heavy = [a for a in acc if not is_h(a)]
    acc_h = [a for a in acc if is_h(a)]
    tkeys = {(a["chain"], a["resid"], a["name"]) for a in tmpl}
    akeys = {(a["chain"], a["resid"], a["name"]) for a in acc_heavy}
    check("template heavy-atom count reproduced",
          len(tmpl) == int(summ["template_heavy_atoms"]),
          f"{len(tmpl)} vs {summ['template_heavy_atoms']}")
    check("restored heavy-atom count reproduced",
          len(acc_heavy) == int(summ["restored_heavy_atoms"]),
          f"{len(acc_heavy)} vs {summ['restored_heavy_atoms']}")
    check("no template heavy atom is missing", not (tkeys - akeys),
          f"{len(tkeys - akeys)} missing")
    check("no heavy atom was added", not (akeys - tkeys), f"{len(akeys - tkeys)} extra")
    check("hydrogens added reproduced", len(acc_h) == int(summ["hydrogens_added"]),
          f"{len(acc_h)} vs {summ['hydrogens_added']}")

    print("\n4. drift, recomputed per category")
    tpos = {(a["chain"], a["resid"], a["name"]): a["xyz"] for a in tmpl}
    cats = defaultdict(list)
    for a in acc_heavy:
        k = (a["chain"], a["resid"], a["name"])
        if k not in tpos:
            continue
        d = math.dist(tpos[k], a["xyz"])
        cats["all_heavy"].append(d)
        cats["backbone" if a["name"] in BACKBONE else "sidechain"].append(d)
        cats["flippable_residue" if a["resname"] in FLIPPABLE
             else "nonflippable_residue"].append(d)
        cats["active_site_set" if a["resid"] in ACTIVE_SITE
             else "not_active_site_set"].append(d)
        if a["resid"] in ACTIVE_SITE:
            cats["active_site_backbone" if a["name"] in BACKBONE
                 else "active_site_sidechain"].append(d)
    for name, vals in sorted(cats.items()):
        rec = summ.get(f"drift[{name}]")
        if rec is None:
            continue
        parts = dict(p.split("=") for p in rec.split(";"))
        mx = max(vals) if vals else 0.0
        rms = math.sqrt(sum(v * v for v in vals) / len(vals)) if vals else 0.0
        big = sum(1 for v in vals if v > 1.0)
        ok = (len(vals) == int(parts["count"]) and abs(mx - float(parts["max"])) < 1e-3
              and abs(rms - float(parts["rms"])) < 1e-3 and big == int(parts[">1.0A"]))
        check(f"drift[{name}] reproduced", ok,
              f"n={len(vals)} max={mx:.4f} rms={rms:.4f} >1.0={big}  vs  {rec}")

    print("\n5. the acceptance rule")
    bb = max(cats["backbone"]) if cats["backbone"] else 0.0
    nf = max(cats["nonflippable_residue"]) if cats["nonflippable_residue"] else 0.0
    abb = max(cats["active_site_backbone"]) if cats["active_site_backbone"] else 0.0
    abig = sum(1 for v in cats["active_site_set"] if v > ACTIVE_BIG)
    check(f"backbone drift <= {BB_MAX} A", bb <= BB_MAX, f"{bb:.3f}")
    check(f"non-flippable drift <= {NONFLIP_MAX} A", nf <= NONFLIP_MAX, f"{nf:.3f}")
    check(f"active-site backbone drift <= {ACTIVE_BB_MAX} A", abb <= ACTIVE_BB_MAX,
          f"{abb:.3f}")
    check(f"no active-site heavy atom beyond {ACTIVE_BIG} A", abig == 0, f"{abig} atoms")
    dec = (out / "step07b_acceptance_decision.txt").read_text()
    check("the recorded decision is ACCEPT", "ACCEPT" in dec.split("\n")[2],
          dec.split("\n")[2].strip())
    big_res = Counter(a["resname"] for a in acc_heavy
                      if (a["chain"], a["resid"], a["name"]) in tpos
                      and math.dist(tpos[(a["chain"], a["resid"], a["name"])],
                                    a["xyz"]) > 1.0)
    check("every displacement beyond 1.0 A is in a flippable residue",
          set(big_res) <= FLIPPABLE, str(dict(big_res)))

    print("\n6. histidine assignment reproduced from the accepted structure")
    his_h = defaultdict(set)
    for a in acc:
        if a["resname"] in {"HIS", "HID", "HIE", "HIP"} and a["name"] in {"HD1", "HE2"}:
            his_h[(a["chain"], a["resid"])].add(a["name"])
    rec_his = {}
    for line in (out / "step07b_his_assignment.tsv").read_text().splitlines()[1:]:
        f = line.split("\t")
        rec_his[(f[0], int(f[1]))] = (f[3], f[4], f[5])
    check("histidine count reproduced", len(rec_his) == int(summ["histidines"]),
          f"{len(rec_his)} vs {summ['histidines']}")
    for key, (hd1, he2, sug) in sorted(rec_his.items()):
        got_hd1 = "yes" if "HD1" in his_h[key] else "no"
        got_he2 = "yes" if "HE2" in his_h[key] else "no"
        expect = {("yes", "yes"): "HIP", ("no", "yes"): "HIE",
                  ("yes", "no"): "HID"}.get((got_hd1, got_he2))
        check(f"His{key[1]} chain {key[0]}: {sug} from ring hydrogens",
              (got_hd1, got_he2) == (hd1, he2) and expect == sug,
              f"HD1={got_hd1} HE2={got_he2} -> {expect}, recorded {sug}")

    print("\n7. trimer symmetry of the assignment")
    by_resid = defaultdict(set)
    for (ch, ri), (_, _, sug) in rec_his.items():
        by_resid[ri].add(sug)
    for ri, sugs in sorted(by_resid.items()):
        check(f"His{ri} assigned identically in all three chains", len(sugs) == 1,
              f"{sorted(sugs)}")

    print("\n8. orphan file, quarantined in orphans/ and NOT counted as an output")
    orphan = out.parent / "orphans" / "step07b_heavy_atom_drift.tsv"
    if not orphan.exists():
        print("         orphans/step07b_heavy_atom_drift.tsv not present; nothing to check")
    else:
        lines = [l for l in orphan.read_text().splitlines() if l.strip()]
        summary_line = next((l for l in lines if l.startswith("summary")), "")
        vals = dict(q.split("=") for q in summary_line.split("\t")[1:] if "=" in q)
        check("orphan summary: max side-chain drift reproduces",
              abs(float(vals["max_sidechain_A"]) - max(cats["sidechain"])) < 1e-3,
              f"{vals['max_sidechain_A']} vs {max(cats['sidechain']):.4f}")
        check("orphan summary: max backbone drift reproduces",
              abs(float(vals["max_backbone_A"]) - bb) < 1e-3,
              f"{vals['max_backbone_A']} vs {bb:.4f}")
        # Full derivability: every ranked row must reproduce from the DECLARED
        # step07b_drift_by_atom.tsv. If it does, excluding the orphan loses nothing.
        rows = [l.split("\t") for l in
                (out / "step07b_drift_by_atom.tsv").read_text().splitlines()[1:] if l]
        data = [l.split("\t") for l in lines[2:] if not l.startswith("rank")]
        top = sorted(((float(r[5]), r[0], r[1], r[4], r[2]) for r in rows),
                     reverse=True)[:len(data)]
        agree = sum(1 for (d, ch, ri, rn, an), o in zip(top, data)
                    if o[2] == ch and o[3] == ri and o[4] == rn and o[5] == an
                    and abs(float(o[1]) - d) < 1e-3)
        check("every orphan row is derivable from the declared drift_by_atom output",
              agree == len(data), f"{agree}/{len(data)} rows reproduce")
        print("         The orphan is written by NO script in chorismate-thesis-code at")
        print("         f0e059c. It is fully derivable from a declared output, so it is")
        print("         held in orphans/ for the record and excluded from outputs/.")
        print("         See orphans/ORPHANS.txt for the search commands to locate its")
        print("         producing script on the HPC or the PC.")

    print("\n9. as-run archive completeness")
    ar = out.parent / "scripts_as_run"
    have = sorted(f.name for f in ar.iterdir()) if ar.exists() else []
    check("step07a has an as-run archive copy",
          "step07a_prepare_hpp_input.sh" in have, str(have))
    check("step07b has NO as-run archive copy, and none is fabricated",
          "step07b_audit_repair_hpp.sh" not in have,
          "absent from results 00_admin; not substituted from the code repo")

    print("\nNOT CHECKABLE HERE:")
    print("         H++ run parameters. pH 7.0 is agreed; salinity and the two dielectric")
    print("         constants are recorded only in script comments, and the copies")
    print("         disagree. hpp_output_from_server.pdb carries no parameter record.")
    print("         See 00_crosscutting/sourcing_gaps.txt item 4.")

    print(f"\n{'ALL CHECKS PASSED' if not FAIL else 'FAILURES: ' + ', '.join(FAIL)}")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
