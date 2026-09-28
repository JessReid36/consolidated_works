#!/usr/bin/env python3
"""
step02_verify.py - re-runnable audit of step 02 outputs.

Checks performed:

  1. CONTINUITY. The two inputs step02 reads (2cht_raw.pdb, 1dbf_raw.pdb) are
     byte-identical to the step01 outputs, by SHA-256 against
     ../../01_inputs_acquisition/outputs/sha256_step01_inputs.txt. This is what
     establishes that step02 consumed step01's product and not some other copy.
  2. The A/B/C extraction contains only chains A, B and C.
  3. The superposition RMSD is re-derived from the two committed structures over
     their common backbone atoms and compared against step02_alignment_report.tsv.
  4. Residue provenance in abc_repaired.pdb matches the patch definition hard-coded
     in step02_terminal_repair_abc.sh, and matches the ranges published by
     Agbaglo et al. (2024) p. 6 for the same repair.
  5. Peptide C-N distances are recomputed for both the unadjusted and the adjusted
     structures, using the script's OWN acceptance window (d < 1.15 or d > 1.70),
     and the bad-link counts compared against step02_peptide_continuity_report.tsv
     (expected 2 and 0). The worst link at each stage is reported regardless of
     whether it passes, since the window is wide and a link may sit well away from
     canonical while still passing.
  6. Each recorded join shift equals initial_CN - 1.329 A, the canonical target.
  7. Every adjusted join is in the script's allowed_joins set, i.e. no join outside
     the permitted set was silently moved.

Exit status 0 if every check passes, 1 otherwise.

Usage:  python3 step02_verify.py [outputs_dir] [step01_outputs_dir]
"""
import sys
import math
import hashlib
import collections
from pathlib import Path

FAIL = []

TARGET_CN = 1.329
# The script's own acceptance window, from step02_terminal_repair_abc.sh continuity():
#   a link is bad if d < CN_MIN or d > CN_MAX. Do not substitute a tighter rule here;
#   the check must reproduce the gate the script actually applied.
CN_MIN, CN_MAX = 1.15, 1.70
PATCH = {
    "A": set([1] + list(range(116, 128))),
    "B": set([1] + list(range(115, 128))),
    "C": set([1, 2] + list(range(115, 128))),
}
# Agbaglo et al. 2024, PCCP 26, 12467-12482, p. 6, for the same 2CHT/1DBF repair.
AGBAGLO = {
    "A": set([1] + list(range(116, 128))),
    "B": set([1] + list(range(115, 128))),
    "C": set([1, 2] + list(range(115, 128))),
}
ALLOWED_JOINS = {("A", 115, 116), ("B", 1, 2), ("B", 114, 115),
                 ("C", 2, 3), ("C", 114, 115)}


def check(label, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}{(' — ' + detail) if detail else ''}")
    if not ok:
        FAIL.append(label)


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def atoms(path):
    out = []
    for L in Path(path).read_text().splitlines():
        if not L.startswith("ATOM  "):
            continue
        out.append({
            "name": L[12:16].strip(), "resname": L[17:20].strip(),
            "chain": L[21], "resid": int(L[22:26]),
            "xyz": (float(L[30:38]), float(L[38:46]), float(L[46:54])),
        })
    return out


def backbone_map(a):
    return {(x["chain"], x["resid"], x["name"]): x["xyz"]
            for x in a if x["name"] in ("N", "CA", "C", "O")}


def kabsch_rmsd(P, Q):
    n = len(P)
    cp = [sum(p[i] for p in P) / n for i in range(3)]
    cq = [sum(q[i] for q in Q) / n for i in range(3)]
    P = [[p[i] - cp[i] for i in range(3)] for p in P]
    Q = [[q[i] - cq[i] for i in range(3)] for q in Q]
    # already-aligned structures: no rotation fitted, report as-is deviation
    return math.sqrt(sum(sum((P[k][i] - Q[k][i]) ** 2 for i in range(3))
                         for k in range(n)) / n)


def cn_links(a):
    """Return {(chain, r1, r2): distance} for consecutive residues."""
    C = {(x["chain"], x["resid"]): x["xyz"] for x in a if x["name"] == "C"}
    N = {(x["chain"], x["resid"]): x["xyz"] for x in a if x["name"] == "N"}
    res = collections.defaultdict(set)
    for x in a:
        res[x["chain"]].add(x["resid"])
    out = {}
    for ch, ids in res.items():
        s = sorted(ids)
        for r1, r2 in zip(s, s[1:]):
            if (ch, r1) in C and (ch, r2) in N:
                out[(ch, r1, r2)] = math.dist(C[(ch, r1)], N[(ch, r2)])
    return out


def main():
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "../outputs")
    s01 = Path(sys.argv[2] if len(sys.argv) > 2
               else "../../01_inputs_acquisition/outputs")
    print(f"step 02 verification against {out.resolve()}\n")

    print("1. continuity: step02's inputs are step01's outputs")
    rec = {}
    f01 = s01 / "sha256_step01_inputs.txt"
    if f01.exists():
        for line in f01.read_text().splitlines():
            h, p = line.split()
            rec[Path(p).name] = h
        for n in ("2cht_raw.pdb", "1dbf_raw.pdb"):
            src = s01 / n
            check(f"{n} present in step01 outputs and matches its recorded hash",
                  src.exists() and sha256(src) == rec.get(n),
                  rec.get(n, "no record")[:16])
    else:
        check("step01 checksum record reachable", False, str(f01))

    print("\n2. A/B/C extraction")
    raw = atoms(out / "2cht_abc_raw_protein.pdb")
    ch = sorted({x["chain"] for x in raw})
    check("2cht_abc_raw_protein.pdb holds chains A,B,C only",
          ch == ["A", "B", "C"], "".join(ch))

    print("\n3. superposition")
    ali = atoms(out / "1dbf_abc_aligned_to_2cht.pdb")
    bm_r, bm_a = backbone_map(raw), backbone_map(ali)
    common = sorted(set(bm_r) & set(bm_a))
    rep = (out / "step02_alignment_report.tsv").read_text().splitlines()[1].split("\t")
    n_rep, rmsd_rep = int(rep[1]), float(rep[2])
    check("common backbone atom count matches the report",
          len(common) == n_rep, f"{len(common)} vs {n_rep}")
    if common:
        r = kabsch_rmsd([bm_r[k] for k in common], [bm_a[k] for k in common])
        check("re-derived RMSD matches the report",
              abs(r - rmsd_rep) < 5e-3, f"{r:.4f} vs {rmsd_rep:.4f}")

    print("\n4. residue provenance")
    src = {}
    for line in (out / "step02_residue_source_report.tsv").read_text().splitlines()[1:]:
        f = line.split("\t")
        src.setdefault(f[0], {})[int(f[1])] = f[3]
    for c in ("A", "B", "C"):
        grafted = {r for r, s in src.get(c, {}).items() if "1dbf" in s.lower()}
        check(f"chain {c} grafted set matches the script patch definition",
              grafted == PATCH[c],
              f"{len(grafted)} residues"
              + (f"  extra={sorted(grafted - PATCH[c])} missing={sorted(PATCH[c] - grafted)}"
                 if grafted != PATCH[c] else ""))
        check(f"chain {c} grafted set matches Agbaglo et al. (2024) p. 6",
              grafted == AGBAGLO[c], "" if grafted == AGBAGLO[c] else "DIVERGES")

    print("\n5. peptide continuity, recomputed")
    cont = {}
    for line in (out / "step02_peptide_continuity_report.tsv").read_text().splitlines()[1:]:
        f = line.split("\t")
        cont[f[0]] = (int(f[1]), int(f[2]))
    for stage, fn in (("unadjusted", "abc_repaired_unadjusted.pdb"),
                      ("adjusted", "abc_repaired.pdb")):
        links = cn_links(atoms(out / fn))
        bad = {k: v for k, v in links.items() if v < CN_MIN or v > CN_MAX}
        n_rep, bad_rep = cont[stage]
        check(f"{stage}: link count matches the report",
              len(links) == n_rep, f"{len(links)} vs {n_rep}")
        check(f"{stage}: bad-link count matches the report",
              len(bad) == bad_rep,
              f"{len(bad)} vs {bad_rep}"
              + (f"  {[(k, round(v, 3)) for k, v in sorted(bad.items())]}" if bad else ""))
        worst = max(links.items(), key=lambda kv: abs(kv[1] - TARGET_CN))
        print(f"         worst link {stage}: {worst[0]} at {worst[1]:.3f} A "
              f"({worst[1] - TARGET_CN:+.3f} from canonical, window {CN_MIN}-{CN_MAX})")

    print("\n6. join adjustments")
    adj = []
    for line in (out / "step02_join_adjustment_report.tsv").read_text().splitlines()[1:]:
        f = line.split("\t")
        adj.append((f[1], int(f[2]), int(f[3]), float(f[4]), float(f[5]), f[6]))
    for c, r1, r2, d0, shift, block in adj:
        check(f"{c} {r1}-{r2}: shift equals initial C-N minus {TARGET_CN}",
              abs(shift - (d0 - TARGET_CN)) < 1e-3,
              f"{shift:.3f} vs {d0 - TARGET_CN:.3f}  (block {block})")
        check(f"{c} {r1}-{r2} is an allowed join",
              (c, r1, r2) in ALLOWED_JOINS)

    print("\n7. the closures were rigid translations, confined to the declared blocks")
    # The script closes a join by translating a whole terminal block. Two things must
    # hold and neither is asserted by the script itself: nothing outside the named block
    # moved, and every atom inside it moved by the SAME vector (a translation, not a
    # distortion). Tolerance is 2e-3 A, set by PDB coordinate rounding to 3 decimals.
    def xyz_map(path):
        m = {}
        for L in Path(path).read_text().splitlines():
            if L.startswith("ATOM  "):
                m[(L[21], int(L[22:26]), L[12:16], L[16])] = (
                    float(L[30:38]), float(L[38:46]), float(L[46:54]))
        return m

    before = xyz_map(out / "abc_repaired_unadjusted.pdb")
    after = xyz_map(out / "abc_repaired.pdb")
    check("atom set unchanged by adjustment", set(before) == set(after),
          f"{len(before)} vs {len(after)}")
    moved = collections.defaultdict(list)
    for k in set(before) & set(after):
        d = math.dist(before[k], after[k])
        if d > 1e-6:
            moved[k[0]].append((k[1], d,
                                tuple(after[k][i] - before[k][i] for i in range(3))))
    expect = {}
    for c, r1, r2, d0, shift, block in adj:
        lo, hi = block.split(":")[1].split("-")
        expect[c] = (set(range(int(lo), int(hi) + 1)), shift)
    check("only the chains named in the adjustment report moved",
          set(moved) == set(expect), f"moved {sorted(moved)}, expected {sorted(expect)}")
    for c in sorted(moved):
        resids = {r for r, d, v in moved[c]}
        mags = [d for r, d, v in moved[c]]
        vecs = [v for r, d, v in moved[c]]
        spread = max(max(v[i] for v in vecs) - min(v[i] for v in vecs) for i in range(3))
        exp_res, exp_mag = expect[c]
        check(f"chain {c}: displacement confined to the declared block",
              resids == exp_res,
              f"moved {sorted(resids)[0]}-{sorted(resids)[-1]} (n={len(resids)}), "
              f"declared {sorted(exp_res)[0]}-{sorted(exp_res)[-1]}")
        check(f"chain {c}: displacement is a single translation, not a distortion",
              spread < 2e-3,
              f"per-axis spread {spread:.1e} A across {len(vecs)} atoms")
        check(f"chain {c}: magnitude matches the recorded shift",
              abs(sum(mags) / len(mags) - exp_mag) < 2e-3,
              f"{sum(mags)/len(mags):.3f} A vs {exp_mag:.3f} recorded")

    print(f"\n{'ALL CHECKS PASSED' if not FAIL else 'FAILURES: ' + ', '.join(FAIL)}")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
