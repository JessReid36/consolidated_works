#!/usr/bin/env bash
# Step 08b - GAFF-type chorismate and build its frcmod on the HPC, using the
# AM1-BCC charges derived off-HPC in step 08a. No sqm here: antechamber runs in
# read-charge mode (-c rc), so this needs only the (working) HPC Antechamber/GAFF.
#
# Follows the reference protocol's Antechamber + GAFF direction. AMBER22 is used
# because AMBER18 is unavailable on the HPC; the GAFF (not gaff2) parameter set is
# selected explicitly to stay close to the paper. One canonical CHA template
# (cha_a) is used for all three ligand copies.
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

root="$HOME/system_development"
indir="$root/02_preparation/accepted_preprotonation"
outdir="$root/03_amber/ligand_gaff"
admin="$root/00_admin"
canonical="$indir/cha_a.mol2"
charges="$outdir/charges_am1bcc.dat"        # from step 08a (scp'd up)
mkdir -p "$outdir" "$admin"

echo "=== step 08b: input presence ==="
for f in "$indir/cha_a.mol2" "$indir/cha_b.mol2" "$indir/cha_c.mol2" "$charges"; do
  [[ -s "$f" ]] || { echo "FAIL missing/empty: $f"; exit 1; }
  echo "PASS $f"
done

echo
echo "=== step 08b: verify cha_a/b/c share atom names + bond topology, and charge file lines up ==="
python3 - "$indir/cha_a.mol2" "$indir/cha_b.mol2" "$indir/cha_c.mol2" "$charges" "$admin/step08b_ligand_identity_check.tsv" <<'PY'
import sys
from collections import Counter
a,b,c,chargefile,report = sys.argv[1:6]

def read(path):
    atoms=[]; bonds=[]; sec=None
    for l in open(path):
        if l.startswith("@<TRIPOS>"): sec=l.strip(); continue
        if sec=="@<TRIPOS>ATOM" and len(l.split())>=9:
            p=l.split(); atoms.append((p[1],p[5],float(p[8])))
        elif sec=="@<TRIPOS>BOND" and len(l.split())>=4:
            p=l.split(); bonds.append((p[1],p[2],p[3]))
    return atoms,bonds

A=read(a); B=read(b); C=read(c)
na=len(A[0])
if na==0: sys.exit("FAIL cha_a has no atoms")
names_a=[x[0] for x in A[0]]; bonds_a=A[1]
dup=[n for n,k in Counter(names_a).items() if k>1]
if dup: sys.exit(f"FAIL cha_a duplicate atom names: {dup}")

rows=[]
for lbl,D in (("cha_a",A),("cha_b",B),("cha_c",C)):
    nm=[x[0] for x in D[0]]
    names_match = nm==names_a
    bonds_match = D[1]==bonds_a
    net=sum(x[2] for x in D[0])
    rows.append((lbl,len(D[0]),len(D[1]),f"{net:.4f}",names_match,bonds_match))
    if not names_match: sys.exit(f"FAIL {lbl}: atom names differ from cha_a")
    if not bonds_match: sys.exit(f"FAIL {lbl}: bond topology differs from cha_a")

with open(report,"w") as f:
    f.write("label\tatoms\tbonds\tnet_charge\tnames_match_cha_a\tbonds_match_cha_a\n")
    for r in rows: f.write("\t".join(map(str,r))+"\n")

qs=[float(x) for x in open(chargefile) if x.strip()]
if len(qs)!=na: sys.exit(f"FAIL charge file has {len(qs)} values, canonical CHA has {na} atoms")
if abs(sum(qs)+2.0)>1e-4: sys.exit(f"FAIL charge file sums to {sum(qs):.4f}, expected -2")
print(f"PASS cha_a/b/c share names+bonds ({na} atoms); AM1-BCC charge file lines up (sum {sum(qs):+.4f})")
print("     note: cha_b/c charges may differ from cha_a (known); one canonical template is used")
PY

echo
echo "=== step 08b: load AMBER22 tools ==="
# module load sources amber.sh, which references several possibly-unset vars
# (PYTHONPATH, PERL5LIB, ...); relax -u just around it so set -u does not abort.
set +u
export PERL5LIB="${PERL5LIB:-}"
module load app/amber22/22
set -u
echo "AMBERHOME=${AMBERHOME:-UNSET}"
for exe in antechamber parmchk2; do
  command -v "$exe" >/dev/null 2>&1 || { echo "FAIL missing executable: $exe"; exit 1; }
  echo "PASS $exe -> $(command -v "$exe")"
done

echo
echo "=== step 08b: GAFF atom typing with AM1-BCC read charges ==="
cd "$outdir"
[[ "$charges" -ef ./charges_am1bcc.dat ]] || cp -f "$charges" ./charges_am1bcc.dat
rm -f cha_gaff.mol2 cha.frcmod antechamber_cha_gaff.log parmchk2_cha_gaff.log \
      ANTECHAMBER* ATOMTYPE.INF NEWPDB.PDB PREP.INF
if ! antechamber -i "$canonical" -fi mol2 -o cha_gaff.mol2 -fo mol2 \
        -at gaff -c rc -cf charges_am1bcc.dat -nc -2 -rn CHA \
        > antechamber_cha_gaff.log 2>&1; then
  echo "FAIL antechamber failed"; tail -20 antechamber_cha_gaff.log; exit 1
fi
[[ -s cha_gaff.mol2 ]] || { echo "FAIL no cha_gaff.mol2"; tail -20 antechamber_cha_gaff.log; exit 1; }
echo "PASS antechamber produced cha_gaff.mol2"

echo
echo "=== step 08b: parmchk2 (GAFF) ==="
if ! parmchk2 -i cha_gaff.mol2 -f mol2 -o cha.frcmod -s gaff > parmchk2_cha_gaff.log 2>&1; then
  echo "FAIL parmchk2 returned non-zero"; tail -20 parmchk2_cha_gaff.log; exit 1
fi
# NOTE: parmchk2 can succeed with an EMPTY log; validate cha.frcmod itself, not the log.
[[ -s cha.frcmod ]] || { echo "FAIL parmchk2 did not produce cha.frcmod"; tail -20 parmchk2_cha_gaff.log; exit 1; }
echo "PASS parmchk2 produced cha.frcmod"

echo
echo "=== step 08b: audit GAFF mol2 + frcmod ==="
python3 - "$canonical" cha_gaff.mol2 cha.frcmod charges_am1bcc.dat "$admin/step08b_cha_gaff_audit.tsv" <<'PY'
import sys
from collections import Counter
canon, gaff, frcmod, chargefile, audit = sys.argv[1:6]

def read(path):
    atoms=[]; bonds=0; sec=None
    for l in open(path):
        if l.startswith("@<TRIPOS>"): sec=l.strip(); continue
        if sec=="@<TRIPOS>ATOM" and len(l.split())>=9:
            p=l.split(); atoms.append({"name":p[1],"type":p[5],"q":float(p[8])})
        elif sec=="@<TRIPOS>BOND" and len(l.split())>=4:
            bonds+=1
    return atoms,bonds

can_atoms,_ = read(canon)
atoms,bonds = read(gaff)
na = len(can_atoms)

if len(atoms)!=na: sys.exit(f"FAIL cha_gaff atoms {len(atoms)}, expected {na}")
if bonds!=na:      sys.exit(f"FAIL cha_gaff bonds {bonds}, expected {na}")
names=[a["name"] for a in atoms]
dup=[n for n,c in Counter(names).items() if c>1]
if dup: sys.exit(f"FAIL duplicate atom names in cha_gaff: {dup}")
net=sum(a["q"] for a in atoms)
if abs(net+2.0)>1e-4: sys.exit(f"FAIL cha_gaff net charge {net:.4f}, expected -2")
types=sorted(set(a["type"] for a in atoms))
dot=[t for t in types if "." in t]
if dot: sys.exit(f"FAIL Tripos dot-types remain: {dot}")

qfile=[float(x) for x in open(chargefile) if x.strip()]
maxdiff=max(abs(a["q"]-q) for a,q in zip(atoms,qfile))
if maxdiff>1e-4: sys.exit(f"FAIL read-charge did not preserve AM1-BCC charges (max diff {maxdiff:.6f})")

frc=open(frcmod,errors="replace").read()
attn=frc.count("ATTN"); missing=frc.upper().count("MISSING")

with open(audit,"w") as f:
    f.write("item\tvalue\n")
    for k,v in [("atoms",len(atoms)),("bonds",bonds),("net_charge",f"{net:.6f}"),
                ("duplicate_atom_names","none"),("gaff_atom_types",",".join(types)),
                ("tripos_dot_types_remaining","none"),
                ("max_charge_diff_vs_am1bcc",f"{maxdiff:.8f}"),
                ("charge_model","AM1-BCC (derived off-HPC in step 08a)"),
                ("frcmod_ATTN_count",attn),("frcmod_MISSING_word_count",missing)]:
        f.write(f"{k}\t{v}\n")

if attn or missing:
    print(f"REVIEW frcmod has ATTN={attn} MISSING={missing} (inspect cha.frcmod)")
print(f"PASS atoms={len(atoms)} bonds={bonds} net={net:+.6f} types={','.join(types)}")
print(f"     charges preserved vs AM1-BCC (max diff {maxdiff:.6f}); frcmod ATTN={attn} MISSING={missing}")
PY

echo
echo "=== step 08b: method note ==="
cat > "$admin/step08_method_note.txt" <<'EOF'
Step 08 method note (chorismate GAFF parameterisation):

Reference protocol:
  Chorismate parameterised with Antechamber + GAFF before AMBER MD.

Charge model (deviation recorded):
  Charges are AM1-BCC - the standard Antechamber/GAFF scheme, and the method
  implied by the protocol (which states "Antechamber/GAFF" without naming the
  charge model). AM1-BCC needs sqm, which cannot run on this HPC
  (missing libopenblas.so.0). Charges were therefore derived AM1-BCC on a local
  Ubuntu workstation with a self-contained conda AmberTools (step 08a), then read
  in here via antechamber -c rc -cf charges_am1bcc.dat. GAFF atom typing and
  parmchk2 run natively on the HPC AMBER22 (no sqm required).
  Provenance of the derivation: 00_admin/step08a_provenance.txt.

Implementation:
  AMBER18 unavailable -> AMBER22 tools; GAFF (not gaff2) selected explicitly.
  One canonical CHA template (cha_a) used for all three ligand copies; cha_b/cha_c
  share atom names and bond topology (their charges may differ but are not used).
    antechamber -i cha_a.mol2 -fi mol2 -o cha_gaff.mol2 -fo mol2 \
                -at gaff -c rc -cf charges_am1bcc.dat -nc -2 -rn CHA
    parmchk2 -i cha_gaff.mol2 -f mol2 -o cha.frcmod -s gaff
EOF
cat "$admin/step08_method_note.txt"

echo
echo "=== step 08b: checksums ==="
sha256sum "$outdir/charges_am1bcc.dat" "$outdir/cha_gaff.mol2" "$outdir/cha.frcmod" \
  > "$admin/sha256_step08_ligand_gaff.txt"
cat "$admin/sha256_step08_ligand_gaff.txt"

echo
echo "=== step 08b: outputs ==="
for f in "$outdir/cha_gaff.mol2" "$outdir/cha.frcmod" \
         "$admin/step08b_ligand_identity_check.tsv" "$admin/step08b_cha_gaff_audit.tsv" \
         "$admin/step08_method_note.txt"; do
  [[ -s "$f" ]] || { echo "FAIL missing/empty output: $f"; exit 1; }
  echo "PASS $f"
done
echo
echo "STEP 08b DONE - cha_gaff.mol2 + cha.frcmod ready for tleap (step 09)"
