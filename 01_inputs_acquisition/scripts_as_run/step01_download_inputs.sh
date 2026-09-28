#!/usr/bin/env bash
# Step 01 — acquire and provenance-lock raw inputs for system_development.
#   Protein scaffold : B. subtilis chorismate mutase, PDB 2CHT (chains A-L;
#                      each active site holds a transition-state analogue, TSA,
#                      NOT the chorismate substrate).
#   Terminal-repair template : PDB 1DBF.
#   Chorismate substrate coords : LigA-C from the CP2K QM/MM example set.
set -euo pipefail

root="$HOME/system_development"
structdir="$root/01_inputs/structures"
ligdir="$root/01_inputs/ligands"
admindir="$root/00_admin"
mkdir -p "$structdir" "$ligdir" "$root/01_inputs/papers" "$admindir"

fetch() {  # url dest
  curl -fL --retry 3 --retry-delay 2 -o "$2" "$1"
  [[ -s "$2" ]] || { echo "empty download: $2" >&2; exit 1; }
}

fetch https://files.rcsb.org/download/2CHT.pdb "$structdir/2cht_raw.pdb"
fetch https://files.rcsb.org/download/1DBF.pdb "$structdir/1dbf_raw.pdb"

# Ligand geometries come from a tutorial repo, so record the exact commit used;
# the SHA-256 checksums below are the primary identity anchor.
tmp="$(mktemp -d)"; trap 'rm -rf "$tmp"' EXIT
git -c pack.threads=1 clone --depth 1 \
  https://github.com/cp2k/cp2k-examples.git "$tmp/cp2k-examples"
cp2k_commit="$(git -C "$tmp/cp2k-examples" rev-parse HEAD)"
for f in LigA LigB LigC; do
  src="$tmp/cp2k-examples/qm_mm/$f.mol2"
  [[ -s "$src" ]] || { echo "missing in cp2k repo: $src" >&2; exit 1; }
  cp -f "$src" "$ligdir/$(tr '[:upper:]' '[:lower:]' <<<"$f").mol2"
done

# Structural inventory: confirm chain set and identify bound heteros.
# TSA marks the active-site location; it is the crystallographic TS analogue,
# distinct from the chorismate reactant introduced later.
for f in "$structdir/2cht_raw.pdb" "$structdir/1dbf_raw.pdb"; do
  a=$(awk '/^ATOM  /{n++} END{print n+0}' "$f")
  h=$(awk '/^HETATM/{n++} END{print n+0}' "$f")
  chains=$(awk '/^(ATOM  |HETATM)/{c=substr($0,22,1); if(c!=" ") print c}' "$f" \
             | sort -u | tr -d '\n')
  printf "%s  ATOM=%s HETATM=%s chains=%s\n" "$(basename "$f")" "$a" "$h" "$chains"
done

echo "--- 2CHT non-water heteros (residue chain resid : count) ---"
awk '/^HETATM/ && substr($0,18,3)!="HOH"{
  key=substr($0,18,3)" "substr($0,22,1)" "substr($0,23,4); n[key]++
} END{for(k in n) print k, ":", n[k]}' "$structdir/2cht_raw.pdb" | sort

# Ligand check: molecule name, atom and bond counts. 24 atoms is consistent with
# the chorismate dianion (C10H8O6) — confirm element composition and net charge
# from the mol2 before it sets the QM-region charge.
for f in liga ligb ligc; do
  awk -v file="$f" '
    /^@<TRIPOS>MOLECULE/{getline name}
    /^@<TRIPOS>ATOM/{sec="atom"; next}
    /^@<TRIPOS>BOND/{sec="bond"; next}
    /^@<TRIPOS>/{sec=""}
    sec=="atom" && NF>=6{a++}
    sec=="bond" && NF>=4{b++}
    END{printf "%s  name=%s atoms=%d bonds=%d\n", file, name, a+0, b+0}
  ' "$ligdir/$f.mol2"
done

# Provenance: role manifest, checksums, retrieval date, cp2k commit.
{
  printf "role\tpath\n"
  printf "bs_cm_scaffold\t%s\n"            "$structdir/2cht_raw.pdb"
  printf "terminal_repair_template\t%s\n"  "$structdir/1dbf_raw.pdb"
  printf "chorismate_ref_a\t%s\n"          "$ligdir/liga.mol2"
  printf "chorismate_ref_b\t%s\n"          "$ligdir/ligb.mol2"
  printf "chorismate_ref_c\t%s\n"          "$ligdir/ligc.mol2"
} > "$admindir/step01_download_manifest.tsv"

sha256sum "$structdir"/2cht_raw.pdb "$structdir"/1dbf_raw.pdb \
          "$ligdir"/lig{a,b,c}.mol2 > "$admindir/sha256_step01_inputs.txt"

{
  date -u +"retrieved_utc=%Y-%m-%dT%H:%M:%SZ"
  printf "cp2k_examples_commit=%s\n" "$cp2k_commit"
} > "$admindir/step01_provenance.txt"

echo "STEP 01 DONE"
