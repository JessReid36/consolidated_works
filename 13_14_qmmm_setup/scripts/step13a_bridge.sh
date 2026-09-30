#!/usr/bin/env bash
# Step 13a - generate the ORCAFF force-field bridge from the AMBER prmtop (once; shared
# by all QM/MM jobs). ORCA 6.0.1 orca_mm -convff -AMBER. The AMBER file-arg mapping is
# ambiguous in the usage text, so this tries the documented forms and reports which
# produces the .ORCAFF.prms bridge. Fast, local (login node), read-mostly.
set -uo pipefail
ORCA=/home/apps2/ORCA/6.0.1
export LD_LIBRARY_PATH="$ORCA/lib:${LD_LIBRARY_PATH:-}"
root="$HOME/system_development"
prmtop="$root/03_amber/tleap_build/complex_solvated.prmtop"
work="$root/05_qmmm/13_bridge"
[[ -s "$prmtop" ]] || { echo "FAIL missing prmtop"; exit 1; }
rm -rf "$work"; mkdir -p "$work"; cd "$work"
cp "$prmtop" ./complex_solvated.prmtop
p=complex_solvated.prmtop

try () {  # label + args...
  local label="$1"; shift
  echo "=================================================================="
  echo ">>> attempt: $label"
  echo ">>> cmd: orca_mm -convff -AMBER $*"
  "$ORCA/orca_mm" -convff -AMBER "$@" > "convff_${label}.log" 2>&1
  echo "    exit=$?"
  echo "    --- tail of log ---"; tail -8 "convff_${label}.log" | sed 's/^/    /'
  echo "    --- ORCAFF/prms produced? ---"; ls -la *.ORCAFF.prms *.prms 2>/dev/null | sed 's/^/    /' || echo "    (none)"
}

echo "ORCA: $ORCA   version dir present: $([[ -x "$ORCA/orca_mm" ]] && echo yes)"
try "A_single"   "$p"
[[ -n "$(ls *.ORCAFF.prms 2>/dev/null)" ]] || try "B_double" "$p" "$p"
[[ -n "$(ls *.ORCAFF.prms 2>/dev/null)" ]] || try "C_verbose" -verbose "$p"

echo "=================================================================="
echo "RESULT: ORCAFF files now present:"
ls -la *.ORCAFF.prms 2>/dev/null || echo "  (none — need a different input form; the logs above show ORCA's complaint)"
echo "STEP 13a DONE"
