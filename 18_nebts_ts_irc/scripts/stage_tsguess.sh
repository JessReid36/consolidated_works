#!/usr/bin/env bash
# stage_tsguess.sh - copy the highest-energy window of a completed scan to
# tsguess.pdb, for use as the transition-state guess of the band calculation.
#
# WHY THIS IS A SEPARATE SCRIPT
#   The band needs three structures: the two optimised endpoints and a guess at
#   the saddle. The scan produces all three, and its maximum is where the
#   reference calculation's guess came from: 18_nebts/tsguess.pdb is
#   byte-identical to 16_scan/win_11.pdb, and window 11 is that scan's energy
#   maximum.
#
#   Keeping this in its own file, copied into each frame directory the way
#   run_scan.sh already is, means there is one definition of how the guess is
#   chosen rather than one per generated job. The scan stage did not drift
#   between the reference and the ensemble because it copies its driver; the
#   band stage did drift because its input was written from scratch.
#
# USAGE
#   Run from a frame's scan directory, after the scan has completed:
#       bash stage_tsguess.sh [output_path]
#   The default output is ../tsguess.pdb, which is where the band input expects
#   it.
#
# EXIT STATUS
#   0 on success, 1 if no maximum can be identified or the window file is
#   missing. The job that calls this should stop on a non-zero status rather
#   than leave the band to start from nothing.

set -uo pipefail
OUT="${1:-../tsguess.pdb}"

if [[ ! -s scan_energies.tsv ]]; then
  echo "TSGUESS_MISSING: no scan_energies.tsv"
  exit 1
fi

# scan_energies.tsv columns: window, achieved reaction coordinate, energy in Eh.
# The maximum energy is the least negative, so a plain numeric sort on column
# three and taking the last row gives the window at the top of the profile.
read -r W E < <(sort -k3 -g scan_energies.tsv | tail -1 | awk '{print $1, $3}')

if [[ -z "${W:-}" ]]; then
  echo "TSGUESS_MISSING: could not identify a maximum in scan_energies.tsv"
  exit 1
fi

if [[ ! -s "win_${W}.pdb" ]]; then
  echo "TSGUESS_MISSING: window $W has no win_${W}.pdb"
  exit 1
fi

cp "win_${W}.pdb" "$OUT" || { echo "TSGUESS_MISSING: copy failed"; exit 1; }

# Report the guess against TWO references, because the first window is not the
# reactant minimum and using it alone is misleading.
#
#   dE/win01        height above the first RESTRAINED window. This is what the
#                   original guard used. Window 1 is optimised under restraints
#                   from the free reactant, so on frames where the scan relaxes
#                   downhill at the start it sits well above the reactant basin
#                   and this number understates the barrier badly.
#   dE/reactant-min height above the lowest window at or before the maximum.
#                   This is the forward barrier along the scan.
#
#   Measured 2026-09-23 on frame 23268: dE/win01 = 4.84, which tripped the old
#   guard, while dE/reactant-min = 14.83, which is an ordinary barrier. Across
#   all 44 frames dE/win01 spans 4.84 to 30.31 while the corrected reference is
#   far tighter. Do NOT reference the global minimum either: on most frames it
#   lies at window 20, the product end, so that gives the REVERSE barrier.
#
# The guard now fires on dE/reactant-min.

read -r D1 DMIN WMIN < <(awk -v W="$W" -v H=627.5094740631 '
  { w[NR]=$1; e[NR]=$3; if ($1 == W) iw=NR; n=NR }
  END{
    if (iw == "") iw = n
    emin=e[1]; imin=1
    for (i=1; i<=iw; i++) if (e[i] < emin) { emin=e[i]; imin=i }
    printf "%.2f %.2f %s", (e[iw]-e[1])*H, (e[iw]-emin)*H, w[imin]
  }' scan_energies.tsv)

echo "TSGUESS window $W, $DMIN kcal/mol above the reactant-side minimum (window $WMIN)"
echo "TSGUESS_REF  $D1 kcal/mol above window 01, the old reference, kept for comparison"

LOW=$(awk -v r="$DMIN" 'BEGIN{print (r<8.0) ? 1 : 0}')
if [[ "$LOW" == "1" ]]; then
  echo "TSGUESS_LOW: the guess is less than 8 kcal/mol above the reactant-side minimum."
  echo "  The band may relax back into the reactant basin rather than climbing."
  echo "  Check the converged path for a strictly decreasing profile."
fi
exit 0
