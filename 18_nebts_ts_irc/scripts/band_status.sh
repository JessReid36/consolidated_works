#!/usr/bin/env bash
# Which live band jobs are ready to cancel before the transition-state stage.
# Ready = band converged AND a PATH SUMMARY written AND the climbing-image
# geometry on disk. Those three are all that anything downstream needs.
ENS="$HOME/system_development/05_qmmm/19_ensemble"
Q=$(qstat -u "$USER" 2>/dev/null)
kill_list=""
printf '%-7s %-16s %-4s %-7s %-12s %-7s %-8s %s\n' \
       frame jobid st 'neb.out' updated barrier hessian verdict
printf -- '---------------------------------------------------------------------------------------\n'
for d in "$ENS"/frame_*; do
  f=$(basename "$d"); f=${f#frame_}
  read -r jid state < <(printf '%s\n' "$Q" | awk -v p="cm19_n${f}" '
    $1 ~ /^[0-9]+\./ { n=$4; sub(/\*$/,"",n)
                       if (n ~ /^cm19_n/ && index(p,n)==1) { print $1, $(NF-1); exit } }')
  [ -n "$jid" ] || continue
  out="$d/neb.out"; ci="$d/neb_NEB-CI_converged.QMRegion.xyz"; hes="$d/neb.appr.hess"
  if [ -f "$out" ]; then
    lines=$(wc -l < "$out"); upd=$(date -r "$out" +'%m-%d %H:%M')
    conv=$(grep -c "THE NEB OPTIMIZATION HAS CONVERGED" "$out")
    psum=$(grep -c "PATH SUMMARY" "$out")
  else
    lines=0; upd="--"; conv=0; psum=0
  fi
  if [ -s "$ci" ]; then cig="yes"; else cig="no"; fi
  if [ -s "$hes" ]; then hs=$(du -h "$hes" 2>/dev/null | cut -f1); else hs="-"; fi
  bar="-"
  if [ "$psum" -gt 0 ]; then
    bar=$(awk '/PATH SUMMARY/{p=1} p && /<= CI/{print $4; exit}' "$out")
    [ -n "$bar" ] || bar="-"
  fi
  if [ "$conv" -gt 0 ] && [ "$psum" -gt 0 ] && [ "$cig" = yes ] && [ "$hs" != "-" ]; then
    verdict="READY - qdel (saddle already started, barrier secured)"; kill_list="$kill_list $jid"
  elif [ "$conv" -gt 0 ] && [ "$psum" -gt 0 ] && [ "$cig" = yes ]; then
    verdict="READY - qdel"; kill_list="$kill_list $jid"
  elif [ "$conv" -gt 0 ]; then verdict="converged, CI file not written yet"
  elif [ "$state" = Q ];    then verdict="queued, not started"
  else                            verdict="band running"
  fi
  printf '%-7s %-16s %-4s %-7s %-12s %-7s %-8s %s\n' \
         "$f" "${jid%%.*}" "$state" "$lines" "$upd" "$bar" "$hs" "$verdict"
done
echo
if [ -n "$kill_list" ]; then
  echo "Ready to cancel before the saddle stage:"
  echo "    qdel$kill_list"
  echo "  (barrier and CI geometry are already on disk for these)"
else
  echo "Nothing ready to cancel yet."
fi
