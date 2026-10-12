#!/bin/bash
# Per-set job status: finished, running, pending, failed, from lists/jobs.txt,
# one set per line:  <set name> <jobid>[,<jobid>...]
#
#   bash slurm/status.sh                 # uses lists/jobs.txt
#   bash slurm/status.sh my_jobs.txt
#
# "done" is Slurm's COMPLETED, i.e. the job exited cleanly. Whether each run
# reached 50 years is for check_block.py to confirm.
set -uo pipefail
REV="$(cd "$(dirname "$0")/.." && pwd)"
TABLE="${1:-$REV/lists/jobs.txt}"
[[ -f "$TABLE" ]] || { echo "no job table: $TABLE"; exit 1; }
count () { grep -cE "$1" <<< "$2" || true; }
printf "%-14s %6s %6s %6s %6s %6s\n" set done run pend failed other
tot=(0 0 0 0 0)
while read -r name ids; do
    [[ -z "$name" || "$name" == \#* ]] && continue
    q=$(squeue -r -h -j "$ids" -o %t 2>/dev/null)
    a=$(sacct -X -n -P -j "$ids" -o JobID,State 2>/dev/null | grep -v '\[' | cut -d'|' -f2 | awk '{print $1}')
    run=$(count '^R$' "$q"); pend=$(count '^PD$' "$q")
    done=$(count '^COMPLETED$' "$a"); fail=$(count '^(FAILED|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL)$' "$a")
    other=$(count '^(CANCELLED|PREEMPTED|BOOT_FAIL|DEADLINE)' "$a")
    printf "%-14s %6d %6d %6d %6d %6d\n" "$name" "$done" "$run" "$pend" "$fail" "$other"
    tot=($((tot[0]+done)) $((tot[1]+run)) $((tot[2]+pend)) $((tot[3]+fail)) $((tot[4]+other)))
done < "$TABLE"
printf "%-14s %6d %6d %6d %6d %6d\n" total "${tot[@]}"
(( tot[3] )) && echo "failed tasks: sacct -X -j <jobid> -o JobID%16,State,ExitCode,Elapsed | grep -vE 'COMPLETED|RUNNING|PENDING'"
exit 0
