#!/usr/bin/env bash
# Drug-response shape (user, 2026-10-07: "feel free to use Hill and linear"): the main run carries the Stage 2b
# drug model the selection rule picks; this variant forces the 3-compartment Hill kill (chain3_hill, with the
# same Norton–Simon coupling and weights) and reruns the control on it, so the reports can show both shapes.
#   bash sensitivity/run_hill_variant.sh          (after a full pipeline run; ~1 h on 20 threads)
# Point fits only (no bootstrap). Scenario B's crowding comes from the main run's untreated-coculture fit
# (outputs/zoo_analyses/analyses_state.jls), the only part not refitted here.
cd "$(dirname "$0")/.." || exit 1
D=outputs/sensitivity/hill; mkdir -p "$D"
python - <<'PYEOF'
import re
s = open('stages.toml', encoding='utf-8').read()
s2 = re.sub(r'(?m)^(id = "s2"\r?\n)', r'\1only = ["chain3_hill"]   # variant: Hill kill forced in Stage 2b\n', s, count=1)
assert s2 != s
open('outputs/sensitivity/hill/stages.toml', 'w', encoding='utf-8', newline='\n').write(s2)
PYEOF
LOG=$D/progress.log
echo "$(date '+%H:%M') hill fits" >> $LOG
STAGES_TOML=$D/stages.toml ZOO_OUT=$D BOOT_B=0 julia -t 16 --project=. run_zoo_fits.jl > $D/fits.log 2>&1 || { echo "fits failed" >> $LOG; exit 1; }
for s in run_control.jl run_control_fasting.jl; do
  echo "$(date '+%H:%M') hill $s" >> $LOG
  ZOO_STATE=$D/zoo_state.jls AN_STATE=outputs/zoo_analyses/analyses_state.jls CONTROL_OUT=$D REDO=1 julia -t 16 --project=. $s > $D/$s.log 2>&1 || echo "$s failed" >> $LOG
done
echo "$(date '+%H:%M') hill done" >> $LOG
