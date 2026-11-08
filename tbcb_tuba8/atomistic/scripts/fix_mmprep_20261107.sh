#!/bin/bash
# The MM/GBSA topologies (mmprep/) were written on 2026-11-04, before four builds were rebuilt all-L on 2026-11-05
# (c_tbce_a8_m0, h_9M1M_Tuba8, x_9M1L_tbce_a1a, x_9M1L_tbce_a8): x_9M1L_tbce_a1a no longer matched its trajectories
# (15147 vs 15144 atoms) and the other three matched only in atom count.  Rewrite those four from the current builds
# (same receptor chains as before; the old directories go to superseded_20261107/mmprep), then run the mmgbsa stage
# for all ten systems with a12b, which now also checks residue/atom names atom by atom.
# Usage: sbatch -p gpu -G 1 -c 4 --mem=24G -t 1:30:00 -J ana_mmgbsa -o <log> --wrap "bash -l scripts/fix_mmprep_20261107.sh"
module load chemistry
module load py-openmm/8.1.1_py312
export PYTHONPATH=$PYTHONPATH:$HOME/fep_pylib
cd $HOME/tbcb_atomistic
M=$SCRATCH/tbcb_atomistic/mmprep; B=$SCRATCH/tbcb_atomistic/build; SUP=$SCRATCH/tbcb_atomistic/superseded_20261107/mmprep
mkdir -p $SUP
for s in c_tbce_a8_m0 h_9M1M_Tuba8 x_9M1L_tbce_a1a x_9M1L_tbce_a8; do
  src=$M/$s; [ -d $SUP/$s ] && src=$SUP/$s
  rec=$(python3 -c "import json; print(','.join(json.load(open('$src/prep.json'))['receptor_chains']))")
  [ -d $SUP/$s ] || mv $M/$s $SUP/$s
  python3 -W ignore scripts/a12_mmpbsa.py prep $B/$s $M/$s --receptor $rec > $SCRATCH/tbcb_atomistic/logs/mmprep_$s.log 2>&1
  echo "$s receptor $rec prep exit $?: $(tail -1 $SCRATCH/tbcb_atomistic/logs/mmprep_$s.log | cut -c1-200)"
done
bash scripts/pipeline_stage.sh mmgbsa
