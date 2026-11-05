#!/bin/bash
# Rebuild the builds that contained D residues (2026-11-05), with the fixed builders: a01/a07 rebuild any side chain
# that PDBFixer put on the D side (enforce_L), a07 uses a proper rotation, a02 holds every N-CA-C-CB dihedral near
# its L value in both minimisation stages and refuses to write a structure in which a residue changed handedness.
# Output goes to build_L/<sid>; nothing in build/ is touched (the swap is a separate, checked step).
# Marlowe CPU array: sbatch -p preempt -A marlowe-m000251 -c 8 --mem=16G -t 1:30:00 --exclude=n18,n28 \
#                      --array=0-10 -o logs/rebuildL_%A_%a.out scripts/rebuild_L_20261105.sh
# then (after the array): scripts/rebuild_L_20261105.sh finalize   (a8_dom from the new a8, strict chirality check)
cd /scratch/m000251/zongminl/tbcb_mech/atomistic
PY="/scratch/m000251/zongminl/mm/envs/amb/bin/python -W ignore"
export OPENMM_CPU_THREADS=${SLURM_CPUS_PER_TASK:-8} OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
B=build_L; mkdir -p $B
AF=inputs/af3
if [ "$1" = "finalize" ]; then
  $PY scripts/a16_resolvate.py $B/a8 $B/a8_dom --pad 2.5 && echo "a8_dom done"
  $PY scripts/a29_chirality_check.py --strict $B/*/minimized.cif.gz > $B/a29_strict.txt 2>&1
  echo "a29 strict exit $?" >> $B/a29_strict.txt
  cat $B/a29_strict.txt
  exit 0
fi
case ${SLURM_ARRAY_TASK_ID} in
  0) $PY scripts/a01_build.py a8 'Tuba8' --out $B && $PY scripts/a02_minimize.py $B/a8 ;;
  1) $PY scripts/a01_build.py a8_S114I 'Tuba8:S114I' --out $B && $PY scripts/a02_minimize.py $B/a8_S114I ;;
  2) $PY scripts/a01_build.py a1a_loop8 'Tuba1a:loop=Tuba8' --out $B && $PY scripts/a02_minimize.py $B/a1a_loop8 ;;
  3) $PY scripts/a01_build.py a8_A40K 'Tuba8:A40K' --out $B && $PY scripts/a02_minimize.py $B/a8_A40K ;;
  4) $PY scripts/a01_build.py a8_I43G 'Tuba8:I43G' --out $B && $PY scripts/a02_minimize.py $B/a8_I43G ;;
  5) $PY scripts/a01_build.py a1a_loop8_I114S 'Tuba1a:loop=Tuba8,I114S' --out $B && $PY scripts/a02_minimize.py $B/a1a_loop8_I114S ;;
  6) $PY scripts/a07_build_complex.py h_9M1M_Tuba8 exp inputs/9M1M.cif.gz --alpha a --partners b --keep-beta b --thread Tuba8 --out $B \
       && $PY scripts/a02_minimize.py $B/h_9M1M_Tuba8 --chains-a a --chains-b b ;;
  7) $PY scripts/a07_build_complex.py x_9M1L_tbce_a1a exp inputs/9M1L.cif.gz --alpha a --partners E --thread Tuba1a --out $B \
       && $PY scripts/a02_minimize.py $B/x_9M1L_tbce_a1a --chains-a E --chains-b a ;;
  8) $PY scripts/a07_build_complex.py x_9M1L_tbce_a8 exp inputs/9M1L.cif.gz --alpha a --partners E --thread Tuba8 --out $B \
       && $PY scripts/a02_minimize.py $B/x_9M1L_tbce_a8 --chains-a E --chains-b a ;;
  9) $PY scripts/a07_build_complex.py c_tbce_a8_m0 af3 $AF/fold_tbce_tuba8/fold_tbce_tuba8_model_0.cif --alpha A --partners B --out $B \
       && $PY scripts/a02_minimize.py $B/c_tbce_a8_m0 --chains-a B --chains-b A ;;
  10) mkdir -p $B/t_a1aYF_K188A_cplx && cp -p build/t_a1aYF_K188A_cplx/solvated.cif.gz build/t_a1aYF_K188A_cplx/build.json $B/t_a1aYF_K188A_cplx/ \
       && $PY scripts/a02_minimize.py $B/t_a1aYF_K188A_cplx --chains-a A --chains-b B ;;
esac
rc=$?
echo "task ${SLURM_ARRAY_TASK_ID} exit $rc"
exit $rc
