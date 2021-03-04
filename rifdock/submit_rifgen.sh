export PYTHON_PATH=/wynton/home/kortemme/krivacic/software/anaconda3/bin/python
export RUNDIR=/wynton/home/kortemme/krivacic/intelligent_design/helix_matcher/test_files/boundary/
qsub -h -cwd -o logs/ -e errs/ -t 107 -l h_rt=6:00:00 -l mem_free=6G -b y \
-N rifgen $PYTHON_PATH \
run_rifgen.py $RUNDIR --sge
