export PYTHON_PATH=/wynton/home/kortemme/krivacic/software/anaconda3/bin/python
export RUNDIR=$1
export num_tasks=$(echo $(ls $RUNDIR | grep patch | wc -l))
qsub -cwd -o logs/ -e logs/ -t 1-$num_tasks -l h_rt=6:00:00 -l mem_free=6G -l scratch=200M -b y \
-N rifdock_full $PYTHON_PATH \
run_sge.py $RUNDIR --sge
