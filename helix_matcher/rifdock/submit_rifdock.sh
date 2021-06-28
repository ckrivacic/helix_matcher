export PYTHON_PATH=/wynton/home/kortemme/krivacic/software/anaconda3/bin/python
export RUNDIR=$1
export num_tasks=$(echo $(ls $RUNDIR | grep patch | wc -l))
qsub -cwd -o logs/ -e logs/ -t $num_tasks -l h_rt=6:00:00 -l mem_free=6G -b y \
-N rifgen $PYTHON_PATH \
run_rifdock.py $RUNDIR --sge
