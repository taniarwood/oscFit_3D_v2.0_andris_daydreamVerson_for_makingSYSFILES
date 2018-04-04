#!/bin/zsh
#(otherwise the default shell would be used)
#$ -S /bin/zsh
#
#(the cpu time for this job, less than 30min to start quickly)
#$ -l h_cpu=03:29:00 
#
#(the maximum memory usage of this job)
#$ -l h_rss=2G
#
#(stderr and stdout are merged together to stdout)
#$ -j y
#
###(send mail on job's begin, end and abort bea) bea
### -m ea
#
# execute job from current directory and not relative to your home directory
#$ -cwd
#
# send output files to the trash, COMMENT TO SEE THE OUTPUT AND DEBUG
#$ -o /dev/null
###$ -o /lustre/fs15/group/icecube/yanezjua/analysis/fitter_logs
#
#$ -P z_nuastr
###$ -P icecube
#

SCRIPT=$1
NR_JOBS=$2
INDIR=$3
TASK=`expr $SGE_TASK_ID - 1`
python $SCRIPT $TASK $NR_JOBS $INDIR
