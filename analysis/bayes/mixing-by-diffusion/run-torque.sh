# FOR SUBMITING JOB TO HAL COMPUTER CLUSTER
#
# walltime : maximum wall clock time (hh:mm:ss)
#PBS -l walltime=48:00:00
#
# join stdout and stderr # commented out
# PBS -j oe 
#
# spool output immediately # commented out
# PBS -k oe
#
# specify directory for stdout and stderr files
#PBS -o ./
#PBS -e ./
#
# specify queue
#PBS -q batch
#
#   nodes: number of nodes
#   ppn: how many cores per node to use
# This is to ask for 1 CPU only
#PBS -l nodes=1:ppn=1
#
# Specify memory
#PBS -l mem=4gb
#
# export all my environment variables to the job
##PBS -V
#
# job name (default = name of script file)
#PBS -N diffusion_r1

if [ -n "$PBS_O_WORKDIR" ]; then
    cd $PBS_O_WORKDIR
fi

# Activate the environment where Assaytools is installed
source activate py35

# Run my program
echo "Started running Assaytools..."
quickmodel --inputs 'inputs_p38_singlet_diffusion' --type 'singlet' --nsamples 100000
echo "Done."
