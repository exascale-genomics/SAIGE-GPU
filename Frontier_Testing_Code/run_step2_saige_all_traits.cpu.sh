#!/bin/bash
#SBATCH -A bif154
#SBATCH -J saige_step2_cpu
#SBATCH -o run_step2_saige_all_traits.cpu.log
#SBATCH -N 15
#SBATCH -q debug
#SBATCH -t 30:00
#SBATCH -C nvme

##################################################################################################################
############################################## Load Modules and Conda Env ########################################
##################################################################################################################

##### START OF SBCAST AND CONDA-UNPACK #####
#Set variables
ENV_NAME=parsl_bash
TAR_FILE=${ENV_NAME}.tar.gz
TAR_DIR=/lustre/orion/bif154/proj-shared/mconery/tools
NNODES=${SLURM_NNODES}
# Move a copy of the env to the NVMe on each node
echo "copying ${ENV_NAME} to each node in the job"
sbcast -pf ${TAR_DIR}/${TAR_FILE} /mnt/bb/${USER}/${TAR_FILE}
if [ ! "$?" == "0" ]; then
    # CHECK EXIT CODE. When SBCAST fails, it may leave partial files on the compute nodes, and if you continue to launch srun,
    # your application may pick up partially complete shared library files, which would give you confusing errors.
    echo "SBCAST failed!"
    exit 1
fi

# Untar the environment file (only need 1 task per node to do this)
srun -N ${NNODES} --ntasks-per-node 1 mkdir -p /mnt/bb/${USER}/${ENV_NAME}
echo "untaring ${ENV_NAME}"
srun -N ${NNODES} --ntasks-per-node 1 tar -xzf /mnt/bb/${USER}/${TAR_FILE} -C /mnt/bb/${USER}/${ENV_NAME}

# Unpack the env
source activate /mnt/bb/${USER}/${ENV_NAME}
srun -N ${NNODES} --ntasks-per-node 1 conda-unpack
##### END OF SBCAST AND CONDA-UNPACK #####

##################################################################################################################
################################################## Execute SAIGE Step 2 ##########################################
##################################################################################################################

#Run the script
saige_step2_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/run_step2_saige_all_traits.cpu.py"
python $saige_step2_script
