'''
runs_step2_all_traits.cpu.py

The purpose of this script is to run step 2 of SAIGE for all the traits we 
ran step 1 for in the CPU-based test. We are using parsl to run the scheduler.

***This needs to be run on Python/3.9 as there's an incompatability with ***
***other versions of python.                                             ***
'''

import parsl
from parsl.app.app import bash_app
from parsl.config import Config
from parsl.providers import LocalProvider
from parsl.executors import HighThroughputExecutor
from parsl.launchers import SrunLauncher
from parsl.addresses import address_by_hostname
import glob
import zmq
import os

# Directories (adjust as needed)
plink_dir = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_genotypes"
pheno_dir = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_phenotypes"
out_dir = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/CPU_test"
saige_step2_script = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/run_step2_saige_single_trait.cpu.sh"
pop = "AFR"
trait_type = "binary"

# Find all rda and varianceRatio files for the population 
rda_files = glob.glob(f"{out_dir}/{pop}/*.step1_output.rda")
ratio_files = glob.glob(f"{out_dir}/{pop}/*.step1_output.varianceRatio.txt")
#Get base file names
rda_files = [os.path.basename(x).rstrip(".step1_output.rda") for x in rda_files]
ratio_files = [os.path.basename(x).rstrip(".step1_output.varianceRatio.txt") for x in ratio_files]
#Find intersection of sets
shared_files = list(set(rda_files) & set(ratio_files))

#Get bgen and variant list file locations
bgen_files = glob.glob(f"{plink_dir}/{pop}/chr*.plink.maf.bgen")
var_files = glob.glob(f"{plink_dir}/{pop}/chr*.plink.maf.var_list.txt")
#Get chromosomes and set as bgen file names
chromosomes_bgen = [int(os.path.basename(x).rstrip(".plink.maf.bgen").lstrip("chr")) for x in bgen_files]
chromosomes_var = [int(os.path.basename(x).rstrip(".plink.maf.var_list.txt").lstrip("chr")) for x in var_files]
chromosomes= list(set(chromosomes_bgen) & set(chromosomes_var))
bgen_files = {chromosomes[x]:f"{plink_dir}/{pop}/chr{chromosomes[x]}.plink.maf.bgen" for x in range(len(chromosomes))}
var_files = {chromosomes[x]:f"{plink_dir}/{pop}/chr{chromosomes[x]}.plink.maf.var_list.txt" for x in range(len(chromosomes))}

#Cast the sample file to the burst buffer (Consider making this more robust in future)
sample_file=f"/mnt/bb/mconery/{pop}.txt"
os.system(f"sbcast {plink_dir}/{pop}/sample_list.txt {sample_file}")

#Loop the shared_files and bgen_files to make all the command options
step_2_opts=[]
for chrom in bgen_files.keys():
    for output_prefix in shared_files:
        #Check if output file exists
        if not os.path.exists(f"{out_dir}/{pop}/{output_prefix}.chr{chrom}.step2.txt"):
            #Set file locations
            bgen_file=bgen_files[chrom]
            var_file=var_files[chrom]
            model_file=f"{out_dir}/{pop}/{output_prefix}.step1_output.rda"
            variance_ratio_file=f"{out_dir}/{pop}/{output_prefix}.step1_output.varianceRatio.txt"
            #Create options/inputs to shell script
            step_2_opts.append(f"{bgen_file} {var_file} {sample_file} {model_file} {variance_ratio_file} {chrom} {output_prefix} {out_dir}/{pop}")

# Get the number of nodes:
node_raw = os.getenv("SLURM_NODELIST")
node_list = [f"frontier{x}" for x  in node_raw.lstrip("frontier[").rstrip("]").split(",")]
num_nodes = len(node_list)

# Configuration for Frontier
config = Config(
    executors=[
        HighThroughputExecutor(
            label="frontier_htex",
            cores_per_worker=1.0,  # One worker per task
            max_workers_per_node=1,  # One worker per node
            provider=LocalProvider(
                # Number of nodes job
                nodes_per_block=num_nodes,
                launcher=SrunLauncher(overrides='-c 56'),
                init_blocks=1,
                max_blocks=1,
            ),
        )
    ],
    usage_tracking=True,
)

# Load the configuration
parsl.clear()
parsl.load(config)

# Define the Parsl bash app for step 2
@bash_app
def run_saige_step2(step_2_opt, saige_step2_loc = saige_step2_script):
    return f"""
    bash {saige_step2_loc} {step_2_opt}
    """

#Print update message
print(f"NOTE: Submitting {len(step_2_opts)} tasks...")

# Submit all tasks
tasks=[]
for i, params in enumerate(step_2_opts):
    task_future = run_saige_step2(
        step_2_opt=params,
    )
    tasks.append(task_future)

print("NOTE: All tasks submitted. Waiting for completion...")

# Wait for all tasks to complete and collect results
results = []
for i, task in enumerate(tasks):
    try:
        result = task.result()
        print(f"SUCCESS: Task {i} completed successfully")
        results.append(result)
    except Exception as e:
        print(f"ERROR: Task {i} failed with error: {e}")
        results.append(None)

#Clean-up temp files
temp_files=glob.glob(f"{out_dir}/{pop}/*[0-9]")
for file in temp_files:
    os.remove(file)

print(f"SUCCESS: Workflow completed. {len([r for r in results if r is not None])} tasks succeeded.")

# Clean up
parsl.clear()

