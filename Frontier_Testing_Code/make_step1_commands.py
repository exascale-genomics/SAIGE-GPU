#!/usr/bin/env python3

import os
import glob

# Set directory paths and variables
plink_dir = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_genotypes"
pheno_dir = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_phenotypes"
out_dir = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/step1_results"
path_to_saige="/lustre/orion/bif154/proj-shared/arodriguez/tools/SAIGE-GPU/src/SAIGE"
saige_step1_script = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/run_step1_saige_single_trait.sh"
output_path = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/step1_results/step1_commands.txt"

populations = ["AFR", "EUR"]
trait_types = ["binary", "continuous"]

cmds = []

causal_flag = False #Flag to turn off making commands for causal phenotypes

for pop in populations:
    plink_file_post = f"/mnt/bb/mconery/{pop}"
    pop_out_dir = os.path.join(out_dir, pop)

    for trait_type in trait_types:
        pattern = os.path.join(pheno_dir, pop, f"Pheno*.{trait_type}*txt")
        phenotype_files = glob.glob(pattern)
        if causal_flag == False:
            phenotype_files = [x for x in phenotype_files if 'causal' not in x]

        for file in phenotype_files:
            temp = os.path.basename(file)
            prefix = temp[:-4] if temp.endswith('.txt') else temp
            step1_output_rda = os.path.join(pop_out_dir, f"{prefix}.step1_output.rda")
            step1_output_variance = os.path.join(pop_out_dir, f"{prefix}.step1_output.varianceRatio.txt")

            if not (os.path.exists(step1_output_rda) and os.path.exists(step1_output_variance)):
                cmd = f"{saige_step1_script} {plink_file_post} {file} {prefix} {pop_out_dir} {path_to_saige} {trait_type}"
                cmds.append(cmd)

with open(output_path, 'w') as f:
    for cmd in cmds:
        f.write(cmd + '\n')
