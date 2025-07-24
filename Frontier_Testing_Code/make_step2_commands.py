#!/usr/bin/env python3

import os
import glob
import csv
from concurrent.futures import ThreadPoolExecutor

# ------------------ Adjustable Parameters -----------------------
output_path = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/step2_results/step2_commands.txt"
plink_dir = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_genotypes"
pheno_dir = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_phenotypes"
causal_var_dir = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_phenotypes/causal_variants"
inp_dir  = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/step1_results"
out_dir  = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/step2_results"
saige_step2_script = "/lustre/orion/bif154/proj-shared/arodriguez/tools/SAIGE-GPU/src/SAIGE/extdata/step2_SPAtests.R"
timing_file = "timing2.csv"
populations = ["AFR", "EUR"]
number_traits_per_job = 15
causal_var_suffix = ".effect_size_0.05.500000_space100000.txt"
minMAF=0
minMAC=0.5
LOCO="TRUE"
is_Firth_beta="TRUE"
allele_order="ref-first"
is_fastTest="FALSE"
pCutoffforFirth=0.05
num_cores = os.cpu_count()

# ------------------ Helper Functions -----------------------

def chunked_generator(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def get_chromosome_from_file(file_path):
    with open(file_path, 'r') as f:
        first_line = f.readline().strip()
        if not first_line:
            return None
        variant_id = first_line.split()[0]
        chromosome = variant_id.split(':')[0]
        return [os.path.basename(file_path).replace(causal_var_suffix, ""), chromosome]

def create_commands_for_chromosome(chromo, chromo_file_prefixes, bgen_file, var_list, pop, inp_dir, out_dir, saige_step2_script, plink_dir, num_cores, number_traits_per_job):
    cmds = []
    os.makedirs(os.path.join(out_dir, f"{pop}_manifests"), exist_ok=True)
    chunk_counter = 1
    sample_file = f"{plink_dir}/{pop}/sample_list.txt"
    for sublist in chunked_generator(chromo_file_prefixes, number_traits_per_job):
        #Make chunk directory if it doesn't exist yet
        os.makedirs(f"{out_dir}/{pop}/chunk_{chunk_counter}",  exist_ok=True)
        manifest_file = os.path.join(out_dir, f"{pop}_manifests", f"chr{chromo}.chunk_{chunk_counter}.manifest.txt")
        manifest_list = [
            [
                f"{inp_dir}/{pop}/{prefix}.step1_output.rda",
                f"{inp_dir}/{pop}/{prefix}.step1_output.varianceRatio.txt",
                f"{out_dir}/{pop}/chunk_{chunk_counter}/{prefix}.chr{chromo}.step2.txt"
            ] for prefix in sublist
        ]
        with open(manifest_file, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerows(manifest_list)
        chunk_command = (
            f"start_time=$(date +%s.%N); Rscript {saige_step2_script} --bgenFile={bgen_file} --bgenFileIndex={bgen_file}.bgi "
            f"--manifestFile={manifest_file} --sampleFile={sample_file} --chrom={chromo} --AlleleOrder=ref-first "
            f"--minMAF={minMAF} --minMAC={minMAC} --LOCO={LOCO} --is_fastTest={is_fastTest} --is_Firth_beta={is_Firth_beta} "
            f"--idstoIncludeFile={var_list} "
            f"--AlleleOrder={allele_order} --pCutoffforFirth={pCutoffforFirth} --nThreads={int(num_cores/2)}; end_time=$(date +%s.%N); "
            f'duration=$(echo "$end_time - $start_time" | bc -l); echo "chunk_{chunk_counter},chr{chromo},$duration" >> {out_dir}/{pop}/timing2.csv; '
            f'mv {out_dir}/{pop}/chunk_{chunk_counter}/*.chr{chromo}.step2.txt {out_dir}/{pop}'
        )
        #Add timing command
        cmds.append(chunk_command)
        chunk_counter += 1
    return cmds

# ------------------ Output Directories -----------------------

for pop in populations:
    os.makedirs(os.path.join(out_dir, pop), exist_ok=True)
    os.makedirs(os.path.join(out_dir, pop + '_manifests'), exist_ok=True)

cmds = []

# ------------------ Main Pipeline -----------------------

for pop in populations:
    # Collect file and chromosome info
    rda_files = glob.glob(f"{inp_dir}/{pop}/*.step1_output.rda")
    ratio_files = glob.glob(f"{inp_dir}/{pop}/*.step1_output.varianceRatio.txt")
    causal_var_files = glob.glob(f"{causal_var_dir}/Pheno*.txt")
    rda_files = [os.path.basename(x).replace(".step1_output.rda", "") for x in rda_files]
    ratio_files = [os.path.basename(x).replace(".step1_output.varianceRatio.txt", "") for x in ratio_files]
    causal_var_files = {".".join(os.path.basename(x).split(".")[:2]):x for x in causal_var_files}
    shared_files = list(set(rda_files) & set(ratio_files))
    causal_dict = {
        x: os.path.join(causal_var_dir, ".".join(x.split(".")[:2]) + causal_var_suffix)
        for x in shared_files if ".".join(x.split(".")[:2]) in causal_var_files
    }
    with ThreadPoolExecutor(max_workers=num_cores) as executor:
        causal_chromosomes = list(executor.map(get_chromosome_from_file, causal_dict.values()))
        causal_chromosomes = [x for x in causal_chromosomes if x is not None]
    chromosomes = sorted(list(set([x[1] for x in causal_chromosomes])))
    chromo_dict = {os.path.basename(".".join(x[0].split(".")[:2])): x[1] for x in causal_chromosomes}
    null_files = [x for x in shared_files if "causal" not in x]
    causal_files = [x for x in shared_files if "causal" in x]
    if null_files != []:
        chromosomes = [str(x) for x in range(1,23,1)]
    chromo_to_file_dict = {
        chromo: [x for x in causal_files if chromo_dict.get(".".join(x.split(".")[:2])) == chromo] + null_files
        for chromo in chromosomes
    }
    chromo_to_file_dict = {
        chromo: [
            prefix for prefix in chromo_to_file_dict[chromo]
            if not os.path.exists(os.path.join(out_dir, f"{prefix}.chr{chromo}.step2.txt"))
        ]
        for chromo in chromosomes
    }
    bgen_files_raw = glob.glob(f"{plink_dir}/{pop}/chr*.plink.maf.bgen")
    var_files_raw = glob.glob(f"{plink_dir}/{pop}/chr*.plink.maf.var_list.txt")
    chromosomes_bgen = [os.path.basename(x).replace(".plink.maf.bgen", "").replace("chr", "") for x in bgen_files_raw]
    chromosomes_var = [os.path.basename(x).replace(".plink.maf.var_list.txt", "").replace("chr", "") for x in var_files_raw]
    chr_intersection = sorted(list(set(chromosomes_bgen) & set(chromosomes_var) & set(chromosomes)))
    bgen_files = {chrom: f"{plink_dir}/{pop}/chr{chrom}.plink.maf.bgen" for chrom in chr_intersection}
    var_files = {chrom: f"{plink_dir}/{pop}/chr{chrom}.plink.maf.var_list.txt" for chrom in chr_intersection}
    
    # ------------ Parallel Command Generation ------------
    with ThreadPoolExecutor(max_workers=min(num_cores, len(chr_intersection))) as executor:
        future_cmds = []
        for chromo in chr_intersection:
            chromo_file_prefixes = chromo_to_file_dict.get(chromo, [])
            bgen_file = bgen_files[chromo]
            var_list  = var_files[chromo]
            if chromo_file_prefixes:
                future_cmds.append(
                    executor.submit(
                        create_commands_for_chromosome,
                        chromo, chromo_file_prefixes, bgen_file, var_list, pop, inp_dir, out_dir,
                        saige_step2_script, plink_dir, num_cores, number_traits_per_job
                    )
                )
        for future in future_cmds:
            cmds.extend(future.result())

# ------------------ Write Commands to Output -----------------------
with open(output_path, 'w') as f:
    for cmd in cmds:
        f.write(cmd + '\n')
