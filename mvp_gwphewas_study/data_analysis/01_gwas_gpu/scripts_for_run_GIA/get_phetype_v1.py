#!/usr/bin/python
import sys, os, argparse
#import pandas 

# python ./get_phetype.py --phecode Phe_454_1 
# for i in `cat /gpfs/alpine/med112/proj-shared/data/phenotype_data/processed/PhenotypeList.Binary.ASN.HARE_ANC.filter10cases.txt`; do echo $i; if [ ! -d /ccs/home/arodriguez/med112/task0101113/output/HARE_ANC_Run/$i/ASN ]; then python /ccs/home/arodriguez/med112/task0101113/batch/pre-gwas/scripts_for_run/get_phetype_v1.py --phecode $i --ancestry ASN; fi; done

parser = argparse.ArgumentParser(description='Create Phenotype files and step1 submit files.')
parser.add_argument('--skip-phenotype', dest='skip_phenotype', action='store_true',
                    help='skip the phenotype file creation. Just create the submit files.')
parser.add_argument('--skip-submit', dest='skip_submit', action='store_true',
                    help='skip the step1 submit file creation. Just create the phenotype files.')
parser.add_argument('--phecode', dest='phecode', action='store', help='pheCode to create files for')
parser.add_argument('--ancestry', dest='ancestry', action='store', help='which ancestry to create submit files for')
args = parser.parse_args()

# goal is to create the phenotype and covariates files for a phecode for each ethnic group

# this part is always the same for each phecode
#path = "/gpfs/alpine/proj-shared/med112/task0107528/task0105060"
#path = "/gpfs/alpine/med112/proj-shared/task0110254"
#path = "/gpfs/alpine/med112/proj-shared/phenotype_data/"
path = "/gpfs/alpine/med112/proj-shared/data/phenotype_data"
PhewasCode_files = ["/gpfs/alpine/med112/proj-shared/data/phenotype_data/processed/GIA/gwphewas_binary_traits.pheno"]
quantitative_files = ["/gpfs/alpine/med112/proj-shared/data/phenotype_data/processed/GIA/gwphewas_quantitative_traits.pheno"]
#pcs_file = "/gpfs/alpine/proj-shared/med112/task0107885/Release4.mvpcoreid.pcs"
pcs_path = "/gpfs/alpine/med112/proj-shared/GIA/output/projections/smartpca_loadings_projection/ancestry_pca"
pc_files = { "EUR" : os.path.join(pcs_path, "EUR/release4.mvp.gia.EUR.pca.eigenvec"),
             "AFR" : os.path.join(pcs_path, "AFR/release4.mvp.gia.AFR.pca.eigenvec") ,
             "AMR" : os.path.join(pcs_path, "AMR/release4.mvp.gia.AMR.pca.eigenvec"),
             "EAS" : os.path.join(pcs_path, "EAS/release4.mvp.gia.EAS.pca.eigenvec")}
#hare_file = "/gpfs/alpine/proj-shared/med112/task0107885/Release4.mvpcoreid.hare"
hare_path = "/gpfs/alpine/med112/proj-shared/GIA/output/projections/smartpca_loadings_projection"
hare_files = { "EUR" : os.path.join(hare_path, "release4.mvp.gia.EUR.txt"),
               "AFR" : os.path.join(hare_path, "release4.mvp.gia.AFR.txt") ,
               "AMR" : os.path.join(hare_path, "release4.mvp.gia.AMR.txt"),
               "EAS" : os.path.join(hare_path, "release4.mvp.gia.EAS.txt") }

core_file = "/gpfs/alpine/proj-shared/med112/task0107885/MVP_Core_Opsdflt_CoreDemo_v19_2_toSummit.csv"
sex_file = "/gpfs/alpine/proj-shared/med112/task0107885/MVP_Core_Opsdflt_CoreDemo_v19_2_toSummit.GIA.csv"
#output_path = "/gpfs/alpine/proj-shared/med112//task0101113/output/HARE_ANC_Run"
output_path = "/gpfs/alpine/proj-shared/med112//task0101113/output/GIA_ANC_Run"
saige_path = "/gpfs/alpine/proj-shared/med112/task0101113/tools/saige-20210825/SAIGE/extdata"
geno_path = "/gpfs/alpine/proj-shared/med112/task0101113/output/pheCodes/inputs/genotypes/prune2"
#R_path="/ccs/home/arodriguez/med112/task0101113/tools/R-2/R-4.0.3/bin"
R_path = ""

# separate into hare groups
hare_file = hare_files[args.ancestry]
fhare = open(hare_file, "r")
#fhare.readline()   # get rid of header line
#hare = {"EAS" : [], "EUR" : [], "AFR" : [], "AMR" : []}
hare = {}
for line in fhare:
    values = line.rstrip("\n").split(" ")
    hare[values[0]] = [args.ancestry]

fhare.close()

# get the age and gender information
fcore = open(core_file, "r")
core_header = fcore.readline().rstrip("\n").split("\t")  # get rid of header line
for line in fcore:
    values = line.rstrip("\n").split("|")
    if values[0] in hare:
        #print(values[0], values[2], values[4])
        hare[values[0]].extend([values[2]])

fcore.close()

fsex = open(sex_file, "r")
for line in fsex:
    values = line.rstrip("\n").split(" ")
    if values[0] in hare:
        hare[values[0]].extend([values[1]])

fsex.close()

# get the pcs into the hare groups
pcs_file = pc_files[args.ancestry]
fpcs = open(pcs_file, "r")
pc_header = fpcs.readline().rstrip("\n").split("\t")[2:]   # get rid of header line
#pcs = {}
for line in fpcs:
    values = line.rstrip("\n").split("\t")
    hare[values[0]].extend(values[2:])

fpcs.close()

# get the gender specific phenotypes
female_binary_file = "/gpfs/alpine/med112/proj-shared/data/phenotype_data/processed/PhenotypeList.Binary.Female.Traits.v2.txt"
male_binary_file = "/gpfs/alpine/med112/proj-shared/data/phenotype_data/processed/PhenotypeList.Binary.Male.Traits.v2.txt"
female_quant_file = "/gpfs/alpine/med112/proj-shared/data/phenotype_data/processed/PhenotypeList.Quantitative.Female.Traits.txt"

gender_specific = {}
ffb_h = open(female_binary_file, "r")
for i in ffb_h:
    gender_specific[i.rstrip("\n")] = "F"
ffb_h.close()

fmb_h = open(male_binary_file, "r")
for i in fmb_h:
    gender_specific[i.rstrip("\n")] = "M"
fmb_h.close()

ffq_h = open(female_quant_file, "r")
for i in ffq_h:
    gender_specific[i.rstrip("\n")] = "F"
ffq_h.close()

### this is unique to each phecode

phecode = args.phecode
#out_file = sys.argv[2]
index = None
phecodeFlag = None
for i in PhewasCode_files:
    #fh = open("%s/%s" % (path, i), "r")
    fh = open(i, "r")
    header = fh.readline().replace("\"", "").rstrip("\n").rstrip("\r").split("\t")
    if phecode in header:
        # get the index
        index = header.index(phecode)
        print(header[index])
        # get fields 1, 2 and index columns
        #print index
        phecode_values = []
        for line in fh:
            values = line.rstrip("\n").rstrip("\r").replace("\"", "").split("\t")
            phecode_values.append([values[0], values[index]])
        # exit for loop
        fh.close()
        phecodeFlag = "binary"
        break
    fh.close()

if index is None:
    for i in quantitative_files:
        fh = open(i, "r")
        #header = fh.readline().replace("\"", "").rstrip("\n").split("\t")
        header = fh.readline().replace("\"", "").rstrip("\n").rstrip("\r").split("\t")
        if phecode in header:
            # get the index
            index = header.index(phecode)
            print(header[index])
            # get fields 1, 2 and index columns
            #print index
            phecode_values = []
            for line in fh:
                #values = line.rstrip("\n").replace("\"", "").split("\t")
                values = line.rstrip("\n").rstrip("\r").replace("\"", "").split("\t")
                phecode_values.append([values[0], values[index]])
            # exit for loop
            fh.close()
            phecodeFlag = "quantitative"
            break
        fh.close()

print(phecodeFlag)

if index is None:
    print("!!!ERROR: The phenotype %s cannot be found." % phecode)
    sys.exit()

#fo = open(out_file, "w")
#pheFiles = {"EAS" : [], "EUR" : [], "AFR" : [], "AMR" : []}
pheFiles = {args.ancestry : []}

for row in phecode_values:
    pheValue = row[1]
    subject_type = phecodeFlag

    if row[0] in hare:
        if hare[row[0]][0] in pheFiles:
            a = [row[0], pheValue]
            a.extend(hare[row[0]][1:])
            pheFiles[hare[row[0]][0]].append(a)

pheHeader = ["MVPCore_ID", "Phenotype", "Age", "Gender"]
pheHeader.extend(pc_header)

for group in pheFiles:
    phecode_path = "%s/%s/%s/inputs" % (output_path, phecode, group)
    out_file = "%s/PhenoFile.%s.%s.txt" % (phecode_path, phecode, group)
    if (not os.path.exists(phecode_path)):
        os.makedirs(phecode_path)
    fo = open(out_file, "w")
    fo.write("%s\n" % "\t".join(pheHeader))
    for line in pheFiles[group]:
        fo.write("%s\n" % ("\t".join(line)))
    fo.close()
    print("Printing pheno for %s - %s at: %s" % (phecode, group, phecode_path))

## create submit files
for group in pheFiles:
    phecode_path = "%s/%s/%s/inputs" % (output_path, phecode, group)
    submit_path = "%s/%s/%s/submit/step1" % (output_path, phecode, group)
    submit_file = "%s/submit.step1.%s.%s.txt" % (submit_path, phecode, group)
    submit2_path = "%s/%s/%s/submit/step2" % (output_path, phecode, group)
    if (not os.path.exists(submit_path)):
        os.makedirs(submit_path)
    if (not os.path.exists(submit2_path)):
        os.makedirs(submit2_path)

    if group == "EAS" or group == "AMR" or group == "AFR":
        walltime = "2:00"
        batch = "batch"
        gpu_qty = 6
        gpu_conf = "~/gemm/nvblas.conf"
    elif group == "EUR":
        walltime = "24:00"
        batch = "batch-hm"
        gpu_qty = 6
        gpu_conf = "~/gemm/nvblas.conf"

    fs = open(submit_file, "w")
    fs.write("#!/bin/bash\n#BSUB -nnodes 1\n#BSUB -W %s\n#BSUB -P MED112\n#BSUB -q %s\n" % (walltime, batch))
    fs.write("#BSUB -o %s.%s.stdout\n#BSUB -e %s.%s.stderr\n#BSUB -J %s.%s\n" % (phecode, group, phecode, group, phecode, group))
    fs.write("#BSUB -alloc_flags \"nvme\"\n\n")
    fs.write("phenoID=\"%s\"\ngroup=\"%s\"\n" % (phecode, group))
    fs.write("saige_path=\"%s\"\n" % (saige_path))
    fs.write("geno_path=\"%s\"\n" % (geno_path))
    step1_outpath = "%s/%s/%s/step1" % (output_path, phecode, group)
    step2_outpath = "%s/%s/%s/step2" % (output_path, phecode, group)
    if (not os.path.exists(step1_outpath)):
        os.makedirs(step1_outpath)
    if (not os.path.exists(step2_outpath)):
        os.makedirs(step2_outpath)

    fs.write("output_path=\"%s\"\n" % (step1_outpath))
    nthreads = 42
    nvme_path = "/mnt/bb/arodriguez"
    fs.write("nvme_path=%s\n\n" % (nvme_path))
    fs.write("module load gcc/10.2.0\n")
    fs.write("module load cuda/11.1.1\n")
    fs.write("module load openblas\n")
    fs.write("module load bzip2\n")
    fs.write("module load hdf5\n")
    geno_name = "20200917.GenotypeData.Release4.mvpcoreid.maf.hw_dropped.ld"
    fs.write("jsrun -n 1 cp -r %s/%s.fam %s/.\n" % (geno_path, geno_name, nvme_path))
    fs.write("jsrun -n 1 cp -r %s/%s.bim %s/.\n" % (geno_path, geno_name, nvme_path))
    fs.write("jsrun -n 1 cp -r %s/%s.bed %s/.\n" % (geno_path, geno_name, nvme_path))
    phenoFile = "%s/PhenoFile.%s.%s.txt" % (phecode_path, phecode, group)
    fs.write("jsrun -n 1 cp -r %s %s/.\n" % (phenoFile, nvme_path))
    out_prefix = "phewas.ld.maf.%s.%s.out" % (phecode, group)
    out_grm = "%s.grm" % (out_prefix)
    log_path = "%s/phewas.ld.maf.%s.%s.log" % (step1_outpath, phecode, group)

    gender_specific_params = ""
    gender_code = "NA"
    if phecode in gender_specific:
        gender_type = gender_specific[phecode]
        if gender_type == "F":
            gender_specific_params = " --FemaleOnly TRUE --FemaleCode 2 --sexCol Gender "
            gender_code = "F"
        elif gender_type == "M":
            gender_specific_params = " --MaleOnly TRUE --MaleCode 1 --sexCol Gender "
            gender_code = "M"

    fs.write("jsrun -n1 -c42 -g%s --smpiargs=\"-disable_gpu_hooks\" -bpacked:42 -EOMP_NUM_THREADS=42 -E NVBLAS_CONFIG_FILE=%s -E LD_PRELOAD=$OLCF_CUDA_ROOT/lib64/libnvblas.so Rscript %s/step1_fitNULLGLMM.R --plinkFile=%s/%s --phenoFile=%s/%s --outputPrefix=%s/%s --phenoCol=Phenotype --covarColList=Age,Gender,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --sampleIDColinphenoFile=MVPCore_ID --traitType=%s --nThreads=%s --minMAFforGRM 0.05 --maxiterPCG 500 --maxiter 25 --LOCO FALSE %s > %s\n" % (gpu_qty, gpu_conf, saige_path, nvme_path, geno_name, nvme_path, os.path.basename(phenoFile), nvme_path, out_prefix, phecodeFlag, nthreads, gender_specific_params, log_path))
    fs.write("jsrun -n 1 cp -r %s/%s.rda %s/.\n" % (nvme_path, out_prefix, step1_outpath))
    fs.write("jsrun -n 1 cp -r %s/%s_30markers.SAIGE.results.txt %s/.\n" % (nvme_path, out_prefix, step1_outpath))
    fs.write("jsrun -n 1 cp -r %s/%s.varianceRatio.txt %s/.\n" % (nvme_path, out_prefix, step1_outpath))
    #fs.write("jsrun -n 1 cp -r %s/%s %s/.\n" % (nvme_path, out_grm, step1_outpath))
    fs.close()
    print("Submit files at: %s" % (submit_file))

    # create RDS file for step1
    #scripts_path = "/gpfs/alpine/proj-shared/med112/task0101113/batch/pre-gwas/scripts_for_run"
    scripts_path = "/gpfs/alpine/proj-shared/med112/task0101113/batch/pre-gwas/gwPheWAS-Summit/scripts_for_run_GIA/"
    cmd = "Rscript %s/create_step1_rds.r %s %s %s/step1.%s.%s.rds %s %s" % (scripts_path, phecode, group, step1_outpath, phecode, group, phecodeFlag, gender_code)
    os.system(cmd)

    submit1_path = "%s/%s/%s/submit/step1" % (output_path, phecode, group)
    submit1_file = "%s/submit.step1.rds.%s.%s.txt" % (submit2_path, phecode, group)
    s1_nodes = 13
    s1_cpus = s1_nodes * 6
    s1_walltime = "2:00"

    fs1 = open(submit1_file, "w")
    fs1.write("#!/bin/bash\n")
    fs1.write("#BSUB -nnodes %s\n" % str(s1_nodes))
    fs1.write("#BSUB -W %s\n" % s1_walltime)
    fs1.write("#BSUB -q batch\n")
    fs1.write("#BSUB -P MED112\n")
    fs1.write("#BSUB -o %s.s1.%s.rds.stdout\n" % (phecode, group))
    fs1.write("#BSUB -e %s.s1.%s.rds.stderr\n" % (phecode, group))
    fs1.write("#BSUB -J %s.s1.%s\n\n" % (phecode, group))
    fs1.write("module load gcc/10.2.0\n")
    fs1.write("module load cuda/11.1.1\n")
    fs1.write("module load openblas\n")
    fs1.write("module load bzip2\n")
    fs1.write("module load hdf5\n")

    fs1.write("jsrun -n%s -a1 -g1 -r6 %s/Rscript %s/submit_step1_with_nvme.r %s/step1.%s.%s.rds %s %s\n" % (str(s1_cpus), R_path, scripts_path,  step1_outpath, phecode, group, phecode, group))
    fs1.close()

    # create RDS file for step2
    cmd = "Rscript %s/create_step2_rds.r %s %s %s/step2.%s.%s.rds" % (scripts_path, phecode, group, step2_outpath, phecode, group)
    os.system(cmd)

    submit2_path = "%s/%s/%s/submit/step2" % (output_path, phecode, group)
    submit2_file = "%s/submit.step2.%s.%s.txt" % (submit2_path, phecode, group)
    s2_nodes = 92
    s2_cpus = s2_nodes * 42
    s2_walltime = "12:00"

    fs2 = open(submit2_file, "w")
    fs2.write("#!/bin/bash\n")
    fs2.write("#BSUB -nnodes %s\n" % str(s2_nodes))
    fs2.write("#BSUB -W %s\n" % s2_walltime)
    fs2.write("#BSUB -q batch\n")
    fs2.write("#BSUB -P MED112\n")
    fs2.write("#BSUB -o %s.s2.%s.stdout\n" % (phecode, group))
    fs2.write("#BSUB -e %s.s2.%s.stderr\n" % (phecode, group))
    fs2.write("#BSUB -J %s.s2.%s\n\n" % (phecode, group))
    fs2.write("module load gcc/10.2.0\n")
    fs2.write("module load cuda/11.1.1\n")
    fs2.write("module load openblas\n")
    fs2.write("module load bzip2\n")
    fs2.write("module load hdf5\n")

    fs2.write("jsrun -n%s -a1 -c1 -r42 Rscript %s/submit_step2_with_nvme.r %s/step2.%s.%s.rds %s %s\n" % (str(s2_cpus), scripts_path,  step2_outpath, phecode, group, phecode, group))
    fs2.close()
