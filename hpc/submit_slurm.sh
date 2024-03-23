#!/bin/bash
# SBATCH --job-name=COVID_Model_Fits_%a    # Job name
# SBATCH --mail-type=END,FAIL,ARRAY_TASKS          # Mail events (NONE, BEGIN, END, FAIL, ALL)
# SBATCH --mail-user=gregory.barnsley@lshtm.ac.uk # Where to send mail
# SBATCH --ntasks=10                    # Run on a single core
# SBATCH --array=1-7
# SBATCH --mem=80gb                     # Job memory request
# SBATCH --time=24:00:00               # Time limit hrs:min:sec
# SBATCH --output=log/COVID_Model_Fits_%a.log   # Standard output and error log
pwd; date

source ~/.bashrc

source activate R

echo $SLURM_ARRAY_TASK_ID

iso3c=$(awk -F',' -v i=1 -v j=$SLURM_ARRAY_TASK_ID 'NR==j {print $i}' "bundles.csv")
id=$(awk -F',' -v i=2 -v j=$SLURM_ARRAY_TASK_ID 'NR==j {print $i}' "bundles.csv")

if [ -f "derived/${id}" ]; then
        echo "${iso3c} already run"
else
#check if it's been unpacked
    id_nozip=$(echo ${id} | sed 's/.zip//')
    if [ -d "derived/${id_nozip}" ]; then
        rm -rf "derived/${id_nozip}"
    fi
    Rscript "script.r" $id
fi
