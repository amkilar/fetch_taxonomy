#!/bin/bash
#SBATCH -c 12                # Number of cores for the main job
#SBATCH --mem=32G            # Memory pool for the main job
#SBATCH -t 07-00:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -J taxonomy          # Name of the main batch job
#SBATCH -p eddy              # Partition to submit to
#SBATCH -o /n/eddy_lab/users/akilar/fetch_taxonomy/%A_%a.out
#SBATCH -e /n/eddy_lab/users/akilar/fetch_taxonomy/%A_%a.err


mkdir -p /n/eddy_lab/users/akilar/fetch_taxonomy/Verterbrates

#module add Mambaforge
module add Miniforge3/24.7.1-fasrc01
mamba activate /n/home10/akilar/software/env_snakemake

snakemake --snakefile /n/eddy_lab/users/akilar/fetch_taxonomy/fetch_taxonomy.smk \
    --cores 1 \
    --use-conda \
    --config HOME_DIR=/n/eddy_lab/users/akilar/fetch_taxonomy \
    INPUT_GENOMES=/n/eddy_lab/data/RNAhub_genomes/Vertebrate_reference_genomes/genomes \
    OUTPUT_TAXONOMY=/n/eddy_lab/users/akilar/fetch_taxonomy/Verterbrates \
    OUTPUT_TAXONOMY_NAME=Verterbrates  \
    --unlock
    #ncbi_api_key=55d74cdd9fb0d170b8cdb6fa59056bd57309 \


    snakemake --snakefile /n/eddy_lab/users/akilar/fetch_taxonomy/fetch_taxonomy.smk \
    --cores 12 \
    --use-conda \
    --config HOME_DIR=/n/eddy_lab/users/akilar/fetch_taxonomy \
    INPUT_GENOMES=/n/eddy_lab/data/RNAhub_genomes/Vertebrate_reference_genomes/genomes \
    OUTPUT_TAXONOMY=/n/eddy_lab/users/akilar/fetch_taxonomy/Verterbrates \
    OUTPUT_TAXONOMY_NAME=Verterbrates  \
    ncbi_api_key=55d74cdd9fb0d170b8cdb6fa59056bd57309 \

