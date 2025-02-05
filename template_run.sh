#!/bin/bash
#SBATCH -c 12               # Number of cores for the main job
#SBATCH --mem=8G            # Memory pool for the main job
#SBATCH -J test             # Name of the main batch job

module add Mambaforge
mamba activate snakemake

snakemake --snakefile /path/to/repository/fetch_taxonomy/fetch_taxonomy.smk \
    --cores 12 \
    --use-conda \
    --config HOME_DIR=/path/to/repository/fetch_taxonomy/ \
    INPUT_GENOMES=/path/to/directory/with/genome/assemblie \
    OUTPUT_TAXONOMY=/path/where/you/want/to/see/your/output/table \
    OUTPUT_TAXONOMY_NAME=table_name  # for example: plants
