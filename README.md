# Taxonomy Retrieval Pipeline

This Snakemake pipeline automates the process of fetching Taxonomic IDs (TaxIDs) and detailed taxonomy information for a set of genome assemblies from NCBI Datasets. The final output is a comprehensive TSV file containing taxonomy information for each assembly.

---

## Pipeline Overview

1. **List Genome Assemblies**: Identifies genome directories in the input folder.
2. **Fetch TaxIDs**: Retrieves TaxIDs for each genome using NCBI Datasets.
3. **Fetch Taxonomy**: Fetches detailed taxonomy information using the retrieved TaxIDs.
4. **Clean Output**: Processes and cleans the taxonomy data.
5. **Merge Results**: Combines all taxonomy information into a final TSV file.

---

## Usage

Run the pipeline with the following command:

```bash
snakemake --snakefile path/to/fetch_taxonomy.smk \
    --cores 8 \
    --use-conda \
    --config  HOME_DIR=/path/to/your/home \
  INPUT_GENOMES=/path/to/genome/directories \
  OUTPUT_TAXONOMY=/path/to/save/output \
  OUTPUT_TAXONOMY_NAME=desired_output_filename
```

---

## Pipeline Structure

### 1. **`checkpoint list_accessions`**
Lists all genome assembly directories in `INPUT_GENOMES` that start with `GCA_` and writes them to `assemblies.txt`.

### 2. **`fetch_taxid`**
For each assembly, fetches the TaxID and organism name from NCBI Datasets.

### 3. **`fetch_taxonomy`**
Uses the retrieved TaxID to fetch detailed taxonomy information (e.g., kingdom, phylum, class).

### 4. **`clean_taxonomy_output`**
Cleans the taxonomy output by removing JSON artifacts and formatting it as a TSV.

### 5. **`merge_results`**
Merges all cleaned taxonomy files into a final comprehensive TSV file.

### 6. **`rule all`**
Specifies the final output of the pipeline to ensure all steps are executed.

---

## Output

The final output file will be located at:

```
/output/path/{OUTPUT_TAXONOMY_NAME}_taxonomy_table.tsv
```

This TSV file will contain the following columns:

- `tax_id`
- `scientific_name`
- `rank`
- `kingdom`
- `phylum`
- `class`
- `order`
- `family`

---

## Troubleshooting

- **No directories found error:**
  - Ensure that the `INPUT_GENOMES` directory contains genome directories starting with `GCA_`.
  - Check for hidden directories like `.snakemake` that might be accidentally included.

- **Missing TaxID or Taxonomy:**
  - Some assemblies might not have corresponding TaxIDs or taxonomy information available from NCBI. The pipeline handles these by filling in `NA` values.

- **FileNotFoundError for `assemblies.txt`:**
  - This usually indicates an issue with the checkpoint. Ensure that the input directory is correctly specified and accessible.


---

## Dependencies

- **Python 3.6+**
- **Snakemake**
- **Mamba** for managing environments
- **jq** for JSON parsing in shell scripts
- NCBI **datasets** CLI tool

