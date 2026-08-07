import os
import pandas as pd

HOME_DIR = config["HOME_DIR"]
INPUT_GENOMES = config["INPUT_GENOMES"]
OUTPUT_TAXONOMY = config["OUTPUT_TAXONOMY"]
OUTPUT_TAXONOMY_NAME = config["OUTPUT_TAXONOMY_NAME"]
NCBI_API_KEY = config.get("ncbi_api_key", "")

os.makedirs(OUTPUT_TAXONOMY, exist_ok=True)

if not NCBI_API_KEY:
    print("Warning: No NCBI API key found! Limiting requests to 5 at a time. Set `ncbi_api_key` in config.yaml to increase the limit to 10.")

###############################################################################
# Resolve accessions ONCE at parse time (no checkpoint needed: the genome
# directories already exist on disk before the pipeline starts, so there is
# nothing "dynamic" here for Snakemake to defer until runtime).
###############################################################################
ACCESSIONS = sorted(
    d for d in os.listdir(INPUT_GENOMES)
    if os.path.isdir(os.path.join(INPUT_GENOMES, d)) and (d.startswith("GCA_") or d.startswith("GCF_"))
)

if not ACCESSIONS:
    raise ValueError(
        f"No GCA_/GCF_ directories found in {INPUT_GENOMES}. "
        "Ensure your genome directories are copied properly."
    )

print(f"Found {len(ACCESSIONS)} assemblies in {INPUT_GENOMES}")

# Also write assemblies.txt for anyone/anything downstream that expects it
# (kept for backward compatibility / manual inspection; no longer a checkpoint output)
_assemblies_path = f"{OUTPUT_TAXONOMY}/assemblies.txt"
os.makedirs(os.path.dirname(_assemblies_path), exist_ok=True)
with open(_assemblies_path, "w") as f:
    for acc in ACCESSIONS:
        f.write(acc + "\n")

wildcard_constraints:
    accession = r"GC[AF]_\d+\.\d+"

###############################################################################
# rule all
###############################################################################
rule all:
    input:
        f"{OUTPUT_TAXONOMY}/{OUTPUT_TAXONOMY_NAME}_taxonomy_table.tsv",
        f"{OUTPUT_TAXONOMY}/taxonomy_creation.log"

###############################################################################
# Rule: fetch_taxid (runs one job per accession)
###############################################################################
rule fetch_taxid:
    """
    For each discovered directory name (treated as an accession),
    fetch the TaxID and organism name from NCBI Datasets.
    """
    output:
        f"{OUTPUT_TAXONOMY}/results/{{accession}}_taxid.tsv"
    params:
        api_key = NCBI_API_KEY
    conda:
        f"{HOME_DIR}/env/ncbi-datasets.yaml"
    shell:
        """
        if [ -n "{params.api_key}" ]; then

            sleep $(awk -v min=0.5 -v max=15 'BEGIN{{srand(); print min+rand()*(max-min)}}')

            echo "Using API key for NCBI Datasets: {params.api_key}"
            datasets summary genome accession {wildcards.accession} --api-key {params.api_key} \
            | jq -r '
                if .reports then
                    (
                    .reports[]
                    | [
                        (.accession // "NA"),
                        (.organism.tax_id // "NA"),
                        (.organism.organism_name // "NA")
                        ]
                    | @tsv
                    )
                else
                    "{wildcards.accession}\tNA\tNA"
                end
            ' > {output}
        else

            sleep $(awk -v min=1 -v max=30 'BEGIN{{srand(); print min+rand()*(max-min)}}')

            datasets summary genome accession {wildcards.accession} \
            | jq -r '
                if .reports then
                    (
                    .reports[]
                    | [
                        (.accession // "NA"),
                        (.organism.tax_id // "NA"),
                        (.organism.organism_name // "NA")
                        ]
                    | @tsv
                    )
                else
                    "{wildcards.accession}\tNA\tNA"
                end
            ' > {output}
        fi
        """


###############################################################################
# Rule: fetch_taxonomy (for each accession, uses the TaxID to get full taxonomy)
###############################################################################
rule fetch_taxonomy:
    """
    Given the TaxID for each accession, fetch the detailed taxonomy info.
    Sleeps in seconds.
    """
    input:
        taxid_info = f"{OUTPUT_TAXONOMY}/results/{{accession}}_taxid.tsv"
    output:
        f"{OUTPUT_TAXONOMY}/results/{{accession}}_taxonomy.tsv"
    params:
        api_key = NCBI_API_KEY
    conda:
        f"{HOME_DIR}/env/ncbi-datasets.yaml"
    shell:
        """
        accession=$(awk -F '\\t' '{{print $1}}' {input.taxid_info})
        tax_id=$(awk -F '\\t' '{{print $2}}' {input.taxid_info})

        if [ "$tax_id" = "NA" ] || [ -z "$tax_id" ]; then

            echo -e "$accession\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA" > {output}

        elif [ -n "{params.api_key}" ]; then

            sleep $(awk -v min=1 -v max=30 'BEGIN{{srand(); print min+rand()*(max-min)}}')

            echo "Using API key for NCBI Datasets: {params.api_key}"
            datasets summary taxonomy taxon $tax_id --api-key {params.api_key}\\
            | jq -r '
                if .reports and (.reports | length > 0) then
                    [
                      "accession\\ttax_id\\tscientific_name\\trank\\tkingdom\\tphylum\\tclass\\torder\\tfamily",
                      (.reports[] | if (.taxonomy // empty) then
                          [
                              "'"$accession"'", 
                              (.taxonomy.tax_id // "NA"),
                              (.taxonomy.current_scientific_name.name // "NA"),
                              (.taxonomy.rank // "NA"),
                              (.taxonomy.classification.kingdom.name // "NA"),
                              (.taxonomy.classification.phylum.name // "NA"),
                              (.taxonomy.classification.class.name // "NA"),
                              (.taxonomy.classification.order.name // "NA"),
                              (.taxonomy.classification.family.name // "NA")
                          ]
                      else
                          ["$accession", "$tax_id", "NA", "NA", "NA", "NA", "NA", "NA"]
                      end | @tsv)
                    ]
                else
                    "$accession\\t$tax_id\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA"
                end
            ' > {output}
        else
            sleep $(awk -v min=0.5 -v max=15 'BEGIN{{srand(); print min+rand()*(max-min)}}')

            datasets summary taxonomy taxon $tax_id \\
            | jq -r '
                if .reports and (.reports | length > 0) then
                    [
                      "accession\\ttax_id\\tscientific_name\\trank\\tkingdom\\tphylum\\tclass\\torder\\tfamily",
                      (.reports[] | if (.taxonomy // empty) then
                          [
                              "'"$accession"'", 
                              (.taxonomy.tax_id // "NA"),
                              (.taxonomy.current_scientific_name.name // "NA"),
                              (.taxonomy.rank // "NA"),
                              (.taxonomy.classification.kingdom.name // "NA"),
                              (.taxonomy.classification.phylum.name // "NA"),
                              (.taxonomy.classification.class.name // "NA"),
                              (.taxonomy.classification.order.name // "NA"),
                              (.taxonomy.classification.family.name // "NA")
                          ]
                      else
                          ["$accession", "$tax_id", "NA", "NA", "NA", "NA", "NA", "NA"]
                      end | @tsv)
                    ]
                else
                    "$accession\\t$tax_id\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA"
                end
            ' > {output}
        fi
        """

###############################################################################
# Rule: clean_taxonomy_output
###############################################################################
rule clean_taxonomy_output:
    """
    Remove the JSON-like artifacts from the raw taxonomy output to produce a clean TSV.
    """
    input:
        raw_taxonomy = f"{OUTPUT_TAXONOMY}/results/{{accession}}_taxonomy.tsv"
    output:
        f"{OUTPUT_TAXONOMY}/results/{{accession}}_taxonomy_cleaned.tsv"
    shell:
        """
        sed 's/\\[//g' {input.raw_taxonomy} | \
        sed 's/\\]//g' | \
        sed 's/\\"//g' | \
        sed 's/\\\\t/\\t/g' | \
        sed 's/,$//g' | \
        sed '/^$/d' | \
        sed 's/^ *//' > {output}
        """

###############################################################################
# Rule: organize by taxonomy
# Use taxonomy file to move genome into a folder corresponding to its taxonomy
###############################################################################

rule organize_by_taxonomy:
    """
    Create symbolic links for genome files based on their taxonomy (Kingdom/Phylum/Class/Order/Family).
    """
    input:
        taxonomy_file = f"{OUTPUT_TAXONOMY}/results/{{accession}}_taxonomy_cleaned.tsv"
    output:
        temp(f"{OUTPUT_TAXONOMY}/organized/{{accession}}/linked.flag")
    run:
        import os
        import pandas as pd
        import glob

        # Validate input and find the correct genome file
        genome_files = glob.glob(f"{INPUT_GENOMES}/{wildcards.accession}/*_genomic.fna*")

        if not genome_files:
            print(f"Genome file for {wildcards.accession} is missing or not a `.fna*` file. Skipping.")
            with open(output[0], 'w') as f:
                f.write(f"Skipped {wildcards.accession} due to missing or incorrect genome file.\n")
            return

        # Load taxonomy information
        taxonomy = pd.read_csv(input.taxonomy_file, sep='\t', header=0)

        # Extract taxonomy levels
        def get_value(column_name):
            return taxonomy[column_name].values[0] if column_name in taxonomy.columns and pd.notna(taxonomy[column_name].values[0]) else f"Unknown_{column_name.capitalize()}"

        kingdom = get_value('kingdom')
        phylum = get_value('phylum')
        class_ = get_value('class')
        order = get_value('order')
        family = get_value('family')

        # Construct taxonomy-based directory
        taxonomy_path = os.path.join(OUTPUT_TAXONOMY, "organized", kingdom, phylum, class_, order, family)
        os.makedirs(taxonomy_path, exist_ok=True)

        # Ensure we only link `.fna` files
        for genome_file in genome_files:
            accession_name = os.path.basename(genome_file)
            symlink_path = os.path.join(taxonomy_path, accession_name)
            if not os.path.exists(symlink_path):
                os.symlink(os.path.abspath(genome_file), symlink_path)

        # Create a flag file to mark symlink creation
        with open(output[0], 'w') as f:
            f.write(f"Symlink created for {wildcards.accession} at {symlink_path}\n")

        print(f"Symlink created for {wildcards.accession}: {symlink_path}")

###############################################################################
# Rule: log_symlink_creation
#    This rule will collect all symlink creation messages and write them to a single log file.
###############################################################################
rule log_symlink_creation:
    """
    Aggregate symlink creation statuses into a single log file.
    """
    input:
        expand(f"{OUTPUT_TAXONOMY}/organized/{{accession}}/linked.flag", accession=ACCESSIONS)
    output:
        log_file = f"{OUTPUT_TAXONOMY}/taxonomy_creation.log"
    run:
        with open(output.log_file, 'w') as logfile:
            for temp_file in input:
                with open(temp_file, 'r') as f:
                    logfile.write(f.read())
        print(f"Symlink creation log written to {output.log_file}")


###############################################################################
# Rule: merge_results
#    Collect all {accession}_taxonomy_cleaned.tsv into one final table
###############################################################################
rule merge_results:
    """
    Merge the cleaned taxonomy files into a final comprehensive TSV table.
    """
    input:
        expand(
            f"{OUTPUT_TAXONOMY}/results/{{accession}}_taxonomy_cleaned.tsv",
            accession=ACCESSIONS
        )
    output:
        f"{OUTPUT_TAXONOMY}/{OUTPUT_TAXONOMY_NAME}_taxonomy_table.tsv"
    run:
        dfs = [pd.read_csv(file, sep='\t') for file in input]
        final_df = pd.concat(dfs, ignore_index=True).drop_duplicates()
        final_df.to_csv(output[0], sep='\t', index=False)
