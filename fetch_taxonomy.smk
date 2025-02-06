import os
import pandas as pd

HOME_DIR = config["HOME_DIR"]
INPUT_GENOMES = config["INPUT_GENOMES"]
OUTPUT_TAXONOMY = config["OUTPUT_TAXONOMY"]
OUTPUT_TAXONOMY_NAME = config["OUTPUT_TAXONOMY_NAME"]

os.makedirs(OUTPUT_TAXONOMY, exist_ok=True)


###############################################################################
# rule all
###############################################################################
rule all:
    input:
        f"{OUTPUT_TAXONOMY}/{OUTPUT_TAXONOMY_NAME}_taxonomy_table.tsv",
        f"{OUTPUT_TAXONOMY}/taxonomy_creation.log"


###############################################################################
# CHECKPOINT: list_accessions
###############################################################################
checkpoint list_accessions:
    """
    Collect all directories in INPUT_GENOMES and write them to assemblies.txt.
    """
    input:
        directory = INPUT_GENOMES
    output:
        assemblies = f"{OUTPUT_TAXONOMY}/assemblies.txt"
    run:
        with open(output[0], 'w') as outfile:
            for dirname in os.listdir(input.directory):
                dirpath = os.path.join(input.directory, dirname)
                if os.path.isdir(dirpath) and dirname.startswith("GCA_"):
                    outfile.write(dirname + '\n')
                    print(f"Found assembly: {dirname}")

        if os.path.getsize(output.assemblies) == 0:
            raise ValueError(
                f"No directories found in {input.directory}. "
                "Ensure your genome directories are copied properly."
            )


###############################################################################
#  Function: get_accessions_from_checkpoint
#    Reads assemblies.txt AFTER the checkpoint is done
###############################################################################
def get_accessions_from_checkpoint(wildcards):
    # "get" the checkpoint object so we can see its outputs
    ck = checkpoints.list_accessions.get()
    with open(ck.output.assemblies) as f:
        return [line.strip() for line in f]


###############################################################################
# Rule: fetch_taxid (runs one job per accession)
###############################################################################
rule fetch_taxid:
    """
    For each discovered directory name (treated as an accession),
    fetch the TaxID and organism name from NCBI Datasets.
    """
    # Because we have a wildcard {accession}, Snakemake will create a job
    # for each item returned by get_accessions_from_checkpoint.
    input:
        # ensures we wait for the checkpoint
        assemblies=lambda wc: checkpoints.list_accessions.get().output.assemblies
    output:
        # one file per accession
        f"{OUTPUT_TAXONOMY}/results/{{accession}}_taxid.tsv"
    conda:
        f"{HOME_DIR}/env/ncbi-datasets.yaml"
    shell:
        """
        # Note: This shell snippet tries to fetch the TaxID for {wildcards.accession}.
        # If your directory name is truly an NCBI assembly accession (e.g. 'GCA_...'),
        # this works. If it's something else, you may need a fallback or custom logic.
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
        """


###############################################################################
# Rule: fetch_taxonomy (for each accession, uses the TaxID to get full taxonomy)
###############################################################################
rule fetch_taxonomy:
    """
    Given the TaxID for each accession, fetch the detailed taxonomy info.
    """
    input:
        taxid_info = f"{OUTPUT_TAXONOMY}/results/{{accession}}_taxid.tsv"
    output:
        f"{OUTPUT_TAXONOMY}/results/{{accession}}_taxonomy.tsv"
    conda:
        f"{HOME_DIR}/env/ncbi-datasets.yaml"
    shell:
        """
        accession=$(awk -F '\\t' '{{print $1}}' {input.taxid_info})
        tax_id=$(awk -F '\\t' '{{print $2}}' {input.taxid_info})

        echo "Processing accession: $accession"
        echo "Extracted TaxID: $tax_id"

        if [ "$tax_id" = "NA" ] || [ -z "$tax_id" ]; then
            echo -e "$accession\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA" > {output}
        else
            datasets summary taxonomy taxon $tax_id \\
            | jq -r '
                if .reports and (.reports | length > 0) then
                    [
                      "tax_id\\tscientific_name\\trank\\tkingdom\\tphylum\\tclass\\torder\\tfamily",
                      (.reports[] | if (.taxonomy // empty) then
                          [
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
        genome_file = lambda wc: checkpoints.list_accessions.get().output.assemblies,
        taxonomy_file = f"{OUTPUT_TAXONOMY}/results/{{accession}}_taxonomy_cleaned.tsv"
    output:
        temp(f"{OUTPUT_TAXONOMY}/organized/{{accession}}/linked.flag")
    run:
        import os
        import pandas as pd
        import glob

        # Validate input and find the correct genome file
        genome_files = glob.glob(f"{INPUT_GENOMES}/{wildcards.accession}/*_genomic.fna")
        
        if not genome_files:
            print(f"Genome file for {wildcards.accession} is missing or not a `.fna` file. Skipping.")
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
        expand(f"{OUTPUT_TAXONOMY}/organized/{{accession}}/linked.flag", accession=get_accessions_from_checkpoint)
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
        # This references the expansions from the checkpoint, so we
        # only run on the actual discovered directories.
        expand(
            f"{OUTPUT_TAXONOMY}/results/{{accession}}_taxonomy_cleaned.tsv",
            accession=get_accessions_from_checkpoint
        )
    output:
        f"{OUTPUT_TAXONOMY}/{OUTPUT_TAXONOMY_NAME}_taxonomy_table.tsv"
    run:
        dfs = [pd.read_csv(file, sep='\t') for file in input]
        final_df = pd.concat(dfs, ignore_index=True).drop_duplicates()
        final_df.to_csv(output[0], sep='\t', index=False)