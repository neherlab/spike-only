data_source = config["data_source"]

if data_source == "open":
    use_open = True
    prefix = "https://data.nextstrain.org/files/ncov/open/100k/"
else:
    use_open = False
    prefix = "s3://nextstrain-ncov-private/100k/"


rule build:
    input:
        "auspice/spike-only.json",


rule download_sequences:
    output:
        sequences="data/sequences.fasta.xz",
    params:
        url=prefix + "sequences.fasta.xz",
        open=use_open,
    shell:
        """
        if [ {params.open} == True ]; then
            curl {params.url} -o {output.sequences}
        else
            aws s3 cp {params.url} {output.sequences}
        fi
        """


rule download_metadata:
    output:
        metadata="data/metadata.tsv.xz",
    params:
        url=prefix + "metadata.tsv.xz",
        open=use_open,
    shell:
        """
        if [ {params.open} == True ]; then
            curl {params.url} -o {output.metadata}
        else
            aws s3 cp {params.url} {output.metadata}
        fi
        """


rule download_lat_longs:
    output:
        "builds/lat_longs.tsv",
    params:
        url="https://raw.githubusercontent.com/nextstrain/ncov/master/defaults/lat_longs.tsv",
    shell:
        """
        curl {params.url} | \
        sed "s/North Rhine Westphalia/North Rhine-Westphalia/g" | \
        sed "s/Baden-Wuerttemberg/Baden-Wurttemberg/g" \
        > {output}
        """


rule align:
    input:
        sequences="data/sequences.fasta.xz",
        reference="resources/reference.fasta",
        annotation="resources/annotation.gff3",
    output:
        alignment="builds/aligned.fasta",
        translation="builds/S.fasta",
        tsv="builds/nextclade.tsv",
    params:
        translation_template=lambda w: "builds/{cds}.fasta",
    shell:
        """
        xzcat {input.sequences} \
        | nextclade3 run --input-ref {input.reference} \
            --input-annotation {input.annotation} \
            --min-seed-cover 0.01 \
            --output-tsv {output.tsv} \
            --output-fasta {output.alignment} \
            --include-reference \
            --output-translations {params.translation_template}
        """


rule mask_terminals:
    input:
        alignment="builds/aligned.fasta",
    output:
        alignment="builds/aligned_masked.fasta",
    shell:
        """
        python3 scripts/mask-alignment.py \
            --alignment {input.alignment} \
            --mask-terminal-gaps \
            --output {output.alignment}
        """


rule index:
    input:
        alignment="builds/aligned_masked.fasta",
    output:
        index="builds/index.tsv",
    shell:
        """
        augur index \
            --sequences {input.alignment} \
            --output {output.index}
        """


rule filter_unknowns:
    """
    # Load into pandas
    # Filter out sequences where sum of columns ?, N, other_IUPAC is >100
    """
    input:
        index="builds/index.tsv",
    output:
        index="builds/index_filtered.tsv",
    run:
        import pandas as pd

        df = pd.read_csv(input.index, sep="\t", low_memory=False)
        df["unknowns"] = df["?"] + df["N"] + df["other_IUPAC"]
        df = df[df.unknowns <= 200]
        df.to_csv(output.index, sep="\t", index=False)


rule filter_bad:
    """
    Using tsv-filter because augur filter struggles with numeric columns
    that have missing data.
    Strains with `'` need to be filtered out because Bio.Phylo does not
    parse newicks with them correctly, see
    https://github.com/biopython/biopython/issues/4536
    https://github.com/biopython/biopython/issues/4537
    """
    input:
        metadata="data/metadata.tsv.xz",
    output:
        metadata="builds/metadata_filtered.tsv.zst",
    shell:
        """
        xzcat {input.metadata} \
        | tsv-filter -H \
            --str-not-in-fld "strain:'" \
            --str-not-in-fld "strain:env" \
            --is-numeric reversion_mutations \
            --is-numeric rare_mutations \
            --is-numeric QC_overall_score \
            --is-numeric clock_deviation \
            --lt reversion_mutations:2 \
            --lt rare_mutations:20 \
            --str-eq QC_frame_shifts:good \
            --str-eq QC_stop_codons:good \
            --str-eq QC_mixed_sites:good \
            --str-eq QC_snp_clusters:good \
            --str-ne QC_missing_data:bad \
            --lt QC_overall_score:100 \
            --lt clock_deviation:18 \
            --gt clock_deviation:-15 \
        | zstd -c > {output.metadata}
        """


rule subsample:
    input:
        metadata="builds/metadata_filtered.tsv.zst",
        index="builds/index_filtered.tsv",
    output:
        metadata="builds/metadata_subsampled.tsv",
    params:
        number=2000,
    shell:
        """
        augur filter \
            --metadata {input.metadata} \
            --sequence-index {input.index} \
            --subsample-max-sequences {params.number} \
            --subsample-seed 42 \
            --group-by country year month \
            --output-metadata {output.metadata}
        """


rule get_strains:
    input:
        metadata="builds/metadata_subsampled.tsv",
    output:
        strains="builds/strains.txt",
    shell:
        """
        tsv-select -H -f strain {input.metadata} >{output.strains}
        echo "Wuhan-Hu-1/2019" >>{output.strains}
        """


rule get_sampled_sequences:
    input:
        sequences="builds/aligned_masked.fasta",
        strain_names="builds/strains.txt",
    output:
        sequences="builds/sequences_sample.fasta",
    shell:
        """
        seqkit grep -f {input.strain_names} <{input.sequences} >{output.sequences}
        """


rule tree:
    input:
        alignment=rules.get_sampled_sequences.output.sequences,
    output:
        tree="builds/tree_raw.nwk",
    params:
        args="'-ninit 2 -n 2 -czb -T 4'",
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --tree-builder-args {params.args} \
            --output {output.tree}
        """


rule fix_iqtree:
    input:
        tree="builds/tree_raw.nwk",
        alignment=rules.get_sampled_sequences.output.sequences,
    output:
        tree="builds/tree_fixed.nwk",
    shell:
        """
        python3 scripts/fix_tree.py \
            --alignment {input.alignment} \
            --input-tree {input.tree} \
            --root "Wuhan-Hu-1/2019" \
            --output {output.tree}
        """


rule refine:
    input:
        tree="builds/tree_fixed.nwk",
        alignment=rules.get_sampled_sequences.output.sequences,
        metadata="builds/metadata_subsampled.tsv",
    output:
        node_data="builds/branch_lengths.json",
        tree="builds/tree.nwk",
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --keep-root \
            --stochastic-resolve \
            --timetree \
            --use-fft \
            --divergence-units mutations \
            --output-node-data {output.node_data} \
            --output-tree {output.tree}
        """


rule ancestral:
    input:
        tree="builds/tree.nwk",
        alignment=rules.get_sampled_sequences.output.sequences,
        annotation="resources/annotation.gb",
        translations="builds/S.fasta",
    output:
        node_data="builds/muts.json",
    params:
        inference="joint",
        translations="builds/%GENE.fasta",
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference} \
            --genes S \
            --annotation {input.annotation} \
            --translations {params.translations} \
            2>&1 | tee {log}
        """


rule colors:
    input:
        ordering="resources/color_ordering.tsv",
        color_schemes="resources/color_schemes.tsv",
        metadata="builds/metadata_subsampled.tsv",
    output:
        colors="builds/colors.tsv",
    shell:
        """
        python3 scripts/assign-colors.py \
            --ordering {input.ordering} \
            --color-schemes {input.color_schemes} \
            --output {output.colors} \
            --metadata {input.metadata} 2>&1
        """


rule auspice_config_to_json:
    input:
        "resources/auspice_config.yaml",
    output:
        "builds/auspice_config.json",
    shell:
        """
        yq -o json {input} >{output}
        """


rule export:
    input:
        tree="builds/tree.nwk",
        node_data="builds/branch_lengths.json",
        ancestral="builds/muts.json",
        auspice_config="builds/auspice_config.json",
        lat_longs=rules.download_lat_longs.output,
        metadata="builds/metadata_subsampled.tsv",
        colors="builds/colors.tsv",
        description="resources/description.md",
    output:
        auspice_json="auspice/spike-only.json",
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --node-data {input.node_data} {input.ancestral} \
            --colors {input.colors} \
            --include-root-sequence-inline \
            --auspice-config {input.auspice_config} \
            --lat-longs {input.lat_longs} \
            --output {output.auspice_json} \
            --metadata {input.metadata} \
            --description {input.description}
        """


rule deploy:
    input:
        "auspice/spike-only.json",
    output:
        touch("builds/spike-only.upload.done"),
    shell:
        """
        nextstrain remote upload nextstrain.org/groups/neherlab/ncov/spike-only {input} 2>&1
        """
