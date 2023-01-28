// input_file in samplesheet can be path to directory with BAM file if previously aligned. 
//then we take some *flag* number of treated and control samples (group - treated/control or num_control_samples in config, input_file - bam, transcriptome - fasta)
// + name of annotation file like "Homo_sapiens.GRCh38.102.gtf" (contents aren't used)
// TODO: check samplesheet for the above rules and !skip_transcriptome-de
// ch_sample // [ sample->"{}_R{}".format(group,replicate), barcode, fasta, gtf, is_transcripts, annotation_str ]

process PIPELINE_TRANSCRIPTOME_DE {
    label 'process_medium'

    container "docker.io/arrid1/pipeline-transcriptome-de:latest"
    //TODO (from bambu)
    input:
    tuple path(fasta), path(gtf)
    val(control_samples)
    val(treated_samples)
    path(control_list)
    path(treated_list)

    output:
    path "tmp.txt"          , emit: ch_tmp_txt
    // path "counts/*"          , emit: ch_ptd_counts
    // path "merged/all_counts.tsv"    , emit: ch_ptd_all_counts
    // path "merged/all_counts_filtered.tsv"    , emit: ch_ptd_all_counts_filtered
    // path "merged/all_gene_counts.tsv"    , emit: ch_ptd_all_gene_counts
    // path "de_analysis/coldata.tsv"    , emit: ch_ptd_coldata
    // path "de_analysis/de_params.tsv"    , emit: ch_ptd_de_params
    // path "de_analysis/results_dge.tsv"    , emit: ch_ptd_results_dge_tsv
    // path "de_analysis/results_dge.pdf"    , emit: ch_ptd_results_dge_pdf
    // path "de_analysis/results_dtu_gene.tsv"    , emit: ch_ptd_results_dtu_gene
    // path "de_analysis/results_dtu_transcript.tsv"    , emit: ch_ptd_results_dtu_transcript
    // path "de_analysis/results_dtu.pdf"    , emit: ch_ptd_results_dtu
    // path "de_analysis/results_dtu_stageR.tsv"    , emit: ch_ptd_results_dtu_stageR
    // path "de_analysis/dtu_plots.pdf"    , emit: ch_ptd_dtu_plots
    
    // path "versions.yml"             , emit: versions

    script:
    """
    echo "skdnvszf" > tmp.txt
    """
    // snakemake --use-conda -j $task.cpus all --config control_samples=$control_samples treated_samples=$treated_samples annotation=$gtf transcriptome=$fasta
    // snakemake -n

    // run_bambu.r \\
    //     --tag=. \\
    //     --ncore=$task.cpus \\
    //     --annotation=$gtf \\
    //     --fasta=$fasta $bams
    // cat <<-END_VERSIONS > versions.yml
    // "${task.process}":
    //     r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    //     bioconductor-bambu: \$(Rscript -e "library(bambu); cat(as.character(packageVersion('bambu')))")
    //     bioconductor-bsgenome: \$(Rscript -e "library(BSgenome); cat(as.character(packageVersion('BSgenome')))")
    // END_VERSIONS
}