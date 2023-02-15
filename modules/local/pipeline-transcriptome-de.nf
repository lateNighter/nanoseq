// input_file in samplesheet can be path to directory with BAM file if previously aligned. 
//then we take some *flag* number of treated and control samples (group - treated/control or num_control_samples in config, input_file - bam, transcriptome - fasta)
// + name of annotation file like "Homo_sapiens.GRCh38.102.gtf" (contents aren't used)
// ch_sample // [ sample->"{}_R{}".format(group,replicate), barcode, fasta, gtf, is_transcripts, annotation_str ]

process PIPELINE_TRANSCRIPTOME_DE {
    label 'process_medium'

    container "tmp_ptd:latest"
    containerOptions "-u \$(id -u):\$(id -g)"
    
    input:
    tuple path(fasta), path(gtf)
    val(control_samples)
    val(treated_samples)
    path(samples_list)
    path(script_dir)

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

    script:
    """
    workdir=\$PWD
    cd ptd/
    snakemake --use-conda -j $task.cpus all > \$workdir/tmp.txt
    
    """
    // snakemake --use-conda -j $task.cpus all --config control_samples=$control_samples treated_samples=$treated_samples annotation=$gtf transcriptome=$fasta
    // snakemake -n
    //echo "skdnvszf" > tmp.txt
    // -s /tmp/repo/pipeline-transcriptome-de/Snakefile

    //chmod a+rwX /tmp/repo/

    // workdir=\$PWD
    // cd ptd/
    // snakemake --use-conda -j $task.cpus all > \$workdir/tmp.txt

    // cd ptd/.snakemake/conda
    // ls -lash > \$workdir/tmp.txt
    
}