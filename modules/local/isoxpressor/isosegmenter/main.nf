
process ISOXPRESSOR_ISOSEGMENTER {
    tag "$chrom"
    label 'process_single'

    // container "isoxpressor:latest"


    input:
    // TODO nf-core: Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/modules/nf-core/bwa/index/main.nf
    // TODO nf-core: Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    tuple val(chrom) , path(fasta)

    output:
    path "isochores/*.csv" , emit: isochores
    path "graphs/*.jpg"    , emit: graphs
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${chrom}"

    """
    mkdir -p isochores graphs
    python /app/IsoXpressor/isoSegmenter/scripts/isoSegmenter.py \
        -i ${fasta} \
        --y_min 1 \
        --y_max 100 \
        -g graphs/${prefix}.jpg \
        --window_size ${params.window_size} \
        -o isochores/${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isoxpressor: 1.0
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${chrom}"

    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isoxpressor: 1.0
    END_VERSIONS
    """
}
