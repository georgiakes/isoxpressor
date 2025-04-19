
process ISO2GTF {
    label 'process_single'

     conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:4.1.3--1':
        'biocontainers/gawk:4.1.3--1' }"

    input:

    path isochores

    output:
    path "*.gtf", emit: gtf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """

    # Read CSV, skipping the header (assuming comma-separated values)
    awk -F, 'NR > 1 {print \$1, \$4, \$5, \$2 ,NR-1}' ${isochores} | while read chr start end class ID; do
        # Format GTF entry
        echo -e "\${chr}\tiso2gtf\texon\t\${start}\t\${end}\t.\t+\t.\tgene_id \\"\${ID}\\"; iso_fam \\"\${class}\\";" >> all_isochores.gtf
    done


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''

    """
    touch all_isochores.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}
