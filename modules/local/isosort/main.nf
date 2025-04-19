process ISOSORT {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coreutils:8.25--1':
        'biocontainers/coreutils:8.25--1' }"

    input:
    path csv

    output:
    path "isochores.sorted.csv", emit: csv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    (head -n 1 ${csv} && tail -n +2 ${csv} | sort -t, -k1,1 -k4,4n) > isochores.sorted.csv


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sort: \$(sort --version | sed '1!d; s/.*sort (GNU coreutils) //')
    END_VERSIONS
    """

    stub:

    """
    touch isochores.sorted.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sort: \$(sort --version | sed '1!d; s/.*sort (GNU coreutils) //')
    END_VERSIONS
    """
}
