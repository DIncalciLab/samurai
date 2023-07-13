process SAMBAMBA_FILTER {

    //module 'apptainer/1.0.2'
    container "quay.io/biocontainers/sambamba:0.8.2--h98b6b92_2"

    label 'process_medium'
    tag "Filter fragment BAM file of ${meta.id}"

    input:
        tuple val(meta), path(bamfile), path(bamindex)
    output:
        tuple val(meta), path("*.bam"), path("*.bam.bai"), emit: filtered_bam
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    //sample_type = meta.sample_type ? "_${meta.sample_type}" : ''
    //destfile = "${meta.samplename}${sample_type}_filtered.bam"

    """

    sambamba view \\
        -h \\
        -t ${task.cpus} \\
        -f bam \\
        $args \\
        ${bamfile} \\
        -o ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sambamba: \$(sambamba -q --version 2>&1 >/dev/null | awk 'NR==2 {print \$2}')
    END_VERSIONS
    """

}