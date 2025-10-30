

process PARABRICKS_MARKDUP {
    tag "${meta.id}"
    label 'process_high'
    label 'process_gpu'
    // needed by the module to work properly can be removed when fixed upstream - see: https://github.com/nf-core/modules/issues/7226
    stageInMode 'copy'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1"

    input:
    tuple val(meta), path(input)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.markdups.bam"),  emit: bam
    path "versions.yml",             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit(1, "Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def num_gpus = task.accelerator ? "--num-gpus ${task.accelerator.request}" : ''

    """

    pbrun \\
        sortbam \\
        --in-bam ${input} \\
        --ref ${fasta} \\
        --out-bam ${prefix}.sorttmp.bam \\
        ‑‑gpusort \\
        ‑‑sort‑order queryname \\
        ${num_gpus} \\
        ${args}

    pbrun \\
        markdups \\
        --gpuwrite \\
        --gpusort  \\
        --in-bam ${prefix}.sorttmp.bam \\
        --ref ${fasta} \\
        --out-bam ${prefix}.markdups.bam \\
        ‑‑out‑duplicate‑metrics ${prefix}.markDuplicates.metrics.txt \\
        ${num_gpus} \\
        ${args2}

    rm -f ${prefix}.sorttmp.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parabricks: \$(parabricks --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.markdups.bam

    # Capture the full version output once and store it in a variable
    pbrun_version_output=\$(pbrun markdups --version 2>&1)

    # Generate compatible_versions.yml
    cat <<EOF > compatible_versions.yml
    "${task.process}":
        pbrun_version: \$(echo "\$pbrun_version_output" | grep "pbrun:" | awk '{print \$2}')
        compatible_with:
        \$(echo "\$pbrun_version_output" | awk '/Compatible With:/,/^---/{ if (\$1 ~ /^[A-Z]/ && \$1 != "Compatible" && \$1 != "---") { printf "  %s: %s\\n", \$1, \$2 } }')
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS

    """
}
