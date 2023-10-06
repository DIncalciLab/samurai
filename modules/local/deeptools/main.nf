process BAMPE_FRAGMENTSIZE {

    label 'process_medium'
    tag "Fragment size graph generation"
    container 'quay.io/biocontainers/deeptools:3.5.1--py_0'

    input:
       val prefix
       path(bamfiles)
       path(bamindexes)
    output:
        path "*_fragment_size_histogram.svg", emit: histogram
        path "*_fragment_size_table.txt", emit: table
        path "*_fragment_size_raw.txt", emit: raw_table
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:

        def filelist = bamfiles.join(" ")
        def args = task.ext.args ?: ''

        """
        bamPEFragmentSize \\
            --bamfiles ${filelist} \\
            -p ${task.cpus} \\
            --histogram ${prefix}_fragment_size_histogram.svg \\
            --table ${prefix}_fragment_size_table.txt \\
            --outRawFragmentLengths  ${prefix}_fragment_size_raw.txt \\
            --plotFileFormat svg \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bamPEFragmentSize: \$( bamPEFragmentSize --version | sed -e "s/bamPEFragmentSize //g" )
        END_VERSIONS
        """

}
