process ICHORCNA_GENERATE_PON {

label "process_low"
tag "Generate panel of normals"

container "quay.io/dincalcilab/ichorcna:0.4.0-2f8cfc8"

input:
    path wigfiles
    path gc_wig
    path map_wig
    path centromere
    path reptime_file

output:
    path "*.rds"       , emit: pon_file
    path "*.txt"       , emit: txt_file
    path "versions.yml", emit: versions

when:
    task.ext.when == null || task.ext.when

script:

def args = task.ext.args ?: ''
def prefix = task.ext.prefix ?: "PoN"
def centro = params.centromere ? "--centromere ${params.centromere}" : ''
def reptime = reptime_file ? "--repTimeWig ${reptime_file}": ''
def ichorcna_script = "createPanelOfNormals.R"

def VERSION = '0.3.2' // Version information not provided by tool on CLI

"""

echo "${wigfiles.join("\n")}" > PoN_wigfiles.txt

${ichorcna_script}\\
    --filelist PoN_wigfiles.txt \\
    --gcWig ${gc_wig} \\
    --mapWig ${map_wig} \\
    ${centro} \\
    ${args} \\
    ${reptime} \\
    --outfile ${prefix}

rm PoN_wigfiles.txt

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    ichorcna: $VERSION
END_VERSIONS
"""

}
