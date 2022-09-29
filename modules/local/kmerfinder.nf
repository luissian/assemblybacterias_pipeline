process KMERFINDER {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::kmerfinder=3.0.2" : null)
    container = "kmerfinder_v3.0.2"
    //container = "singularity-kmerFinder.3.0.2"

    input:
    tuple val(meta), path(reads)
    path kmerfinderDB


    output:
    path "${samplename}*"  , emit: kmer_result
    path "versions.yml", emit: versions

    script:
    // single_end = false
    samplename = "${meta.id}"
    // in_reads = single_end ? "${reads}" : "${reads[0]} ${reads[1]}"
    kmerfinderDB = params.kmer_bacteria_db


    """
    kmerfinder.py \\
    --infile $reads \\
    --output_folder ./ \\
    --db_path $kmerfinderDB/bacteria.ATG \\
    -tax  $kmerfinderDB/bacteria.name \\
    -x

    mv results.txt ${samplename}_results.txt
    mv results.spa ${samplename}_results.spa
    mv data.json ${samplename}_data.json

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        kmerfinder: 3.0.2
    END_VERSIONS
    """
}
