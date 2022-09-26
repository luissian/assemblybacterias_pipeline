process KMERFINDER_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::kmerfinder=3.0.2" : null)
    container = "kmerfinder_v3.0.2"
    //container = "singularity-kmerFinder.3.0.2"

    input:
    tuple val(meta), path(reads)
    path kmerfinderDB


    output:
    path "kmerfinder/*"  , emit: kmerfinder_dir
    path "versions.yml", emit: versions

script:
    single_end = false
    samplename_dir = "${reads}"
    // in_reads = single_end ? "${reads}" : "${reads[0]} ${reads[1]}"
    kmerfinderDB = params.kmer_bacteria_db


    """
    kmerfinder.py \\
    --infile $reads \\
    --output_folder kmerfinder \\
    --db_path $kmerfinderDB/bacteria.ATG \\
    -tax  $kmerfinderDB/bacteria.name \\
    -x


    cat <<-END_VERSIONS > versions.yml
    ${task.process}":
        kmerfincer: \$(echo \$(kmerfinder -v 2>&1))
    END_VERSIONS
    """
}
