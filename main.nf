#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/assemblybacterias
========================================================================================
 nf-core/assemblybacterias Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/assemblybacterias
----------------------------------------------------------------------------------------
*/

log.info Headers.nf_core(workflow, params.monochrome_logs)

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////+
def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run nf-core/assemblybacterias --input 'samplesheet.csv' -profile docker --fasta GRCh39.fna.gz"
    log.info NfcoreSchema.params_help(workflow, params, json_schema, command)
    exit 0
}

////////////////////////////////////////////////////
/* --         VALIDATE PARAMETERS              -- */
////////////////////////////////////////////////////+
if (params.validate_params) {
    NfcoreSchema.validateParameters(params, json_schema, log)
}

//if (params.bacteria_database) { ch_kmerfinder_db = file(params.bacteria_database, checkIfExists: true)} else { exit 1, "kmerfinder Database file not specified!" }
//ch_kmerfinder_db = Channel.fromPath(params.bacteria_database, type:'dir')
ch_kmerfinder_db = file(params.bacteria_database)

if (params.bacteria_taxonomy) { ch_kmerfinder_taxonomy = file(params.bacteria_taxonomy, checkIfExists: true) } else { exit 1, "Kmerfinder taxonmy file does not exist" }


if (params.reference_ncbi_bacteria) {ch_reference_ncbi_bacteria = file(params.reference_ncbi_bacteria, checkIfExists: true)} else {exit 1, "Bacteria reference file does not exist"}
////////////////////////////////////////////////////
/* --     Collect configuration parameters     -- */
////////////////////////////////////////////////////

// SampleSheet input
if (params.input) { ch_input = file(params.input, checkIfExists: true) } else { exit 1, "Input samplesheet file not specified!" }

// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(', ')}"
}

// TODO nf-core: Add any reference files that are needed
// Configurable reference genomes
//
// NOTE - THIS IS NOT USED IN THIS PIPELINE, EXAMPLE ONLY
// If you want to use the channel below in a process, define the following:
//   input:
//   file fasta from ch_fasta
//
//params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
// if (params.fasta) { ch_fasta = file(params.fasta, checkIfExists: true) }

// Check AWS batch settings
if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, 'Specify correct --awsqueue and --awsregion parameters on AWSBatch!'
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, 'Outdir not on S3 - specify S3 Bucket to run on AWSBatch!'
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, 'Specify a local tracedir or run without trace! S3 cannot be used for tracefiles.'
}

// Stage config files
ch_multiqc_config = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs = file("$projectDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$projectDir/docs/images/", checkIfExists: true)

/*
 * Create a channel for input read files
 */
/*
if (params.input_paths) {
    if (params.single_end) {
        Channel
            .from(params.input_paths)
            .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true) ] ] }
            .ifEmpty { exit 1, 'params.input_paths was empty - no input files supplied' }
            .into { ch_read_files_fastp; ch_read_files_fastqc;  ch_read_files_trimming }
    } else {
        Channel
            .from(params.input_paths)
            .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true), file(row[1][1], checkIfExists: true) ] ] }
            .ifEmpty { exit 1, 'params.input_paths was empty - no input files supplied' }
            .into { ch_read_files_fastp; ch_read_files_fastqc;  ch_read_files_trimming }
    }
} else {
    Channel
        .fromFilePairs(params.input, size: params.single_end ? 1 : 2)
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.input}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --single_end on the command line." }
        .into { ch_read_files_fastp; ch_read_files_fastqc; ch_read_files_trimming }
}
*/




////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////
log.info NfcoreSchema.params_summary_log(workflow, params, json_schema)

// Header log info
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = workflow.runName
// TODO nf-core: Report custom parameters here
summary['Input']            = params.input
summary['Fasta Ref']        = params.fasta
if (params.gff)                        summary['GFF'] = params.gff
//summary['Data Type']        = params.single_end ? 'Single-End' : 'Paired-End'
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']   = params.awsregion
    summary['AWS Queue']    = params.awsqueue
    summary['AWS CLI']      = params.awscli
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Profile Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Profile Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config Profile URL']         = params.config_profile_url
summary['Config Files'] = workflow.configFiles.join(', ')
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
    summary['MultiQC maxsize']   = params.max_multiqc_email_size
}

// Check the hostnames against configured profiles
checkHostname()

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-assemblybacterias-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/assemblybacterias Workflow Summary'
    section_href: 'https://github.com/nf-core/assemblybacterias'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.indexOf('.csv') > 0) filename
                      else null
        }

    output:
    file 'software_versions_mqc.yaml' into ch_software_versions_yaml
    file 'software_versions.csv'

    script:
    // TODO nf-core: Get all tools to print their version number here
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    multiqc --version > v_multiqc.txt
    fastp --version > v_fastp.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

/*
 * PREPROCESSING: if specie genome fasta is provided check and uncompress genome fasta file if needed
 */
if (params.fasta){
    if (params.fasta.endsWith('.gz')) {
        process GUNZIP_FASTA {
            label 'error_retry'
            if (params.save_reference) {
                publishDir "${params.outdir}/genome", mode: params.publish_dir_mode
            }

            input:
            path fasta from params.fasta

            output:
            path "$unzip" into ch_fasta

            script:
            unzip = fasta.toString() - '.gz'
            """
            pigz -f -d -p $task.cpus $fasta
            """
        }
    } else {
        ch_fasta = file(params.fasta)
    }
} else {
    //See: https://nextflow-io.github.io/patterns/index.html#_optional_input
    ch_fasta = file('NO_FILE')
}

/*
 * PREPROCESSING: Uncompress gff annotation file
 */
if (params.gff) {
    file(params.gff, checkIfExists: true)
    if (params.gff.endsWith('.gz')) {
        process GUNZIP_GFF {
            label 'error_retry'
            if (params.save_reference) {
                publishDir "${params.outdir}/genome", mode: params.publish_dir_mode
            }

            input:
            path gff from params.gff

            output:
            path "$unzip" into ch_gff

            script:
            unzip = gff.toString() - '.gz'
            """
            pigz -f -d -p $task.cpus $gff
            """
        }
    } else {
        ch_gff = file(params.gff)
    }
} else {
    //See: https://nextflow-io.github.io/patterns/index.html#_optional_input
    ch_gff = file('NO_FILE')
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                     PARSE DESIGN FILE                               -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*
 * PREPROCESSING: Reformat samplesheet and check validity
 */
process CHECK_SAMPLESHEET {
    tag "$samplesheet"
    publishDir "${params.outdir}/", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.endsWith(".tsv")) "preprocess/sra/$filename"
                      else "pipeline_info/$filename"
                }

    input:
    path samplesheet from ch_input

    output:
    path "samplesheet.valid.csv" into ch_samplesheet_reformat
    path "sra_run_info.tsv" optional true

    script:  // These scripts are bundled with the pipeline, in nf-core/viralrecon/bin/
    run_sra = !params.skip_sra && !isOffline()
    """
    awk -F, '{if(\$1 != "" && \$2 != "") {print \$0}}' $samplesheet > nonsra_id.csv
    check_samplesheet.py nonsra_id.csv nonsra.samplesheet.csv
    awk -F, '{if(\$1 != "" && \$2 == "" && \$3 == "") {print \$1}}' $samplesheet > sra_id.list
    if $run_sra && [ -s sra_id.list ]
    then
        fetch_sra_runinfo.py sra_id.list sra_run_info.tsv --platform ILLUMINA --library_layout SINGLE,PAIRED
        sra_runinfo_to_samplesheet.py sra_run_info.tsv sra.samplesheet.csv
    fi
    if [ -f nonsra.samplesheet.csv ]
    then
        head -n 1 nonsra.samplesheet.csv > samplesheet.valid.csv
    else
        head -n 1 sra.samplesheet.csv > samplesheet.valid.csv
    fi
    tail -n +2 -q *sra.samplesheet.csv >> samplesheet.valid.csv
    """
}

// Function to get list of [ sample, single_end?, is_sra?, is_ftp?, [ fastq_1, fastq_2 ], [ md5_1, md5_2] ]
def validate_input(LinkedHashMap sample) {
    def sample_id = sample.sample_id
    def single_end = sample.single_end.toBoolean()
    def is_sra = sample.is_sra.toBoolean()
    def is_ftp = sample.is_ftp.toBoolean()
    def fastq_1 = sample.fastq_1
    def fastq_2 = sample.fastq_2
    def md5_1 = sample.md5_1
    def md5_2 = sample.md5_2

    def array = []
    if (!is_sra) {
        if (single_end) {
            array = [ sample_id, single_end, is_sra, is_ftp, [ file(fastq_1, checkIfExists: true) ] ]
        } else {
            array = [ sample_id, single_end, is_sra, is_ftp, [ file(fastq_1, checkIfExists: true), file(fastq_2, checkIfExists: true) ] ]
        }
    } else {
        array = [ sample_id, single_end, is_sra, is_ftp, [ fastq_1, fastq_2 ], [ md5_1, md5_2 ] ]
    }

    return array
}

/*
 * Create channels for input fastq files
 */
ch_samplesheet_reformat
    .splitCsv(header:true, sep:',')
    .map { validate_input(it) }
    .into { ch_reads_all
            ch_reads_sra }

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                     DOWNLOAD SRA FILES                              -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*
 * STEP 1: Download and check SRA data
 */
if (!params.skip_sra || !isOffline()) {
    ch_reads_sra
        .filter { it[2] }
        .into { ch_reads_sra_ftp
                ch_reads_sra_dump }

    process SRA_FASTQ_FTP {
        tag "$sample"
        label 'process_medium'
        label 'error_retry'
        publishDir "${params.outdir}/preprocess/sra", mode: params.publish_dir_mode,
            saveAs: { filename ->
                          if (filename.endsWith(".md5")) "md5/$filename"
                          else params.save_sra_fastq ? filename : null
                    }

        when:
        is_ftp

        input:
        tuple val(sample), val(single_end), val(is_sra), val(is_ftp), val(fastq), val(md5) from ch_reads_sra_ftp

        output:
        tuple val(sample), val(single_end), val(is_sra), val(is_ftp), path("*.fastq.gz") into ch_sra_fastq_ftp
        path "*.md5"

        script:
        if (single_end) {
            """
            curl -L ${fastq[0]} -o ${sample}.fastq.gz
            echo "${md5[0]}  ${sample}.fastq.gz" > ${sample}.fastq.gz.md5
            md5sum -c ${sample}.fastq.gz.md5
            """
        } else {
            """
            curl -L ${fastq[0]} -o ${sample}_1.fastq.gz
            echo "${md5[0]}  ${sample}_1.fastq.gz" > ${sample}_1.fastq.gz.md5
            md5sum -c ${sample}_1.fastq.gz.md5
            curl -L ${fastq[1]} -o ${sample}_2.fastq.gz
            echo "${md5[1]}  ${sample}_2.fastq.gz" > ${sample}_2.fastq.gz.md5
            md5sum -c ${sample}_2.fastq.gz.md5
            """
        }
    }

    process SRA_FASTQ_DUMP {
        tag "$sample"
        label 'process_medium'
        label 'error_retry'
        publishDir "${params.outdir}/preprocess/sra", mode: params.publish_dir_mode,
            saveAs: { filename ->
                          if (filename.endsWith(".log")) "log/$filename"
                          else params.save_sra_fastq ? filename : null
                    }

        when:
        !is_ftp

        input:
        tuple val(sample), val(single_end), val(is_sra), val(is_ftp) from ch_reads_sra_dump.map { it[0..3] }

        output:
        tuple val(sample), val(single_end), val(is_sra), val(is_ftp), path("*.fastq.gz") into ch_sra_fastq_dump
        path "*.log"

        script:
        prefix = "${sample.split('_')[0..-2].join('_')}"
        pe = single_end ? "" : "--readids --split-e"
        rm_orphan = single_end ? "" : "[ -f  ${prefix}.fastq.gz ] && rm ${prefix}.fastq.gz"
        """
        parallel-fastq-dump \\
            --sra-id $prefix \\
            --threads $task.cpus \\
            --outdir ./ \\
            --tmpdir ./ \\
            --gzip \\
            $pe \\
            > ${prefix}.fastq_dump.log
        $rm_orphan
        """
    }

    ch_reads_all
        .filter { !it[2] }
        .concat(ch_sra_fastq_ftp, ch_sra_fastq_dump)
        .set { ch_reads_all }
}

ch_reads_all
    .map { [ it[0].split('_')[0..-2].join('_'), it[1], it[4] ] }
    .groupTuple(by: [0, 1])
    .map { [ it[0], it[1], it[2].flatten() ] }
    .set { ch_reads_all }

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                     MERGE RESEQUENCED FASTQ                         -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*
 * STEP 2: Merge FastQ files with the same sample identifier
 */
process CAT_FASTQ {
    tag "$sample"

    input:
    tuple val(sample), val(single_end), path(reads) from ch_reads_all

    output:
    tuple val(sample), val(single_end), path("*.merged.fastq.gz") into ch_cat_fastqc,
                                                                       ch_cat_fastp

    script:
    readList = reads.collect{it.toString()}
    if (!single_end) {
        if (readList.size > 2) {
            def read1 = []
            def read2 = []
            readList.eachWithIndex{ v, ix -> ( ix & 1 ? read2 : read1 ) << v }
            """
            cat ${read1.sort().join(' ')} > ${sample}_1.merged.fastq.gz
            cat ${read2.sort().join(' ')} > ${sample}_2.merged.fastq.gz
            """
        } else {
            """
            ln -s ${reads[0]} ${sample}_1.merged.fastq.gz
            ln -s ${reads[1]} ${sample}_2.merged.fastq.gz
            """
        }
    } else {
        if (readList.size > 1) {
            """
            cat ${readList.sort().join(' ')} > ${sample}.merged.fastq.gz
            """
        } else {
            """
            ln -s $reads ${sample}.merged.fastq.gz
            """
        }
    }
}

/*
 * STEP 3 - FastQC
 */

process fastqc {
    tag "$name"
    label 'process_low'
    publishDir "${params.outdir}/fastqc", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      filename.indexOf('.zip') > 0 ? "zips/$filename" : "$filename"
        }

    input:
    tuple val(name), val(single_end), path(reads) from ch_cat_fastqc

    output:
    file '*_fastqc.{zip,html}' into ch_fastqc_results

    script:
    """
    fastqc --quiet --threads $task.cpus $reads
    """
}


/*
 * STEP 4: Fastp adapter trimming and quality filtering
 */

process FASTP {
        tag "$sample"
        label 'process_low'
        publishDir "${params.outdir}/preprocess/fastp", mode: params.publish_dir_mode,
            saveAs: { filename ->
                        if (filename.endsWith(".json")) filename
                        else if (filename.endsWith(".fastp.html")) filename
                        else if (filename.endsWith("_fastqc.html")) "fastqc/$filename"
                        else if (filename.endsWith(".zip")) "fastqc/zips/$filename"
                        else if (filename.endsWith(".log")) "log/$filename"
                        else params.save_trimmed ? filename : null
                    }
        input:
        tuple val(sample), val(single_end), path(reads) from ch_cat_fastp

        output:
        tuple val(sample), val(single_end), path("*.trim.fastq.gz") into ch_fastp_kmerfider,
                                                                    ch_fastp_unicycler
        path "*.json" into ch_fastp_mqc
        //path "*_fastqc.{zip,html}" into ch_fastp_fastqc_mqc
        path "*.{log,fastp.html}"
        path "*.fail.fastq.gz"

        script:
        // Added soft-links to original fastqs for consistent naming in MultiQC
        autodetect = single_end ? "" : "--detect_adapter_for_pe"
        """
        IN_READS='--in1 ${sample}.fastq.gz'
        OUT_READS='--out1 ${sample}.trim.fastq.gz --failed_out ${sample}.fail.fastq.gz'
        if $single_end; then
            [ ! -f  ${sample}.fastq.gz ] && ln -s $reads ${sample}.fastq.gz
        else
            [ ! -f  ${sample}_1.fastq.gz ] && ln -s ${reads[0]} ${sample}_1.fastq.gz
            [ ! -f  ${sample}_2.fastq.gz ] && ln -s ${reads[1]} ${sample}_2.fastq.gz
            IN_READS='--in1 ${sample}_1.fastq.gz --in2 ${sample}_2.fastq.gz'
            OUT_READS='--out1 ${sample}_1.trim.fastq.gz --out2 ${sample}_2.trim.fastq.gz --unpaired1 ${sample}_1.fail.fastq.gz --unpaired2 ${sample}_2.fail.fastq.gz'
        fi

        fastp \\
            \$IN_READS \\
            \$OUT_READS \\
            $autodetect \\
            --cut_front \\
            --cut_tail \\
            --cut_mean_quality $params.cut_mean_quality \\
            --qualified_quality_phred $params.qualified_quality_phred \\
            --unqualified_percent_limit $params.unqualified_percent_limit \\
            --length_required $params.min_trim_length \\
            --trim_poly_x \\
            --thread $task.cpus \\
            --json ${sample}.fastp.json \\
            --html ${sample}.fastp.html \\
            2> ${sample}.fastp.log

        """
}



/*
 * STEP 3 - Kmerfinder
 */
if (params.used_external_reference == false ){
    process kmerfinder {
        tag "$sample"
        label 'process_low'

        publishDir "${params.outdir}/kmerfinder/${sample}", mode: params.publish_dir_mode

        input:
        tuple val(sample), val(single_end), path(reads) from ch_fastp_kmerfider
        file kmerfinderDB from ch_kmerfinder_db
        file kmerfinderTAX from ch_kmerfinder_taxonomy

        output:
        //path "${sample}/*.txt" into ch_kmerfinder_results
        path "${sample}_results.txt" into ch_kmerfinder_results

        script:
        """
        IN_READS='-i ${sample}.trim.fastq.gz'
        if $single_end; then
            [ ! -f  ${sample}.fastq.gz ] && ln -s $reads ${sample}.fastq.gz
        else
            [ ! -f  ${sample}_1.trim.fastq.gz ] && ln -s ${reads[0]} ${sample}_1.trim.fastq.gz
            [ ! -f  ${sample}_2.trim.fastq.gz ] && ln -s ${reads[1]} ${sample}_2.trim.fastq.gz
            IN_READS='-i ${sample}_1.trim.fastq.gz  ${sample}_2.trim.fastq.gz'
        fi

        kmerfinder.py \\
        \$IN_READS -o ${sample} \\
        -db  $kmerfinderDB/bacteria.ATG \\
        -tax $kmerfinderTAX  -x
        mv ${sample}/results.txt ${sample}_results.txt
        """
    }
}

/*
 * STEP 4 - Find common reference from kmerfinder results
 */
if (params.used_external_reference == false ){
    process FIND_COMMON_REFERENCE{
        tag "Find common Reference"
        label 'process_low'
        publishDir "${params.outdir}/reference_download", mode: params.publish_dir_mode

        input:
        path ('kmerfinder_results/') from ch_kmerfinder_results.collect().ifEmpty([])
        // file reference_bacteria_file from ch_reference_ncbi_bacteria

        nucleotide_end='_genomic.fna.gz'
        protein_end='_protein.faa.gz'
        gff_end='_genomic.gff.gz'
        output:
        file 'references_found.tsv'
        file 'bacteria_id' into ch_findcomon_download

        script:
        """
        find_common_reference.py -d kmerfinder_results -o references_found.tsv
        bacteria_id=\$(head -n1 references_found.tsv | cut -f1 -d\$'\t')
        echo \$bacteria_id > 'bacteria_id'
        """

    }
}

process REFERENCE_DOWNLOAD{
    tag "Reference Dowload"
    label 'process_low'
    publishDir "${params.outdir}/reference_download", mode: params.publish_dir_mode

    input:
    file bacteria_file_id from ch_findcomon_download

    file reference_bacteria_file from ch_reference_ncbi_bacteria


    nucleotide_end='_genomic.fna.gz'
    protein_end='_protein.faa.gz'
    gff_end='_genomic.gff.gz'
    output:
    file 'REFERENCES/*_genomic.fna' into ch_reference_fna
    file 'REFERENCES/*_protein.faa' into ch_reference_protein
    file 'REFERENCES/*_genomic.gff' into ch_reference_gff
    // file 'references_found.tsv'


    script:
    """
    bacteriaID=\$(cat ${bacteria_file_id})
    ftp_path=\$(grep "\$bacteriaID" $reference_bacteria_file | cut -f20)
    echo \$ftp_path > path_ftp
    download_reference.py -url \$ftp_path -out_dir REFERENCES
    gunzip REFERENCES/*.gz
    """
}
/*
 * STEP 4 - Download reference L
 */
/*
process UNICYCLER {
	tag "$prefix"
    label 'process_low'
	publishDir path: { "${params.outdir}/unicycler" }, mode: 'copy'

	input:
	//set file(readsR1),file(readsR2) from ch_fastp_unicycler
    tuple val(sample), val(single_end), path(reads) from ch_fastp_unicycler
	output:
	file "${prefix}_assembly.fasta" into ch_unicycler_quast, ch_unicycler_prokka

	script:
	//prefix = readsR1.toString() - ~/(.R1)?(_1)?(_R1)?(_trimmed)?(_paired)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
	"""
	unicycler --threads ${task.cpus} -1 ${sample}_1.trim.fastq.gz -2${sample}_2.trim.fastq.gz  -o .

	"""
}
*/

/*
 * STEP 3 - Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    file output_docs from ch_output_docs
    file images from ch_output_docs_images

    output:
    file 'results_description.html'

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}


/*
 * STEP 2 - MultiQC
 */

process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: params.publish_dir_mode

    input:
    file (multiqc_config) from ch_multiqc_config
    file (mqc_custom_config) from ch_multiqc_custom_config.collect().ifEmpty([])
    // TODO nf-core: Add in log files from your new processes for MultiQC to find!
    file ('fastqc/*') from ch_fastqc_results.collect().ifEmpty([])
    file ('software_versions/*') from ch_software_versions_yaml.collect()
    file workflow_summary from ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")

    output:
    file "*multiqc_report.html" into ch_multiqc_report
    file "*_data"
    file "multiqc_plots"

    script:
    rtitle = ''
    rfilename = ''
    if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
        rtitle = "--title \"${workflow.runName}\""
        rfilename = "--filename " + workflow.runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report"
    }
    custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
    // TODO nf-core: Specify which MultiQC modules to use with -m for a faster run time
    """
    multiqc -f $rtitle $rfilename $custom_config_file .
    """
}
/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/assemblybacterias] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[nf-core/assemblybacterias] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // TODO nf-core: If not using MultiQC, strip out this code (including params.max_multiqc_email_size)
    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = ch_multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList) {
                log.warn "[nf-core/assemblybacterias] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nf-core/assemblybacterias] Could not attach MultiQC report to summary email"
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$projectDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$projectDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, projectDir: "$projectDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$projectDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[nf-core/assemblybacterias] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
            if ( mqc_report.size() <= params.max_multiqc_email_size.toBytes() ) {
              mail_cmd += [ '-A', mqc_report ]
            }
            mail_cmd.execute() << email_html
            log.info "[nf-core/assemblybacterias] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[nf-core/assemblybacterias]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[nf-core/assemblybacterias]${c_red} Pipeline completed with errors${c_reset}-"
    }

}

workflow.onError {
    // Print unexpected parameters - easiest is to just rerun validation
    NfcoreSchema.validateParameters(params, json_schema, log)
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = 'hostname'.execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "${c_red}====================================================${c_reset}\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "${c_red}====================================================${c_reset}\n"
                }
            }
        }
    }
}

// Function to check if running offline
def isOffline() {
    try {
        return NXF_OFFLINE as Boolean
    }
    catch( Exception e ) {
        return false
    }
}
