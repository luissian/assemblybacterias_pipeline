//
// Kmerfinder and download common reference
//


include {KMERFINDER_RUN         } from '../../modules/local/kmerfinder'
// include {DOWNLOAD_REFERENCE } from ''

workflow KMERFINDER_CONTAMINATION {
    take:
    reads           // channel: [ val(meta), [ reads ] ]
    kmer_bacteria_db


    main:
    reads
        .dump(tag:'kmerfinder')
        .set { reads_kmer }


    KMERFINDER_RUN ( reads_kmer , kmer_bacteria_db)



    emit:

    versions = KMERFINDER_RUN.out.versions
}
