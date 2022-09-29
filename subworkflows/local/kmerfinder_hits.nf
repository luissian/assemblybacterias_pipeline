//
// Kmerfinder and download common reference
//


include {KMERFINDER             } from '../../modules/local/kmerfinder'
include {HITS_KMERFINDER    } from '../../modules/local/hits_kmerfinder'
// include {DOWNLOAD_REFERENCE } from ''

workflow KMERFINDER_HITS {
    take:
    reads           // channel: [ val(meta), [ reads ] ]
    kmer_bacteria_db

    main:
    reads
        .dump(tag:'kmerfinder')
        .set { reads_kmer }

    //kmerfinder_result  = Channel.empty()
    KMERFINDER ( reads_kmer , kmer_bacteria_db)
    kmerfinder_result = KMERFINDER.out.kmer_result.collect()

    HITS_KMERFINDER { kmerfinder_result }

    kmer_result    = KMERFINDER.out.kmer_result
    versions = KMERFINDER.out.versions.first()

    emit:
    kmer_result
    versions
}
