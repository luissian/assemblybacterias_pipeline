process HITS_KMERFINDER {

    label 'process_low'

    input:
    path kmerfinder_result

    output:
    path "*.json", emit: kmerfinder_hits

    script:
    """
    kmerfinder_hits.py --path . --output .
    """
}
