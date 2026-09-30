process CREATE_TEST_ADATA {
    tag "create_test_adata"
    label 'process_single'
    container 'ghcr.io/break-through-cancer/btc-containers/scverse@sha256:47a5a7292df74c7d4446609a5fae9676235292a60a1fa46ff02762cb3d10d0dc'

    input:
        val num_adatas
        val with_metadata
        val seeds

    output:
        path "*adata.h5ad",     emit: adata
        path "samplesheet.csv", emit: samplesheet

    script:
    template 'create_adata.py'

}