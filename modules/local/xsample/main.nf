process STAPLE_XSAMPLE {
    tag "cross-sample"
    container "ghcr.io/break-through-cancer/btc-containers/scverse@sha256:47a5a7292df74c7d4446609a5fae9676235292a60a1fa46ff02762cb3d10d0dc"


    input:
    path collected_items, stageAs: "?/*" // stage numbers the files as 0,1,2...

    output:
    path "reports/*.csv",               emit: reports,          optional: true //full csv reports
    path "reports/mqc/*mqc*",           emit: multiqc_files,    optional: true //multiqc reports for viewing and ai
    path "versions.yml",                emit: versions

    script:
    template 'xsample.py'
}