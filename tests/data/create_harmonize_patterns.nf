process CREATE_HARMONIZE_PATTERNS {
    tag "create_patterns_${test_type}"
    label 'process_single'
    container 'ghcr.io/break-through-cancer/btc-containers/scverse@sha256:ed44380c6e6e73fc575b743eba864941c26880053e50c3d70b5c9bfc526c0520'

    input:
        val test_type

    output:
        path "patterns/*.csv", emit: patterns

    script:
    template 'create_harmonize_patterns.py'
}
