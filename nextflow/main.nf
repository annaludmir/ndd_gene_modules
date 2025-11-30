params.config = "config.yaml"
params.genes  = "gene_list.txt"

workflow {

    Channel
        .value(params.config)
        .set { cfg_ch }

    // STEP 1 — Calculate GES scores
    ges_results_ch = ges_score(cfg_ch)

    // STEP 2 — Choose enrichment mode based on config
    if (check_mode(params.config) == "gsea" || check_mode(params.config) == "both") {
        gsea_results_ch = gsea(ges_results_ch)
        plots_gsea(gsea_results_ch)
    }

    if (check_mode(params.config) == "deseq" || check_mode(params.config) == "both") {
        deseq_results_ch = deseq(params.genes, cfg_ch)
        plots_deseq(deseq_results_ch)
    }

}

process ges_score {
    publishDir "results/ges", mode: "copy"

    input:
    path config

    output:
    path "ges_out", emit: ges

    script:
    """
    python3 - << 'EOF'
    from ges_score_calculation import run_ges_pipeline
    run_ges_pipeline("$config")
    EOF

    mkdir ges_out
    cp -r results/ges/* ges_out/
    """
}


process gsea {
    publishDir "results/gsea", mode: "copy"

    input:
    path ges_results

    output:
    path "gsea_out", emit: gsea

    script:
    """
    python3 - << 'EOF'
    from enrichment_pipeline_for_gene_list import run_gsea_pipeline
    run_gsea_pipeline("$params.config", "$params.genes")
    EOF

    mkdir gsea_out
    cp -r results/gsea/* gsea_out/
    """
}

process deseq {
    publishDir "results/deseq", mode: "copy"

    input:
    path genesFile
    path config

    output:
    path "deseq_out", emit: deseq

    script:
    """
    python3 - << 'EOF'
    from enrichment_pipeline_for_gene_list import run_deseq_pipeline
    run_deseq_pipeline("$config", "$genesFile")
    EOF

    mkdir deseq_out
    cp -r results/deseq/* deseq_out/
    """
}

process plots_gsea {
    publishDir "results/plots_gsea", mode: "copy"

    input:
    path results

    script:
    """
    python3 - << 'EOF'
    from create_figs_ges import plot_bar_chart
    plot_bar_chart(results_folder="results/gsea")
    EOF
    """
}

process plots_deseq {
    publishDir "results/plots_deseq", mode: "copy"

    input:
    path results

    script:
    """
    python3 - << 'EOF'
    from create_figs_deseq import plot_bar_chart
    plot_bar_chart(results_folder="results/deseq")
    EOF
    """
}

def check_mode(cfg){
    def parsed = new File(cfg).text.readLines().find{ it.contains("analysis_mode") }
    return parsed.split(":")[1].trim()
}
