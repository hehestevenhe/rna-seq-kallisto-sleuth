log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("tidyverse")
library("sleuth")

so <- sleuth_load(snakemake@input[["so"]])

top_n <- -strtoi(snakemake@params["top_n"])

results <- read_tsv(snakemake@input[["transcripts"]])

dir.create( snakemake@output[[1]] )

# Create plots for top n transcripts
top_transcripts <- results %>%    
    filter(qval <= snakemake@params[["fdr"]]) %>%        
    top_n(top_n,qval) %>%
    select(., c("ext_gene","target_id"))
    drop_na()

for (i in nrow(top_transcripts)) {
    plot_bootstrap(so, top_transcripts[i,"target_id"], color_by = snakemake@params[["color_by"]], units = "tpm")
    transcript_save <- str_replace_all(top_transcript[i, "target_id"], ":", "_")
    ggsave(file = str_c(snakemake@output[[1]], "/", top_transcript[i,"ext_gene"], ".", transcript_save, ".", snakemake@wildcards[["model"]] , ".bootstrap.pdf"))
    }

# Create plots for genes of interest
if ( !is.null(snakemake@params[["genes"]]) ) {
  genes <- tibble( ext_gene = snakemake@params[["genes"]]) %>%
    distinct(ext_gene)
} else {
  # "genes" is null, if the list provided in config.yaml is empty
  genes <- tibble( ext_gene = character() )
}

for (gene in genes){
    transcripts <- results[,!(names(results) %in% "canonical")]  %>%   # Removes the canonical column, non-canonical transcripts are
        filter(ext_gene == gene) %>%                                   # otherwise not visualised due to the drop_na()
        drop_na() %>%
        pull(target_id)

    if ( length( transcripts > 0 ) ) {
        for (transcript in transcripts) {
            plot_bootstrap(so, transcript, color_by = snakemake@params[["color_by"]], units = "tpm")
            ggsave(file = str_c(snakemake@output[[1]], "/", gene, ".", str_replace_all(transcript, ":", "_"), ".", snakemake@wildcards[["model"]] , ".bootstrap.pdf"))
        }
    }
}
