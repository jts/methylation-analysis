source("~/simpsonlab/users/jsimpson/code/methylation-analysis/methylation_plots.R")

quick_plot <- function(data, kmer, title = "") 
{ 
    ggplot(data[c(kmer, gsub("M", "C", kmer))], aes(level_mean, fill=model)) + 
        geom_histogram(aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]), binwidth = 0.25, alpha=0.5, position="identity") + 
        facet_grid(dataset ~ model, scale="free") + 
        ggtitle(title) 
}

load_data <- function(round=0) {
    a <- load_training_data(sprintf("M.SssI.ecoli_e2925.sqk006.sorted.bam.round%d.methyltrain.tsv.sampled", round), "methylated")
    b <- load_training_data(sprintf("NA12878.merged.sorted.bam.round%d.methyltrain.tsv.sampled", round), "human")
    c <- load_training_data(sprintf("loman.ecoli_k12.pcr.sqk006.sorted.bam.round%d.methyltrain.tsv.sampled", round), "pcr")
    d <- load_training_data(sprintf("loman.ecoli_k12.sqk006.sorted.bam.round%d.methyltrain.tsv.sampled", round), "k12")
    all <- rbind(a,b,c,d)
    setkey(all, model_kmer)
    return(all)
}
