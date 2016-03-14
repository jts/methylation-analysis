#
# Helper to plot multiple ggplot figures on a single page
# From http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#
multiplot <- function(..., plotlist=NULL, cols) {
    require(grid)
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    numPlots = length(plots)
    # Make the panel
    plotCols = cols # Number of columns of plots
    plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols

    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
    vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
        curRow = ceiling(i/plotCols)
        curCol = (i-1) %% plotCols + 1
        print(plots[[i]], vp = vplayout(curRow, curCol ))
    }
}

# from https://gist.github.com/hadley/3880527
# i use this because functions like llply(list_dataframes, SUBSET2, model_kmer == kmer)
# don't work with built-in subset
SUBSET2 <- function(dat, ex) {
  ex <- substitute(ex)
  frames <- rev(sys.frames())
  for(f in frames) {
    call <- substitute(substitute(x, f), list(x = ex))
    ex <- eval(call)
  }
  dat[eval(ex, dat, parent.frame()),]
}

#
# Methylation Training plots
#
load_training_data <- function(filename, dataset_name)
{
    require(data.table)
    data <- fread(filename)
    data$dataset = dataset_name
    setkey(data, model_kmer)
    return(data)
}

get_params_for_kmer <- function(in_kmer, params) {
    o = subset(params, kmer == in_kmer)[1,]
    return(o)
}

plot_event_means_for_kmer <- function(m_kmer, unmethylated_data, methylated_data, params, model_short_name = "t.006")
{
    require(plyr)
    require(ggplot2)

    # Convert Ms in the kmer to Cs for the control data
    c_kmer <- gsub("M", "C", m_kmer)

    # subset to the required kmers
    u_subset <- unmethylated_data[c_kmer][model == model_short_name]
    m_subset <- methylated_data[m_kmer][model == model_short_name]

    # get the gaussian parmeters for these kmers
    c_mean = subset(params, kmer == c_kmer)[1,]$level_mean
    c_stdv = subset(params, kmer == c_kmer)[1,]$level_stdv
    
    m_mean = subset(params, kmer == m_kmer)[1,]$level_mean
    m_stdv = subset(params, kmer == m_kmer)[1,]$level_stdv

    if(nrow(u_subset) > 0) {
        u_subset$dataset = paste(u_subset$dataset, c_kmer)
    }

    if(nrow(m_subset) > 0) {
        m_subset$dataset = paste(m_subset$dataset, m_kmer)
    }

    # merge the subsetted data into a single data frame
    all <- rbind(u_subset, m_subset)

    # plot event means
    p <- ggplot(all, aes(level_mean, fill=dataset)) + 
            #geom_density(alpha=0.5) +
            geom_histogram(aes(y = ..density..), alpha=0.4, binwidth=0.1, position="identity") +
            stat_function(fun = dnorm, colour="red", arg=list(mean=m_mean, sd=m_stdv)) + 
            stat_function(fun = dnorm, arg=list(mean=c_mean, sd=c_stdv)) + 
            ggtitle(paste("event mean current", m_kmer))
    return(p)
}

generate_training_plot <- function(outfile, twomer, control_data, methylated_data, params, plot_func)
{
    kmers = make_context_mers(twomer, 4, 0)
    plots = c()
    for(i in 1:length(kmers)) {
        print(kmers[[i]])
        plots[[i]] = plot_func(kmers[[i]], control_data, methylated_data, params)
    }

    n_cols = sqrt(length(plots))
    pdf(outfile, 6 * n_cols, 2 * n_cols)
    multiplot(plotlist=plots, cols=n_cols)
    dev.off()
}

make_training_plots <- function(training_in, control_in)
{
    control_training <- load_training_data(control_in, "unmethylated")
    methylated_training <- load_training_data(training_in, "methylated")
    params <- read.table("t.006.M.SssI.trained", col.names=c("kmer", "level_mean", "level_stdv", "sd_mean", "sd_stdv"))
    
    generate_training_plot("training_plots_abcMG_event_mean.pdf", "MG", control_training, methylated_training, params, plot_event_means_for_kmer)
    generate_training_plot("training_plots_abcGG_event_mean.pdf", "GG", control_training, methylated_training, params, plot_event_means_for_kmer)
    
    #generate_training_plot("training_plots_abcMG_event_stdv.pdf", "MG", control_training, methylated_training, params, plot_event_stdv_for_kmer)
    #generate_training_plot("training_plots_abcGG_event_stdv.pdf", "GG", control_training, methylated_training, params, plot_event_stdv_for_kmer)

    pdf("all_m_kmers.pdf")
    kmers <- make_context_mers("", 6, 0, c('A', 'C', 'G', 'M', 'T'))
    for(i in 1:length(kmers)) {
        curr <- kmers[[i]]
        if(length( grep("M", curr) ) > 0) {
            if(nrow(methylated_training[curr]) > 100) {
                print(curr)
                print(nrow(methylated_training[ curr ]))
                p <- plot_event_means_for_kmer(curr, control_training, methylated_training, params)
                multiplot(p, cols=1)
            }
        }
    }
    dev.off()
}

#
# Human analysis plots
#
read_ont_scores_file <- function(filename) { 
    return(read.table(filename, header=T))
}

read_bisulfite_scores_file <- function(filename) {
    return(read.table(filename, header=T))
}

human_cpg_island_plot <- function(bisulfite_file, nanopore_file, out_file) {
    require(ggplot2)
    require(stringr)

    ont <- read_ont_scores_file(nanopore_file)
    ont$ont_percent_methylated = 100 * (ont$called_sites_methylated / ont$called_sites)

    bisulfite <- read_bisulfite_scores_file(bisulfite_file)
    
    # merge the data sets together on the common Cpg island key
    merged <- merge(ont, bisulfite, by.x="key", by.y="key")
    merged$near_gene = merged$gene.x != "."
    ggplot(merged, aes(bisulfite_percent_methylated, ont_percent_methylated, color=near_gene)) + 
        geom_point() +
        xlab("ENCODE NA12878 percent methylated (bisulfite)") +
        ylab("Nanopore called sites percent methylated (ont)") +
        ggtitle("Methylation signal at CpG Islands")
   ggsave(out_file, width=10,height=10, useDingbats = FALSE)

   # plot histograms
   p1 <- ggplot(ont, aes(ont_percent_methylated, fill=gene != ".")) + geom_histogram(alpha=0.5, position="identity", binwidth=4)
   p2 <- ggplot(bisulfite, aes(bisulfite_percent_methylated, fill=gene != ".")) + geom_histogram(alpha=0.5, position="identity", binwidth=4)

   pdf("histogram.pdf")
   multiplot(plotlist=list(p1, p2), cols=1)
   dev.off()

   plot_ont_hist <- function(data, title) { 
        p <- ggplot(data, aes(ont_percent_methylated, fill=gene != ".")) +
             geom_histogram(alpha=0.5, position="identity", binwidth=4) +
             ggtitle(title)
        return(p)
   }

   plot_bs_hist <- function(data, title) {
      p <- ggplot(data, aes(bisulfite_percent_methylated, fill=gene != ".")) +
           geom_histogram(alpha=0.5, position="identity", binwidth=4) +
           ggtitle(title)
      return(p)
   }

   p_ont_x <- plot_ont_hist(subset(ont, str_sub(key, 1, 4) == "chrX"), "ONT ChrX") 
   p_ont_not_x <- plot_ont_hist(subset(ont, str_sub(key, 1, 4) != "chrX"), "ONT Not ChrX")
   p_bs_x <- plot_bs_hist(subset(bisulfite, str_sub(key, 1, 4) == "chrX"), "Bisulfite ChrX") 
   p_bs_not_x <- plot_bs_hist(subset(bisulfite, str_sub(key, 1, 4) != "chrX"), "Bisulfite Not ChrX")

   pdf("histogram_chromosomes.pdf", 15, 10)
   multiplot(plotlist=list(p_ont_x, p_ont_not_x, p_bs_x, p_bs_not_x), cols=2)
   dev.off()
}

#
# Distance to TSS analysis
#
load_distance_data <- function(filename, tag) {
    data <- read.table(filename, header=T)
    data$dataset = tag
    return(data)
}

distance_to_TSS_plot <- function(out_file) {
    require(ggplot2)

    data1 = load_distance_data("NA12878.pcr.simpson.021616.methylated_sites.distance_to_TSS.table", "ont.pcr")
    data2 = load_distance_data("NA12878.pcr_MSssI.simpson.021016.methylated_sites.distance_to_TSS.table", "ont.pcr_MSssI")
    data3 = load_distance_data("NA12878.native.merged.methylated_sites.distance_to_TSS.table", "ont.native")
    data4 = load_distance_data("bisulfite.distance_to_TSS.table", "bisulfite")
    data_all <- rbind(data1, data2, data3, data4)
    pdf(out_file, 12, 4)
    p <- ggplot(data_all, aes(max_distance, percent_methylated, group=dataset, color=dataset)) + 
            geom_point(size=1) + 
            geom_line(size=0.1) + 
            xlab("Binned distance to TSS") + 
            ylab("Percent Methylated") + 
            ylim(0, 1)
    multiplot(p, cols=1)
    dev.off()
}

#
#
#
call_accuracy_by_threshold <- function() {
    require(ggplot2)
    data <- read.table("accuracy.bythreshold.tsv", header=T)
    
    pdf("accuracy_by_threshold.pdf", 12, 6)
    p1 <- ggplot(data, aes(threshold, accuracy)) + geom_line() + xlim(0, 10) + ylim(0, 1) + xlab("Log-Likelihood ratio threshold"); 
    p2 <- ggplot(data, aes(threshold, called)) + geom_line() + xlim(0, 10) + xlab("Log-Likelihood ratio threshold"); 
    multiplot(p1, p2, cols=2); 
    dev.off()
}

#
# Global methylation analysis
#
global_methylation_plot <- function(outfile, ...) {
    require(ggplot2)
    require(tools)
    require(stringr)

    files = list(...)
    
    pl <- list()
    for(i in 1:length(files[[1]])) {
        
        f <- files[[1]][[i]]
        suf_pos = str_locate(f, ".sorted.bam.methyltest.sites.tsv")

        base = str_sub(f, 0, suf_pos[,1] - 1)
        
        d <- load_site_data(f, base)

        # Calculate posterior probabilty the site is methylated
        d$posterior = 1 / (1 + exp(-d$LL_RATIO))

        # Get genome-wide average methylation
        avg = mean(d$posterior)
        
        # Reset dataset name to include average
        name <- sprintf("%s (avg 5mC: %.2f)", base, avg)
        d$dataset = name

        pl[[i]] <- ggplot(subset(d, N_CPG == 1), aes(posterior)) + geom_histogram(binwidth=0.05) + xlab("P(methylated | D)") + ggtitle(name)
    }
    
    w = 4
    h = length(pl) * w
    pdf(outfile, height=h, width=w)
    multiplot(plotlist=pl, cols=1)
    dev.off()
}

#
# Utility functions for generating kmers
#
make_context_mers <- function(context, n_bases_before, n_bases_after, alphabet=c('A', 'C', 'G', 'T')) {

    prefix = make_mers(n_bases_before, alphabet)
    suffix = make_mers(n_bases_after, alphabet)
    out = c()
    for(p in prefix) {
        for(s in suffix) {
            out = c(out, paste(p, context, s, sep=''))
        }
    }
    return(out)
}

add_bases <- function(strings, alphabet) {
    out = c()
    for(s in strings) {
        for(b in alphabet) {
            out = c(out, paste(s,b, sep='')) # ya i know
        }
    }
    return(out)
}

make_mers <- function(n, alphabet) {
    l = c("")
    if(n == 0) {
        return(l)   
    }
    
    for(i in 1:n) {
        l = add_bases(l, alphabet)
    }
    return(l)
}

#
# Main, when called from Rscript
#
args <- commandArgs(TRUE)
command = args[1]

if(! interactive()) {
    if(command == "training_plots") {
        make_training_plots(args[2], args[3])
    } else if(command == "strand_classification_plot") {
        strand_classification_plot(args[2], args[3], args[4])
    } else if(command == "human_cpg_island_plot") {
        human_cpg_island_plot(args[2], args[3], args[4])
    } else if(command == "global_methylation") {
        outfile = args[length(args)]
        global_methylation_plot(outfile, as.vector(args[c(-1, -length(args))]))
    }
}
