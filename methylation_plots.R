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
    data <- read.table(filename, header=T)
    data$dataset = dataset_name
    
    in_model = "r7.3_template_median68pA.model"
    tbl = data.table(subset(data, model == in_model))
    
    setkey(tbl, model_kmer)
    return(tbl)
}

plot_event_means_for_kmer <- function(m_kmer, unmethylated_data, methylated_data, params)
{
    require(plyr)
    require(ggplot2)

    # Convert Ms in the kmer to Cs for the control data
    c_kmer <- gsub("M", "C", m_kmer)

    # subset to the required kmers
    u_subset <- unmethylated_data[c_kmer]
    m_subset <- methylated_data[m_kmer]

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
    p <- ggplot(all, aes(event_mean, fill=dataset)) + 
            #geom_density(alpha=0.5) +
            geom_histogram(aes(y = ..density..), alpha=0.4, binwidth=0.1, position="identity") +
            stat_function(fun = dnorm, colour="red", arg=list(mean=m_mean, sd=m_stdv)) + 
            stat_function(fun = dnorm, arg=list(mean=c_mean, sd=c_stdv)) + 
            ggtitle(paste("event mean current", m_kmer))
    return(p)
}

make_training_plots <- function()
{
    control_training <- load_training_data("ERX708228.ecoli.sorted.bam.methyltrain.0.tsv", "unmethylated")
    methylated_training <- load_training_data("M.SssI.e2925_ecoli.sorted.bam.methyltrain.9.tsv", "methylated")
    params <- read.table("r7.3_template_median68pA.model.methyltrain", col.names=c("kmer", "level_mean", "level_stdv", "sd_mean", "sd_stdv"))
    
    pdf("training_plots.pdf", 48, 16)
    
    # CG in the last position
    kmers = make_context_mers("MG", 3, 0)
    plots = c()
    for(i in 1:length(kmers)) {
        print(kmers[[i]])
        plots[[i]] = plot_event_means_for_kmer(kmers[[i]], control_training, methylated_training, params)
    }

    multiplot(plotlist=plots, cols=8)
    dev.off()
    
    pdf("control_training_plots.pdf", 48, 16)
    
    # GG in the last position as a control
    kmers = make_context_mers("GG", 3, 0)
    plots = c()
    for(i in 1:length(kmers)) {
        plots[[i]] = plot_event_means_for_kmer(kmers[[i]], control_training, methylated_training, params)
    }

    multiplot(plotlist=plots, cols=8)
    dev.off()

    pdf("training_all_m_mers.pdf")
    #kmers = make_context_mers("MG", 3, 0)
    kmers <- make_context_mers("", 5, 0, c('A', 'C', 'G', 'M', 'T'))
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
# Methylation Testing plots
#

load_site_data <- function(filename, dataset_name)
{
    data <- read.table(filename, header=T)
    data$dataset = dataset_name
    return(data)
}

load_read_classification_data <- function(filename, dataset_name)
{
    data <- read.table(filename, header=T)
    data$dataset = dataset_name
    return(data)
}

#
# Plot the histogram of scores from a sites.tsv
#
site_histogram <- function(in_kmer, dataset) {
    require(plyr)
    require(ggplot2)
    
    sub <- subset(dataset, SEQUENCE == in_kmer)
    
    # hack, give at least one row
    sub <- rbind(sub, dataset[1,])

    p1 <- ggplot(sub, aes(LL_RATIO)) + 
            geom_histogram(aes(y = ..density..), alpha=0.4, binwidth=0.1, position="identity") +
            xlim(-10, 10) +
            ggtitle(in_kmer)

    return(p1)
}

site_likelihood_plots <- function(in_file, out_file) {
    data <- load_site_data(in_file, "methylated")
    pdf(out_file, 32, 16)
    
    # CG in the last position
    kmers = make_context_mers("CG", 3, 0)
    plots = c()
    for(i in 1:length(kmers)) {
        plots[[i]] = site_histogram(kmers[[i]], data)
    }
    multiplot(plotlist=plots, cols=8)
    dev.off()
}

read_classification_plot <- function(m_file, c_file, out_file) {
    require(ggplot2)
    m_data <- load_read_classification_data(m_file, "methylated")
    c_data <- load_read_classification_data(c_file, "unmethylated")
    all <- rbind(m_data, c_data)

    p <- ggplot(all, aes(n_cpg, sum_ll_ratio, color=dataset)) + geom_point()

    ggsave(out_file, p)
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

    ont <- read_ont_scores_file(nanopore_file)
    bisulfite <- read_bisulfite_scores_file(bisulfite_file)
    # merge the data sets together on the common Cpg island key
    merged <- merge(ont, bisulfite, by.x="key", by.y="key")
    merged$near_gene = merged$gene.x != "."
    ggplot(merged, aes(bisulfite_percent_methylated, 1 / (1 + exp(-sum_ll_ratio / n_ont_sites)), color=near_gene)) + 
        geom_point() +
        xlab("ENCODE NA12878 percent methylated (bisulfite)") +
        ylab("P(methylated | ONT)") +
        ggtitle("Methylation signal at CpG Islands")


   ggsave(out_file, width=10,height=10)

   # plot histograms
   p1 <- ggplot(ont, aes(100 * (1 / (1 + exp(-sum_ll_ratio / n_ont_sites))), fill=gene != ".")) + geom_histogram(alpha=0.5, position="identity", binwidth=4)
   #p1 <- ggplot(ont, aes(sum_posterior / n_ont_sites, fill=gene != ".")) + geom_density(alpha=0.5)
   p2 <- ggplot(bisulfite, aes(bisulfite_percent_methylated, fill=gene != ".")) + geom_histogram(alpha=0.5, position="identity", binwidth=4)
   pdf("histogram.pdf")
   multiplot(plotlist=list(p1, p2), cols=1)
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
    
    all <- NULL
    for(f in files[[1]]) {
        
        suf_pos = str_locate(f, ".sorted.bam.methyltest.sites.tsv")

        base = str_sub(f, 0, suf_pos[,1] - 1)
        
        d <- load_site_data(f, base)

        # Calculate posterior probabilty the site is methylated
        d$posterior = 1 / (1 + exp(-d$LL_RATIO))

        # Get genome-wide average methylation
        avg = mean(d$posterior)
        
        # Reset dataset name to include average
        d$dataset = sprintf("%s (avg 5mC: %.2f)", base, avg)

        if(is.null(all)) {
            all = d
        } else {
            all = rbind(all, d)
        }
    }

    print(head(all))
    #ggplot(all, aes((1 / (1 + exp(-LL_RATIO))), fill=dataset)) + geom_histogram(position="identity", binwidth = 0.04, alpha=0.5)
    ggplot(all, aes(posterior, fill=dataset)) + geom_density(alpha=0.5) + xlab("P(methylated | D)")
    ggsave(outfile)
}

# Plot multiple 5-mers on a single plot
site_multiplot <- function(kmers, ...) {
    plots = c()
    for(i in 1:length(kmers)) {
        plots[[i]] = site_histogram(kmers[[i]], ...)
    }
    multiplot(plotlist=plots, cols=8)
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
        make_training_plots()
    } else if(command == "site_likelihood_plots") {
        site_likelihood_plots(args[2], args[3])
    } else if(command == "read_classification_plot") {
        read_classification_plot(args[2], args[3], args[4])
    } else if(command == "human_cpg_island_plot") {
        human_cpg_island_plot(args[2], args[3], args[4])
    } else if(command == "global_methylation") {
        outfile = args[length(args)]
        global_methylation_plot(outfile, as.vector(args[c(-1, -length(args))]))
    }
}
