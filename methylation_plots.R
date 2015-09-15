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
    
    in_model = "t"
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
    p <- ggplot(all, aes(level_mean, fill=dataset)) + 
            #geom_density(alpha=0.5) +
            geom_histogram(aes(y = ..density..), alpha=0.4, binwidth=0.1, position="identity") +
            stat_function(fun = dnorm, colour="red", arg=list(mean=m_mean, sd=m_stdv)) + 
            stat_function(fun = dnorm, arg=list(mean=c_mean, sd=c_stdv)) + 
            ggtitle(paste("event mean current", m_kmer))
    return(p)
}

plot_event_stdv_for_kmer <- function(m_kmer, unmethylated_data, methylated_data, params)
{
    require(plyr)
    require(ggplot2)
    require(statmod)

    # Convert Ms in the kmer to Cs for the control data
    c_kmer <- gsub("M", "C", m_kmer)

    # subset to the required kmers
    u_subset <- unmethylated_data[c_kmer]
    m_subset <- methylated_data[m_kmer]

    # get the gaussian parmeters for these kmers
    c_sd_mean = subset(params, kmer == c_kmer)[1,]$sd_meana
    c_sd_stdv = subset(params, kmer == c_kmer)[1,]$sd_stdv
    
    m_sd_mean = subset(params, kmer == m_kmer)[1,]$sd_mean
    m_sd_stdv = subset(params, kmer == m_kmer)[1,]$sd_stdv

    if(nrow(u_subset) > 0) {
        u_subset$dataset = paste(u_subset$dataset, c_kmer)
    }

    if(nrow(m_subset) > 0) {
        m_subset$dataset = paste(m_subset$dataset, m_kmer)
    }

    # merge the subsetted data into a single data frame
    all <- rbind(u_subset, m_subset)

    m_ig_lambda = m_sd_mean ** 3 / m_sd_stdv ** 2
    c_ig_lambda = c_sd_mean ** 3 / c_sd_stdv ** 2
    # plot event means
    p <- ggplot(all, aes(level_stdv, fill=dataset)) + 
            #geom_density(alpha=0.5) +
            geom_histogram(aes(y = ..density..), alpha=0.4, binwidth=0.1, position="identity") +
            stat_function(fun = dinvgauss, colour="red", arg=list(mean=m_sd_mean, shape=m_ig_lambda)) + 
            stat_function(fun = dinvgauss, arg=list(mean=c_sd_mean, shape=c_ig_lambda)) + 
            ggtitle(paste("event current stdv", m_kmer))
    return(p)
}

plot_event_duration_for_kmer <- function(m_kmer, unmethylated_data, methylated_data, params)
{
    require(plyr)
    require(ggplot2)
    require(statmod)

    # Convert Ms in the kmer to Cs for the control data
    c_kmer <- gsub("M", "C", m_kmer)

    # subset to the required kmers
    u_subset <- unmethylated_data[c_kmer]
    m_subset <- methylated_data[m_kmer]

    if(nrow(u_subset) > 0) {
        u_subset$dataset = paste(u_subset$dataset, c_kmer)
    }

    if(nrow(m_subset) > 0) {
        m_subset$dataset = paste(m_subset$dataset, m_kmer)
    }

    # merge the subsetted data into a single data frame
    all <- rbind(u_subset, m_subset)

    p <- ggplot(all, aes(duration, fill=dataset)) + 
            geom_histogram(aes(y = ..density..), alpha=0.4, binwidth=0.02, position="identity") +
            ggtitle(paste("event duration", m_kmer)) +
            xlim(0, 0.2)
    return(p)
}


generate_training_plot <- function(outfile, twomer, control_data, methylated_data, params, plot_func)
{
    pdf(outfile, 48, 16)
    kmers = make_context_mers(twomer, 3, 0)
    plots = c()
    for(i in 1:length(kmers)) {
        print(kmers[[i]])
        plots[[i]] = plot_func(kmers[[i]], control_data, methylated_data, params)
    }

    multiplot(plotlist=plots, cols=8)
    dev.off()
}

make_training_plots <- function(training_in, control_in)
{
    control_training <- load_training_data(control_in, "unmethylated")
    methylated_training <- load_training_data(training_in, "methylated")
    params <- read.table("r7.3_template_median68pA.model.methyltrain", col.names=c("kmer", "level_mean", "level_stdv", "sd_mean", "sd_stdv"))
    
    generate_training_plot("training_plots_abcMG_event_mean.pdf", "MG", control_training, methylated_training, params, plot_event_means_for_kmer)
    generate_training_plot("training_plots_abcGG_event_mean.pdf", "GG", control_training, methylated_training, params, plot_event_means_for_kmer)
    
    generate_training_plot("training_plots_abcMG_event_stdv.pdf", "MG", control_training, methylated_training, params, plot_event_stdv_for_kmer)
    generate_training_plot("training_plots_abcGG_event_stdv.pdf", "GG", control_training, methylated_training, params, plot_event_stdv_for_kmer)

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
    data$dataset = paste(dataset_name, data$complement_model, sep="-")
    return(data)
}

load_strand_classification_data <- function(filename, dataset_name)
{
    data <- read.table(filename, header=T)
    data$dataset = paste(dataset_name, data$model, sep="-")
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
    ggsave(out_file, p, width=20, height=10)

    p2 <- ggplot(all, aes(sum_ll_ratio / n_cpg, color=dataset)) + geom_density(alpha=0.5)
    ggsave("read_classification_density.pdf", width=20, height=10)
}

strand_classification_plot <- function(m_file, c_file, out_file) {
    require(ggplot2)
    m_data <- load_strand_classification_data(m_file, "methylated")
    c_data <- load_strand_classification_data(c_file, "unmethylated")
    all <- rbind(m_data, c_data)

    p <- ggplot(all, aes(n_cpg, sum_ll_ratio, color=dataset)) + geom_point()
    ggsave(out_file, p, width=20, height=10)

    p2 <- ggplot(all, aes(sum_ll_ratio / n_cpg, color=dataset)) + geom_density(alpha=0.5)
    ggsave("strand_classification_density.pdf", width=20, height=10)
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
# Event alignment summary
#
load_eventalign_summary <- function(filename) {
   data <- read.table(filename, header=T)
   data$dataset = filename
   return(data)
}

#
# Site comparison
#
site_comparison_plot <- function(infile, outfile) {
    require(ggplot2)
    data <- read.table(infile, header=T)
    pdf(outfile, 10, 10)
    p <- ggplot(data, aes(bisulfite_percent_methylated, ont_p_methylated)) + geom_hex()
    multiplot(p, cols=1)
    dev.off()
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
        make_training_plots(args[2], args[3])
    } else if(command == "site_likelihood_plots") {
        site_likelihood_plots(args[2], args[3])
    } else if(command == "read_classification_plot") {
        read_classification_plot(args[2], args[3], args[4])
    } else if(command == "strand_classification_plot") {
        strand_classification_plot(args[2], args[3], args[4])
    } else if(command == "human_cpg_island_plot") {
        human_cpg_island_plot(args[2], args[3], args[4])
    } else if(command == "site_comparison_plot") {
        site_comparison_plot(args[2], args[3])
    } else if(command == "global_methylation") {
        outfile = args[length(args)]
        global_methylation_plot(outfile, as.vector(args[c(-1, -length(args))]))
    }
}
