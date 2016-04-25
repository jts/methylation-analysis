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
            ggtitle(paste("event mean current", m_kmer)) + 
            global_theme()

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

make_panel <- function(names, data_sets, param_means, param_stdvs, panel_label)
{
    require(grid)
    palette <- c("red", "blue")

    # merge datasets
    all_data <- NULL
    for(n in names) {
        if(is.null(all_data)) {
            all_data <- data_sets[[n]]
        } else {
            all_data <- rbind(all_data, data_sets[[n]])
        }
    }

    # extract the kmer we are plotting for the title
    # NB all rows have the same kmer
    kmer <- head(all_data, 1)$model_kmer

    # make the histogram
    p <- ggplot(all_data, aes(level_mean, fill=dataset)) +
            geom_histogram(aes(y = ..density..), alpha=0.4, binwidth=0.1, position="identity") +
            scale_fill_manual(values=palette) +
            xlab("Measured event level (pA)") +
            ggtitle(paste("Event distribution for k-mer", kmer))

    # add panel label
    p <- p + annotation_custom(textGrob(label = panel_label, x = 0.05, y = 0.95, gp=gpar(fontsize=20)))

    # add gaussian fits
    idx <- 1
    for(n in names) {
        p <- p + stat_function(fun = dnorm, colour=palette[[idx]], args=list(mean=param_means[[n]], sd=param_stdvs[[n]]))
        idx <- idx + 1
    }

    return(p + global_theme())
}

#
# Emissions figure - examples of distributions
#
make_emissions_figure <- function(outfile,
                                  panelA_file,
                                  panelB1_file,
                                  panelB2_file,
                                  panelC1_file,
                                  panelC2_file)
{
    require(stringr)

    data_sets <- list()
    param_means <- list()
    param_stdvs <- list()

    idx <- 1
    for(current_file in c(panelA_file, panelB1_file, panelB2_file, panelC1_file, panelC2_file)) {

        # Parse the structured name
        name_fields = str_split(current_file, fixed("."))[[1]]
        display_name <- str_replace(str_replace(name_fields[3], "pcr", "PCR"), "_", "+")

        # Load the data
        data_sets[[current_file]] <- load_training_data(current_file, display_name)

        # parse model parameters
        strand = head(data_sets[[current_file]], 1)$model
        kmer = head(data_sets[[current_file]], 1)$model_kmer
        modelname  <- str_c(c(strand, name_fields[2:6], "model"), collapse=".")
        params <- read.table(modelname, col.names=c("kmer", "level_mean", "level_stdv", "sd_mean", "sd_stdv"))

        # save the parameters for this kmer
        kmer_params <- get_params_for_kmer(kmer, params)
        param_means[[current_file]] <- kmer_params[1,]$level_mean
        param_stdvs[[current_file]] <- kmer_params[1,]$level_stdv
    }

    pdf(outfile, 7, 15)
    panel_A <- make_panel(list(panelA_file), data_sets, param_means, param_stdvs, "A")
    panel_B <- make_panel(list(panelB1_file, panelB2_file), data_sets, param_means, param_stdvs, "B")
    panel_C <- make_panel(list(panelC1_file, panelC2_file), data_sets, param_means, param_stdvs, "C")

    multiplot(panel_A, panel_B, panel_C, cols=1)
    dev.off()
}

#
# Make a plot of the k-mer mean differences by position of the methylated base
#
make_mean_shift_by_position_figure <- function(outfile, filename) {
    data <- read.table(filename, header=T)
    data$model = factor(data$model, levels = c("template", "comp.pop1", "comp.pop2"))

    pdf(outfile)
    p <- ggplot(data, aes(difference)) + 
             geom_histogram(binwidth=0.25) + 
             facet_grid(m_pattern ~ model, scales="free_y") + 
             xlim(-4, 4) + 
             global_theme();
    multiplot(p, cols=1)
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

# from http://www.sthda.com/english/wiki/ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page-r-software-and-data-visualization
get_legend<-function(myggplot){
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}

human_cpg_island_plot <- function(bisulfite_file, nanopore_file, out_file) {
    require(ggplot2)
    require(gridExtra)
    require(stringr)

    # Load and remove zero-depth sites
    ont <- subset(read_ont_scores_file(nanopore_file), called_sites > 0)
    ont$ont_percent_methylated = 100 * (ont$called_sites_methylated / ont$called_sites)

    bisulfite <- subset(read_bisulfite_scores_file(bisulfite_file), called_sites > 0)
    bisulfite$bisulfite_percent_methylated = 100 * (bisulfite$called_sites_methylated / bisulfite$called_sites)
    
    # merge the data sets together on the common Cpg island key
    merged <- merge(ont, bisulfite, by.x="key", by.y="key")
    merged$in_promoter = merged$feature.x == "promoter"
    scatter <- ggplot(merged, aes(bisulfite_percent_methylated, ont_percent_methylated, colour=in_promoter)) + 
        geom_point() +
        xlab("Percent methylated (bisulfite)") +
        ylab("Percent methylated (nanopore)") +
        scale_colour_discrete(name="Is CGI in a promoter?") + 
        global_theme() + xlim(0,100) + ylim(0,100)

    legend <- get_legend(scatter)

    scatter <- scatter + theme(legend.position="none")

    hist_top <- ggplot(merged,aes(bisulfite_percent_methylated, fill=in_promoter)) +
                geom_histogram(binwidth=1) + 
                global_theme() + theme(legend.position="none") +
                theme(axis.title.x = element_blank()) 

    hist_right <- ggplot(merged, aes(ont_percent_methylated, fill=in_promoter)) + 
                geom_histogram(binwidth=1,aes(fill=in_promoter)) + 
                coord_flip() +
                global_theme() + theme(legend.position="none") + theme(axis.title.y = element_blank()) 

    empty <- ggplot()+geom_point(aes(1,1), colour="white")+
         theme(axis.ticks=element_blank(), 
               panel.background=element_blank(), 
               axis.text.x=element_blank(), axis.text.y=element_blank(),           
               axis.title.x=element_blank(), axis.title.y=element_blank())

    pdf(out_file, width=10, height=10)
    grid.arrange(hist_top, legend, scatter, hist_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
    dev.off()

    # Write the correlation to a file
    cor_str <- paste("Correlation", cor(merged$bisulfite_percent_methylated, merged$ont_percent_methylated))
    fp <- file(paste(out_file, "correlation", sep="."))
    writeLines(cor_str, fp)
    close(fp)
}

#
# Distance to TSS analysis
#
make_display_name <- function(structured_name) {
    require(stringr)
    fields <- str_split(structured_name, fixed("."))[[1]]

    sample = fields[1]
    treatment = fields[2]

    # Rename treatment
    if(treatment == "pcr_MSssI") {
        treatment = "PCR+M.SssI"
    } else if(treatment == "pcr") {
        treatment = "PCR"
    } else if(treatment == "native") {
        treatment = "natural"
    }

    id = str_c(sample, treatment, sep=" ")

    if(treatment == "bisulfite") {
        return(id)
    } else {
        
        lab = fields[3]
        date = fields[4]
        if(lab == "merged") {
            return(sprintf("%s (%s)", id, lab))
        } else {
            return(sprintf("%s (%s-%s)", id, lab, date))
        }
    }
}

load_distance_data <- function(filename, tag) {
    data <- read.table(filename, header=T)
    data$dataset = make_display_name(filename)
    return(data)
}

load_all_TSS_data <- function(...) {
    in_files = list(...)
    all_data = NULL
    for(i in 1:length(in_files[[1]])) {
        filename = in_files[[1]][[i]]
        data = load_distance_data(filename, filename)
        if(is.null(all_data)) {
            all_data = data
        } else {
            all_data = rbind(all_data, data)
        }
    }
    return(all_data)
}

TSS_distance_plot <- function(out_file, ...) {
    require(ggplot2)
 
    all_data = load_all_TSS_data(...)   

    pdf(out_file, 12, 4)
    p <- ggplot(subset(all_data, chromosome == "autosomes"), aes(distance, percent_methylated, group=dataset, color=dataset)) +
            geom_point(size=1) + 
            geom_line(size=0.1) + 
            xlab("Binned distance to TSS") + 
            ylab("Percent methylated") + 
            ylim(0, 100) +
            global_theme()
    multiplot(p, cols=1)
    dev.off()
}

TSS_distance_plot_by_chromosome <- function(out_file, ...) {
    require(ggplot2)
    
    all_data = load_all_TSS_data(...)   

    #http://stackoverflow.com/questions/19014531/sort-by-chromosome-name
    chrOrder <-paste("chr", c((1:22),"X","Y","M"), sep="")
    all_data$chromosome = factor(all_data$chromosome, levels=chrOrder, ordered=TRUE)
    pdf(out_file, 10, 60)
    p <- ggplot(subset(all_data, chromosome != "all"), aes(distance, percent_methylated, group=dataset, color=dataset)) + 
            geom_point(size=1) + 
            geom_line(size=0.1) + 
            xlab("Binned distance to TSS") + 
            ylab("Percent methylated") + 
            ylim(0, 100) +
            facet_grid(chromosome ~ .) +
            theme_bw()
    multiplot(p, cols=1)
    dev.off()
}

#
#
#
call_accuracy_by_threshold <- function(in_file, out_file) {
    require(ggplot2)
    require(grid)
    data <- read.table(in_file, header=T)
    
    pdf(out_file, 12, 6)
    p1 <- ggplot(data, aes(threshold, 1 - accuracy)) + geom_line() + xlim(0, 10) + ylim(0, 0.20) + xlab("Log likelihood ratio threshold") + ylab("Error rate") + global_theme()
    p2 <- ggplot(data, aes(threshold, called)) + geom_line() + xlim(0, 10) + xlab("Log likelihood ratio threshold") + ylab("Number of calls") + global_theme()
    
    # Add panel labels
    p1 <- p1 + annotation_custom(textGrob(label = "A", x = 0.10, y = 0.95, gp=gpar(fontsize=20)))
    p2 <- p2 + annotation_custom(textGrob(label = "B", x = 0.10, y = 0.95, gp=gpar(fontsize=20)))

    multiplot(p1, p2, cols=2); 
    dev.off()
}

call_accuracy_roc <- function(in_file, out_file) {
    require(ggplot2)
    data <- read.table(in_file, header=T)
    
    pdf(out_file, 6, 6)
    p1 <- ggplot(data, aes(1 - specificity, recall)) + 
        geom_line() + 
        xlim(0, 1) + 
        ylim(0, 1) + 
        xlab("False positive rate") + 
        ylab("True positive rate") + 
        global_theme()
    multiplot(p1, cols=1); 
    dev.off()
}

call_accuracy_by_kmer <- function(in_file, out_file) {
    require(ggplot2)
    data <- read.table(in_file, header=T)
    
    pdf(out_file, 12, 6)
    p1 <- ggplot(data, aes(kmer, 1 - accuracy)) + 
        geom_point() + 
        global_theme() +
        xlab("k-mer preceding CpG site") +
        ylab("Error rate") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    multiplot(p1, cols=1); 
    dev.off()
}

site_likelihood_distribution <- function(out_file, ...) {
    in_files = list(...)
    all_data = NULL
    for(i in 1:length(in_files[[1]])) {
        filename = in_files[[1]][[i]]
        data = read.table(filename, header=T)
        data$dataset = make_display_name(filename)
        if(is.null(all_data)) {
            all_data = data
        } else {
            all_data = rbind(all_data, data)
        }
    }

    pdf(out_file, 9, 12)
    p1 <- ggplot(subset(all_data, NumCpGs == 1), aes(LogLikRatio)) + 
        facet_grid(dataset ~ .) + 
        geom_histogram(binwidth=0.5, aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) +
        xlim(-15, 15) +
        xlab("Log likelihood ratio") +
        ylab("Density") +
        global_theme() 

    multiplot(p1, cols=1); 
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
require(ggplot2)

# set the theme all plots will use
global_theme = theme_bw

args <- commandArgs(TRUE)
command = args[1]

if(! interactive()) {
    if(command == "training_plots") {
        make_training_plots(args[2], args[3])
    } else if(command == "make_emissions_figure") {
        make_emissions_figure(args[2], args[3], args[4], args[5], args[6], args[7])
    } else if(command == "make_mean_shift_by_position_figure") {
        make_mean_shift_by_position_figure(args[2], args[3])
    } else if(command == "human_cpg_island_plot") {
        human_cpg_island_plot(args[2], args[3], args[4])
    } else if(command == "TSS_distance_plot") {
        outfile = args[length(args)]
        TSS_distance_plot(outfile, as.vector(args[c(-1, -length(args))]))
    } else if(command == "TSS_distance_plot_by_chromosome") {
        outfile = args[length(args)]
        TSS_distance_plot_by_chromosome(outfile, as.vector(args[c(-1, -length(args))]))
    } else if(command == "call_accuracy_by_threshold") {
        call_accuracy_by_threshold(args[2], args[3])
    } else if(command == "call_accuracy_roc") {
        call_accuracy_roc(args[2], args[3])
    } else if(command == "call_accuracy_by_kmer") {
        call_accuracy_by_kmer(args[2], args[3])
    } else if(command == "site_likelihood_distribution") {
        outfile = args[length(args)]
        site_likelihood_distribution(outfile, as.vector(args[c(-1, -length(args))]))
    }
}
