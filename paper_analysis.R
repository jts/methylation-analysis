get_gaussian_parameters <- function(model, k) {
    
    row <- model[k]
    mean = row[1,]$level_mean
    stdv = row[1,]$level_stdv
    return(c(mean, stdv))
}

read_input_model <- function(filename) {
    df <- read.table(filename, header=T)
    dt <- data.table(df)
    setkey(dt, kmer)
    return(dt)
}


read_trained_model <- function(filename) {
    df <- read.table(filename, col.names=c("kmer", "level_mean", "level_stdv", "sd_mean", "sd_stdv"))
    dt <- data.table(df)
    setkey(dt, kmer)
    return(dt)
}

gaussian_kl <- function(m1, s1, m2, s2) {
    kl <- log(s2 / s1) + ((s1 ^ 2 + (m1 - m2) ^ 2) / (2 * s2 ^ 2)) - (0.5)
    return(kl)
}

get_all_distributions <- function() {

    dam <- analyze_distributions("dam");
    dcm_a <- analyze_distributions("dcm_a");
    dcm_t <- analyze_distributions("dcm_t");
    cg <- analyze_distributions("CpG");
    return(rbind(cg, dam, dcm_a, dcm_t))
}

template_name = "t.006"
complement_pop1_name = "c.p1.006"
complement_pop2_name = "c.p2.006"
models <- c(template_name, complement_pop1_name, complement_pop2_name)

generate_shifted_kmer_table <- function(data, threshold = 1.0) {

    require(plyr)

    ddply(data, "model_strand", function(x) {
        total <- nrow(x)
        total_methylated <- sum(str_count(x$kmer, "M") > 0)
        shifted <- sum(abs(x$delta) > threshold)
        shifted_methylated <- sum(abs(x$delta) > threshold & str_count(x$kmer, "M") > 0)
        data.frame(total_kmers = total,
                   num_shifted = shifted, 
                   total_methylated = total_methylated,
                   num_methylated_shifted = shifted_methylated)
    })
}

plot_shift_by_model_strand <- function(data) {
    require(ggplot2)
    ggplot(subset(data, num_methylated_sites <= 1), aes(delta, fill=site_positions)) +
        geom_histogram(aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]), binwidth=0.05, position="identity", alpha=0.5) +
        facet_grid(site_positions ~ model_strand)
}


get_df_for_training_set <- function(alphabet_type, model_key) {
    return(rbind(get_df_for_training_set_by_strand(alphabet_type, model_key, template_name),
                 get_df_for_training_set_by_strand(alphabet_type, model_key, complement_pop1_name),
                 get_df_for_training_set_by_strand(alphabet_type, model_key, complement_pop2_name)))
}

get_df_for_training_set_by_strand <- function(alphabet_type, trained_model_key, model_strand) {
    require(stringr)
    require(data.table)

    methylated_symbol = "M"
    
    methylated_pattern = ""
    unmethylated_pattern = ""
    
    # Initialize data dependent on the type of methylation we're looking at
    if(alphabet_type == "dam") {
        methylated_pattern <- "GMTC"
        unmethylated_pattern <- "GATC"
    } else if(alphabet_type == "CpG") {
        methylated_pattern <- "MG"
        unmethylated_pattern <- "CG"
    } else if(alphabet_type == "CpG_negative") {
        methylated_pattern <- "MG"
        unmethylated_pattern <- "CG"
    } else if(alphabet_type == "dcm_a") {
        methylated_pattern <- "CMAGG"
        unmethylated_pattern <- "CCAGG"
    } else if(alphabet_type == "dcm_t") {
        methylated_pattern <- "CMTGG"
        unmethylated_pattern <- "CCTGG"
    } else if(alphabet_type == "nucleotide") {
        methylated_pattern <- ""
        unmethylated_pattern <- ""
    }

    # Determine which DNA base is methylated
    methylated_base <- ""
    offset = 0
    for(i in 1:str_length(methylated_pattern)) {
        b <- str_sub(methylated_pattern, i, i)
        if(b == methylated_symbol) {
            methylated_base <- str_sub(unmethylated_pattern, i, i)
            offset <- i
        }
    }
    match_prefix = str_sub(methylated_pattern, 1, offset - 1)
    match_suffix = str_sub(methylated_pattern, offset + 1)

    # This is the model we compare levels to
    base_model_name <- sprintf("%s.ont.model", model_strand)
    input_model <- read_input_model(base_model_name)

    # This is the trained model
    trained_model <- read_trained_model(sprintf("%s.%s.trained.model", model_strand, trained_model_key))
    
    # Debug print
    print(sprintf("Methylated base: %s prefix: %s suffix: %s base model: %s", 
                  methylated_base, match_prefix, match_suffix, base_model_name))

    kmer <- vector()
    num_methylated_sites <- vector()
    site_positions <- vector()
    kl_score <- vector()
    delta <- vector()
    masked_kmers <- vector()

    # Iteraate over the kmers
    trained_kmers <- trained_model$kmer
    for(test_kmer in trained_kmers) {

        # Transform the kmer into the unmethylated version
        unmethylated_kmer <- gsub(methylated_symbol, methylated_base, test_kmer)
        
        # Look up the gaussian parameters from the input and trained model
        p1 <- get_gaussian_parameters(trained_model, test_kmer)
        p2 <- get_gaussian_parameters(input_model, unmethylated_kmer)

        m1 <- p1[1]
        s1 <- p1[2]
        
        m2 <- p2[1]
        s2 <- p2[2]
        
        # Calculate summary stats
        kl <- gaussian_kl(m1, s1, m2, s2)
        d <- m1 - m2

        num_sites <- 0
        num_valid_sites <- 0
        is_valid_kmer <- 1
        site_str <- ""
        
        # Check if this kmer contains a methylation site at an invalid location
        # This is necessary because the model contains all possible k^sigma strings
        # and some of these strings do not contain a valid recognition site
        masked_kmer = test_kmer
        for(i in 1:str_length(test_kmer)) {
            s = str_sub(test_kmer, i, i)
            if(s == methylated_symbol) {
                num_sites <- num_sites + 1
                
                # Compare the methylation site to the expected pattern
                is_valid_pattern = 1

                for(j in 1:str_length(methylated_pattern)) {
                    idx_t <- i - offset + j
                    idx_p <- j    
                    base_t <- str_sub(test_kmer, idx_t, idx_t)
                    base_p <- str_sub(methylated_pattern, idx_p, idx_p)
                    if(base_t != "" && base_t != base_p) {
                        is_valid_pattern = 0
                        is_valid_kmer = 0
                    }
                }

                if(is_valid_pattern) {
                    num_valid_sites <- num_valid_sites + 1
                    site_str <- paste(site_str, i, sep=",")
                }
            } else {
                substr(masked_kmer, i, i) <- 'a'
            }
        }

        # Also check whether the kmer contains an unmethylated
        # recognition site, which is also invalid
        if(methylated_pattern != "" && str_count(test_kmer, unmethylated_pattern) > 0) {
            is_valid_kmer = 0
        }

        # If the kmer is valid, record stats
        if(is_valid_kmer) {
            kmer <- c(kmer, test_kmer)
            num_methylated_sites <- c(num_methylated_sites, num_valid_sites)
            kl_score <- c(kl_score, kl)
            site_positions <- c(site_positions, site_str)
            delta <- c(delta, d)
            masked_kmers <- c(masked_kmers, masked_kmer)
        }
    }

    df = data.frame(kmer, num_methylated_sites, kl_score, site_positions, delta, masked_kmers)
    df$type <- alphabet_type
    df$model_strand = model_strand
    return(df)
}

analyze_accuracy <- function()
{
    # Load site comparison file
    data <- read.table("NA12878.merged.site_comparison.tsv", header=T)

    depth_threshold <- 5

    # Filter down to high-depth sites that are binary methylated/unmethylated
    binary_sites <- subset(data, num_sites == 1 & bisulfite_depth >= depth_threshold & (bisulfite_percent_methylated == 0.0 | bisulfite_percent_methylated == 100.0))
    print(dim(binary_sites))

    # Make classifications of the bisulfite and ONT data
    ont_threshold = 0.5
    binary_sites$ont_call <- ifelse(binary_sites$ont_p_methylated < ont_threshold, "unmethylated", "methylated")
    
    bisulfite_threshold = 50.0
    binary_sites$bisulfite_call <- ifelse(binary_sites$bisulfite_percent_methylated < bisulfite_threshold, "unmethylated", "methylated")

    print(head(binary_sites))
    print(table(binary_sites$bisulfite_call, binary_sites$ont_call))
    print(summary(table(binary_sites$bisulfite_call, binary_sites$ont_call)))
    return(binary_sites)
}

context_accuracy <- function(df)
{
    
    get_kmer_at_idx <- function(x, idx) {
        k <- 6
        return(str_sub(x["context"], idx, idx + k - 1))
    }

    df$k1 <- apply(df, 1, get_kmer_at_idx, idx=1)
    df$k2 <- apply(df, 1, get_kmer_at_idx, idx=2)
    df$k3 <- apply(df, 1, get_kmer_at_idx, idx=3)
    df$k4 <- apply(df, 1, get_kmer_at_idx, idx=4)
    df$k5 <- apply(df, 1, get_kmer_at_idx, idx=5)
    df$k6 <- apply(df, 1, get_kmer_at_idx, idx=6)
    return(df)
}

context_accuracy_plot <- function(d)
{
    require(plyr)
    d$correct = d$ont_call == d$bisulfite_call
    return(ddply(d, "k1", summarise, accuracy = (mean(correct))))
}
