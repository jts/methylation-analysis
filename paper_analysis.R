get_gaussian_parameters <- function(model, k) {
    
    row <- model[k]
    mean = row[1,]$level_mean
    stdv = row[1,]$level_stdv
    return(c(mean, stdv))
}

read_model <- function(filename) {
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

analyze_distributions <- function(type) {
    require(stringr)
    require(data.table)

    methylated_symbol = "M"
    
    methylated_pattern = ""
    unmethylated_pattern = ""

    trained_model_key = ""
    input_model_key = ""
    
    if(type == "dam") {
        methylated_pattern <- "GMTC"
        unmethylated_pattern <- "GATC"
        trained_model_key = "dam"
        input_model_key = "dam"
    } else if(type == "CpG") {
        methylated_pattern <- "MG"
        unmethylated_pattern <- "CG"
        trained_model_key = "M.SssI"
        input_model_key = "cpg"
    } else if(type == "dcm_a") {
        methylated_pattern <- "CMAGG"
        unmethylated_pattern <- "CCAGG"
        trained_model_key = "dcm"
        input_model_key = "dcm"
    } else if(type == "dcm_t") {
        methylated_pattern <- "CMTGG"
        unmethylated_pattern <- "CCTGG"
        trained_model_key = "dcm"
        input_model_key = "dcm"
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

    print(paste("methylated base:", methylated_base, sep=" "))

    match_prefix = str_sub(methylated_pattern, 1, offset - 1)
    match_suffix = str_sub(methylated_pattern, offset + 1)
    print(sprintf("match -- prefix: %s suffix: %s", match_prefix, match_suffix))

    input_model <- read_model(sprintf("t.006.pcr.trained.%s_expanded.model", input_model_key))
    trained_model <- read_model(sprintf("t.006.%s.trained", trained_model_key))

    trained_kmers <- trained_model$kmer

    kmer <- vector()
    num_methylated_sites <- vector()
    site_positions <- vector()
    kl_score <- vector()
    delta <- vector()

    for(test_kmer in trained_kmers) {
        unmethylated_kmer <- gsub(methylated_symbol, methylated_base, test_kmer)
        
        p1 <- get_gaussian_parameters(trained_model, test_kmer)
        p2 <- get_gaussian_parameters(input_model, unmethylated_kmer)

        m1 <- p1[1]
        s1 <- p1[2]
        
        m2 <- p2[1]
        s2 <- p2[2]
        
        kl <- gaussian_kl(m1, s1, m2, s2)
        d <- m1 - m2

        num_sites <- 0
        num_valid_sites <- 0
        is_valid_kmer <- 1
        site_str <- ""
        
        # Check if this kmer contains a methylation site at an invalid location
        # This is necessary because the model contains all possible k^sigma strings
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
            }
        }

        # Also check whether the kmer contains an unmethylated
        # recognition site, which is also invalid
        if(str_count(test_kmer, unmethylated_pattern) > 0) {
            is_valid_kmer = 0
        }

        # Is this a valid k-mer over the methylation alphabet?
        if(is_valid_kmer) {
            kmer <- c(kmer, test_kmer)
            num_methylated_sites <- c(num_methylated_sites, num_valid_sites)
            kl_score <- c(kl_score, kl)
            site_positions <- c(site_positions, site_str)
            delta <- c(delta, d)
        }
    }

    df = data.frame(kmer, num_methylated_sites, kl_score, site_positions, delta)
    df$type <- type
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
