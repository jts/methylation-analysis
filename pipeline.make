SHELL=/bin/bash -o pipefail

##################################################
#
# Step 0. Preamble: set up paths, variables and
#         install software.
#
##################################################

# do not leave failed files around
.DELETE_ON_ERROR:

# do not delete intermediate files
.SECONDARY:

# Parameters to control execution
THREADS=4

# Path to this file
SCRIPT_DIR:=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

# Install samtools
samtools.version:
	git clone --recursive https://github.com/samtools/htslib.git
	cd htslib; git checkout 1.2.1; make
	git clone --recursive https://github.com/samtools/samtools.git
	cd samtools; git checkout 1.2; make
	-cd samtools; git log | head -1 > ../$@

# Install bwa
bwa.version:
	git clone https://github.com/lh3/bwa.git
	cd bwa; make
	-cd bwa; git log | head -1 > ../$@

# Install python libs
pythonlibs.version:
	pip install biopython >> $@
	git clone https://github.com/arq5x/poretools.git
	pip install -r poretools/requirements.txt
	cd poretools && python setup.py install
	pip freeze >> $@

# Install bedtools
bedtools.version:
	git clone https://github.com/arq5x/bedtools2.git bedtools
	cd bedtools; make
	-cd bedtools; git log | head -1 > ../$@

# Install nanopolish, automatically downloading libhdf5
nanopolish.version:
	git clone --recursive https://github.com/jts/nanopolish.git
	cd nanopolish; git checkout aba1b32; make
	-cd nanopolish; git log | head -1 > ../$@

#
PORETOOLS=poretools

#
# Output targets
#
.DEFAULT_GOAL := all
all: all-nucleotide all-cpg

all-software: pythonlibs.version bedtools.version nanopolish.version samtools.version bwa.version

all-nucleotide: ecoli_er2925.native.timp.102615.alphabet_nucleotide.fofn \
                ecoli_er2925.native.timp.110915.alphabet_nucleotide.fofn \
                ecoli_er2925.pcr.timp.113015.alphabet_nucleotide.fofn \
                ecoli_er2925.pcr_MSssI.timp.113015.alphabet_nucleotide.fofn  \
                ecoli_k12.pcr.loman.250915.alphabet_nucleotide.fofn \
                ecoli_er2925.pcr.timp.021216.alphabet_nucleotide.fofn \
                ecoli_er2925.pcr_MSssI.timp.021216.alphabet_nucleotide.fofn

all-cpg: ecoli_er2925.pcr.timp.113015.alphabet_cpg.fofn \
         ecoli_er2925.pcr_MSssI.timp.113015.alphabet_cpg.fofn \
         ecoli_k12.pcr.loman.250915.alphabet_cpg.fofn \
         ecoli_er2925.pcr.timp.021216.alphabet_cpg.fofn \
         ecoli_er2925.pcr_MSssI.timp.021216.alphabet_cpg.fofn

all-island-plots: NA12878.native.timp.093015.cpg_island_plot.pdf \
                  NA12878.native.simpson.101515.cpg_island_plot.pdf \
                  NA12878.native.simpson.103015.cpg_island_plot.pdf \
                  NA12878.native.merged.cpg_island_plot.pdf \
                  NA12878.pcr.simpson.021616.cpg_island_plot.pdf \
                  NA12878.pcr_MSssI.simpson.021016.cpg_island_plot.pdf

all-accuracy-plots: accuracy.roc.pdf accuracy.by_threshold.pdf accuracy.by_kmer.pdf site.likelihood.distribution.pdf

all-TSS-plots: methylation_by_TSS_distance.pdf methylation_by_TSS_distance_by_chromosome.pdf 

all-plots: all-accuracy-plots all-island-plots all-TSS-plots

##################################################
#
# Step 1. Prepare input data 
#
##################################################

DATA_ROOT=../data

#
# Define variables for each data set
#
ECOLI_K12_NATIVE_DATA=ecoli_k12.native.loman.run1.fasta
ECOLI_K12_PCR_RUN1_DATA=ecoli_k12.pcr.loman.250915.fasta

ECOLI_ER2925_NATIVE_RUN1_DATA=ecoli_er2925.native.timp.102615.fasta
ECOLI_ER2925_NATIVE_RUN2_DATA=ecoli_er2925.native.timp.110915.fasta

ECOLI_ER2925_MSSSI_RUN1_DATA=ecoli_er2925.MSssI.timp.100215.fasta
ECOLI_ER2925_MSSSI_RUN2_DATA=ecoli_er2925.MSssI.timp.100615.fasta

ECOLI_ER2925_PCR_RUN1_DATA=ecoli_er2925.pcr.timp.113015.fasta
ECOLI_ER2925_PCR_RUN2_DATA=ecoli_er2925.pcr.timp.021216.fasta
ECOLI_ER2925_PCR_MSSSI_RUN1_DATA=ecoli_er2925.pcr_MSssI.timp.113015.fasta
ECOLI_ER2925_PCR_MSSSI_RUN2_DATA=ecoli_er2925.pcr_MSssI.timp.021216.fasta

HUMAN_NA12878_NATIVE_RUN1_DATA=NA12878.native.timp.093015.fasta
HUMAN_NA12878_NATIVE_RUN2_DATA=NA12878.native.simpson.101515.fasta
HUMAN_NA12878_NATIVE_RUN3_DATA=NA12878.native.simpson.103015.fasta
HUMAN_NA12878_NATIVE_MERGED_DATA=NA12878.native.merged.fasta

HUMAN_NA12878_PCR_DATA=NA12878.pcr.simpson.021616.fasta
HUMAN_NA12878_PCR_MSSSI_DATA=NA12878.pcr_MSssI.simpson.021016.fasta

# For each data set that we use, define a variable containing its reference

# E.coli
ECOLI_REFERENCE=ecoli_k12.fasta
$(ECOLI_K12_NATIVE_DATA)_REFERENCE=$(ECOLI_REFERENCE)
$(ECOLI_K12_PCR_RUN1_DATA)_REFERENCE=$(ECOLI_REFERENCE)

$(ECOLI_ER2925_NATIVE_RUN1_DATA)_REFERENCE=$(ECOLI_REFERENCE)
$(ECOLI_ER2925_NATIVE_RUN2_DATA)_REFERENCE=$(ECOLI_REFERENCE)

$(ECOLI_ER2925_MSSSI_DATA)_REFERENCE=$(ECOLI_REFERENCE)
$(ECOLI_ER2925_PCR_RUN1_DATA)_REFERENCE=$(ECOLI_REFERENCE)
$(ECOLI_ER2925_PCR_RUN2_DATA)_REFERENCE=$(ECOLI_REFERENCE)
$(ECOLI_ER2925_PCR_MSSSI_RUN1_DATA)_REFERENCE=$(ECOLI_REFERENCE)
$(ECOLI_ER2925_PCR_MSSSI_RUN2_DATA)_REFERENCE=$(ECOLI_REFERENCE)

# Human
HUMAN_REFERENCE=GRCh38.primary_assembly.genome.fa
$(HUMAN_NA12878_NATIVE_RUN1_DATA)_REFERENCE=$(HUMAN_REFERENCE)
$(HUMAN_NA12878_NATIVE_RUN2_DATA)_REFERENCE=$(HUMAN_REFERENCE)
$(HUMAN_NA12878_NATIVE_RUN3_DATA)_REFERENCE=$(HUMAN_REFERENCE)
$(HUMAN_NA12878_NATIVE_MERGED_DATA)_REFERENCE=$(HUMAN_REFERENCE)
$(HUMAN_NA12878_PCR_DATA)_REFERENCE=$(HUMAN_REFERENCE)
$(HUMAN_NA12878_PCR_MSSSI_DATA)_REFERENCE=$(HUMAN_REFERENCE)

# Reference and region file used for model training
TRAINING_REFERENCE=$(ECOLI_REFERENCE)
TRAINING_REGION="Chromosome:50000-3250000"

# for debugging
METHYLTRAIN_EXTRA_OPTS = 

# the name of the output model that we use for the human analysis
TRAINED_MODEL_FOFN=ecoli_er2925.pcr_MSssI.timp.021216.alphabet_cpg.fofn

# The log-likelihood threshold required to make a call
CALL_THRESHOLD=2.5

#
# Download data
#

# In this Makefile we default to wget for portability.
# It is highly recommended that you use aspera instead as the transfer will
# go much faster. If you want to use aspera change this variable to the path
# to a bash script that implements the aspera download. See download_from_aspera.sh
# in this directory for an example that you will have to modify to use your key pair.

#DOWNLOADER=/path/to/download_from_aspera.sh
DOWNLOADER=wget

# Download a tar file from the ENA
$(DATA_ROOT)/%.tar:
	cd $(DATA_ROOT)	&& $(DOWNLOADER) ftp://ftp.sra.ebi.ac.uk/vol1/ERA540/ERA540530/oxfordnanopore_native/$(@F)

# untar it
$(DATA_ROOT)/%: $(DATA_ROOT)/%.tar
	cd $(DATA_ROOT) && tar -xf $*.tar

all-ecoli-k12-data: $(DATA_ROOT)/ecoli_k12.native.loman.240915.fast5 \
                    $(DATA_ROOT)/ecoli_k12.native.loman.280915.fast5 \
                    $(DATA_ROOT)/ecoli_k12.pcr.loman.250915.fast5 \
                    $(DATA_ROOT)/ecoli_k12.pcr.loman.280915.fast5

# symlink to structured names
$(DATA_ROOT)/ecoli_k12.native.loman.240915.fast5: $(DATA_ROOT)/MAP006-1
	cd $(DATA_ROOT) && ln -s MAP006-1/MAP006-1_downloads $(@F)

$(DATA_ROOT)/ecoli_k12.native.loman.280915.fast5: $(DATA_ROOT)/MAP006-2
	cd $(DATA_ROOT) && ln -s MAP006-2/MAP006-2_downloads $(@F)

$(DATA_ROOT)/ecoli_k12.pcr.loman.250915.fast5: $(DATA_ROOT)/MAP006-PCR-1
	cd $(DATA_ROOT) && ln -s MAP006-PCR-1/MAP006-PCR_downloads $(@F)

$(DATA_ROOT)/ecoli_k12.pcr.loman.280915.fast5: $(DATA_ROOT)/MAP006-PCR-2
	cd $(DATA_ROOT) && ln -s MAP006-PCR-2/MAP006-PCR-2_downloads $(@F)

# Reference genomes
$(ECOLI_REFERENCE):
	wget ftp://ftp.ensemblgenomes.org/pub/release-29/bacteria//fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.2.29.dna.genome.fa.gz
	gunzip Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.2.29.dna.genome.fa.gz
	mv Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.2.29.dna.genome.fa $@

$(HUMAN_REFERENCE):
	wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/GRCh38.primary_assembly.genome.fa.gz
	-gunzip GRCh38.primary_assembly.genome.fa.gz

#
# Rules for generating fasta files from FAST5
#

# Pass reads
%.pass.fasta: $(DATA_ROOT)/%.fast5/pass
	$(PORETOOLS) fasta --type 2D $< > $@

# Fail reads
%.fail.fasta: $(DATA_ROOT)/%.fast5/fail
	$(PORETOOLS) fasta --type 2D $< > $@

# all reads
%.fasta: %.pass.fasta %.fail.fasta
	cat $^ > $@

# all human runs 
$(HUMAN_NA12878_NATIVE_MERGED_DATA): NA12878.native.timp.093015.fasta \
                                     NA12878.native.simpson.101515.fasta \
                                     NA12878.native.simpson.103015.fasta
	cat $^ > $@

##################################################
#
# Step 2. Align each data set to a reference
#
##################################################

# Index a reference genome with BWA
%.fasta.bwt: %.fasta bwa.version
	bwa/bwa index $<

%.fa.bwt: %.fa bwa.version
	bwa/bwa index $<

# We use secondary expansion to construct the name of the variable
# containing the reference genome for this sample
.SECONDEXPANSION:
%.sorted.bam: %.fasta $$(%.fasta_REFERENCE) $$(%.fasta_REFERENCE).bwt bwa.version samtools.version
	bwa/bwa mem -t $(THREADS) -x ont2d $($*.fasta_REFERENCE) $< |\
        samtools/samtools view -Sb - |\
        samtools/samtools sort -f - $@

%.bam.bai: %.bam
	samtools/samtools index $<

%.eventalign.summary %.eventalign: %.sorted.bam %.sorted.bam.bai
	nanopolish/nanopolish eventalign --summary $*.eventalign.summary \
                                     --scale-events \
                                     --print-read-names \
                                     -b $*.sorted.bam \
                                     -r $*.fasta \
                                     -g $($*.fasta_REFERENCE) \
                                     -t $(THREADS) > $*.eventalign

##################################################
#
# Step 3. Train models
#
##################################################

# symlink the ont models from the analysis directory
t.006.ont.model: $(SCRIPT_DIR)/models/r7.3_e6_70bps_6mer_template_median68pA.model
	ln -s $^ $@

c.p1.006.ont.model: $(SCRIPT_DIR)/models/r7.3_e6_70bps_6mer_complement_median68pA_pop1.model
	ln -s $^ $@

c.p2.006.ont.model: $(SCRIPT_DIR)/models/r7.3_e6_70bps_6mer_complement_median68pA_pop2.model
	ln -s $^ $@

# Make a file-of-filenames from a model set
%.fofn: t.006.%.model c.p1.006.%.model c.p2.006.%.model
	echo $^ | tr " " "\n" > $@

# this code generates a training rule for a data set/alphabet pair 
define generate-training-rules

$(eval DATASET=$(basename $1))
$(eval ALPHABET=$2)
$(eval PREFIX=$(DATASET).alphabet_$(ALPHABET))

# Training rule
$(PREFIX).fofn: $(DATASET).fasta $(DATASET).sorted.bam $(DATASET).sorted.bam.bai $(TRAINING_REFERENCE).alphabet_$(ALPHABET) ont.alphabet_$(ALPHABET).fofn nanopolish.version
	nanopolish/nanopolish methyltrain -t $(THREADS) $(METHYLTRAIN_EXTRA_OPTS) \
                                      --train-kmers all \
                                      --out-fofn $$@ \
                                      --out-suffix .$(PREFIX).model \
                                      -m ont.alphabet_$(ALPHABET).fofn \
                                      -b $(DATASET).sorted.bam \
                                      -r $(DATASET).fasta \
                                      -g $(TRAINING_REFERENCE).alphabet_$(ALPHABET) \
                                      $(TRAINING_REGION)

# Rules for additional files made by methyltrain
t.006.$(PREFIX).model: $(PREFIX).fofn
c.p1.006.$(PREFIX).model: $(PREFIX).fofn
c.p2.006.$(PREFIX).model: $(PREFIX).fofn
$(DATASET).sorted.bam.methyltrain.tsv: $(PREFIX).fofn
methyltrain.$(PREFIX).model.summary: $(PREFIX).fofn

endef

#
# 3a. Train over a normal nucleotide (ACGT) alphabet
#

# Special case, the expanded reference for the alphabet is just the reference
$(TRAINING_REFERENCE).alphabet_nucleotide: $(TRAINING_REFERENCE)
	ln -s $< $@

# Special case, the expanded model for the nucleotide alphabet is just the ONT model
%.alphabet_nucleotide.model: %.model
	ln -s $< $@

# make the rules using the generation function
$(eval $(call generate-training-rules,$(ECOLI_K12_NATIVE_RUN1_DATA),nucleotide))
$(eval $(call generate-training-rules,$(ECOLI_K12_PCR_RUN1_DATA),nucleotide))
$(eval $(call generate-training-rules,$(ECOLI_ER2925_NATIVE_RUN1_DATA),nucleotide))
$(eval $(call generate-training-rules,$(ECOLI_ER2925_NATIVE_RUN2_DATA),nucleotide))
$(eval $(call generate-training-rules,$(ECOLI_ER2925_PCR_RUN1_DATA),nucleotide))
$(eval $(call generate-training-rules,$(ECOLI_ER2925_PCR_RUN2_DATA),nucleotide))
$(eval $(call generate-training-rules,$(ECOLI_ER2925_PCR_MSSSI_RUN1_DATA),nucleotide))
$(eval $(call generate-training-rules,$(ECOLI_ER2925_PCR_MSSSI_RUN2_DATA),nucleotide))

#
# 3b. Train over a CpG alphabet
#

# Expand the alphabet to include 5-mC k-mers
%.alphabet_cpg.model: %.model
	python $(SCRIPT_DIR)/expand_model_alphabet.py --alphabet cpg $< > $@

# Convert all CG dinucleotides of the reference genome to MG
$(TRAINING_REFERENCE).alphabet_cpg: $(TRAINING_REFERENCE) pythonlibs.version
	python $(SCRIPT_DIR)/methylate_reference.py --recognition cpg $< > $@

# make the rules using the generation function
$(eval $(call generate-training-rules,$(ECOLI_ER2925_MSSSI_DATA),cpg))
$(eval $(call generate-training-rules,$(ECOLI_ER2925_NATIVE_RUN1_DATA),cpg))
$(eval $(call generate-training-rules,$(ECOLI_ER2925_NATIVE_RUN2_DATA),cpg))
$(eval $(call generate-training-rules,$(ECOLI_K12_PCR_RUN1_DATA),cpg))
$(eval $(call generate-training-rules,$(ECOLI_ER2925_PCR_RUN1_DATA),cpg))
$(eval $(call generate-training-rules,$(ECOLI_ER2925_PCR_RUN2_DATA),cpg))
$(eval $(call generate-training-rules,$(ECOLI_ER2925_PCR_MSSSI_RUN1_DATA),cpg))
$(eval $(call generate-training-rules,$(ECOLI_ER2925_PCR_MSSSI_RUN2_DATA),cpg))

#
# 3c. Make training plots
#
training_plots_abcMG_event_mean.pdf: $(MSSSI_TRAINING_BAM).methyltrain.tsv $(PCR_TRAINING_BAM).methyltrain.tsv
	Rscript $(SCRIPT_DIR)/methylation_plots.R training_plots $^

#
# 3d. Make training tables from the summary files
#

# PCR data

# nucleotide
results/pcr.nucleotide.training.table.tex: methyltrain.$(ECOLI_ER2925_PCR_RUN1_DATA:.fasta=.alphabet_nucleotide.model.summary) \
                                           methyltrain.$(ECOLI_ER2925_PCR_RUN2_DATA:.fasta=.alphabet_nucleotide.model.summary) \
                                           methyltrain.$(ECOLI_K12_PCR_RUN1_DATA:.fasta=.alphabet_nucleotide.model.summary)
	mkdir -p results
	python $(SCRIPT_DIR)/generate_training_table.py --treatment pcr --alphabet nucleotide $^ > $@

# cpg
results/pcr.cpg.training.table.tex: methyltrain.$(ECOLI_ER2925_PCR_RUN1_DATA:.fasta=.alphabet_cpg.model.summary) \
                                           methyltrain.$(ECOLI_ER2925_PCR_RUN2_DATA:.fasta=.alphabet_cpg.model.summary) \
                                           methyltrain.$(ECOLI_K12_PCR_RUN1_DATA:.fasta=.alphabet_cpg.model.summary)
	mkdir -p results
	python $(SCRIPT_DIR)/generate_training_table.py --treatment pcr --alphabet cpg $^ > $@

# PCR + M.SssI

# nucleotide
results/pcr_MSssI.nucleotide.training.table.tex: methyltrain.$(ECOLI_ER2925_PCR_MSSSI_RUN1_DATA:.fasta=.alphabet_nucleotide.model.summary) \
                                                 methyltrain.$(ECOLI_ER2925_PCR_MSSSI_RUN2_DATA:.fasta=.alphabet_nucleotide.model.summary)
	mkdir -p results
	python $(SCRIPT_DIR)/generate_training_table.py --treatment pcr_MSssI --alphabet nucleotide $^ > $@

# CpG
results/pcr_MSssI.cpg.training.table.tex: methyltrain.$(ECOLI_ER2925_PCR_MSSSI_RUN1_DATA:.fasta=.alphabet_cpg.model.summary) \
                                          methyltrain.$(ECOLI_ER2925_PCR_MSSSI_RUN2_DATA:.fasta=.alphabet_cpg.model.summary)
	mkdir -p results
	python $(SCRIPT_DIR)/generate_training_table.py --treatment pcr_MSssI --alphabet cpg $^ > $@

# Merge all tables
results/all.training.tables.tex: results/pcr.nucleotide.training.table.tex \
                                 results/pcr_MSssI.nucleotide.training.table.tex \
                                 results/pcr.cpg.training.table.tex \
                                 results/pcr_MSssI.cpg.training.table.tex
	cat $^ > $@
##################################################
#
# Step 4. Human genome analysis
#
##################################################

# Calculate methylation likelihood ratios at every CpG covered by nanopore reads
%.methyltest.sites.bed %.methyltest.reads.tsv %.methyltest.strand.tsv: % %.bai $(TRAINED_MODEL_FOFN)
	$(eval TMP_BAM = $<)
	$(eval TMP_FASTA = $(TMP_BAM:.sorted.bam=.fasta))
	$(eval TMP_REF = $($(TMP_FASTA)_REFERENCE))
	nanopolish/nanopolish methyltest  -t $(THREADS) \
        -m $(TRAINED_MODEL_FOFN) \
        -b $(TMP_BAM) \
        -r $(TMP_FASTA) \
        -g $(TMP_REF)

# Convert a site BED file into a tsv file for R
%.methyltest.sites.tsv: %.methyltest.sites.bed
	cat $< | python $(SCRIPT_DIR)/annotated_bed_to_tsv.py > $@

# Download a bed file summarizing an NA12878 bisulfite experiment from ENCODE
BISULFITE_BED=ENCFF279HCL.bed
$(BISULFITE_BED):
	wget https://www.encodeproject.org/files/ENCFF279HCL/@@download/$(BISULFITE_BED).gz
	gunzip $(BISULFITE_BED).gz

#
# CpG Island Methylation Results
#

CGI_FILE=$(SCRIPT_DIR)/annotations/cgisl_hg38.bed.gz
PROMOTER_FILE=$(SCRIPT_DIR)/annotations/gencode.v24.promoter.bed.gz
CGI_PROMOTER_BED=cpg_islands.promoter.bed

# Annotate the CpG islands with whether they are <= 2kb upstream of a gene
$(CGI_PROMOTER_BED): $(CGI_FILE) $(PROMOTOR_FILE) bedtools.version
	bedtools/bin/bedtools map -o first -c 4 -a <(zcat $(CGI_FILE) | bedtools/bin/bedtools sort) \
                                            -b <(zcat $(PROMOTER_FILE) | bedtools/bin/bedtools sort) | \
                                            sed s/:_/=/ | awk '{ print $$1 "\t" $$2 "\t" $$3 "\t" $$4 ";Feature=" $$5 }' > $@


# Calculate a summary data for each CpG island from the bisulfite data
NA12878.bisulfite_score.cpg_islands: $(CGI_PROMOTER_BED) $(BISULFITE_BED) bedtools.version
	bedtools/bin/bedtools intersect -wb -b $(CGI_PROMOTER_BED) -a $(BISULFITE_BED) | \
        python $(SCRIPT_DIR)/calculate_methylation_at_cpg_islands.py -t bisulfite > $@

# Calculate a summary score for each CpG island from the ONT reads
%.ont_score.cpg_islands: $(CGI_PROMOTER_BED) %.sorted.bam.methyltest.sites.bed bedtools.version
	bedtools/bin/bedtools intersect -wb -b $(CGI_PROMOTER_BED) \
                                        -a $*.sorted.bam.methyltest.sites.bed | \
        python $(SCRIPT_DIR)/calculate_methylation_at_cpg_islands.py -t ont -c $(CALL_THRESHOLD) > $@

%.cpg_island_plot.pdf: NA12878.bisulfite_score.cpg_islands %.ont_score.cpg_islands
	Rscript $(SCRIPT_DIR)/methylation_plots.R human_cpg_island_plot $^ $@

#
# TSS results
#

TSS_FILE=$(SCRIPT_DIR)/annotations/gencode.v24.tss.bed.gz

# Calculate methylation as a function of distance from a TSS for the ONT data
%.methylated_sites.distance_to_TSS.bed: %.sorted.bam.methyltest.sites.bed bedtools.version $(TSS_FILE)
	bedtools/bin/bedtools closest -D b -b <(zcat $(TSS_FILE) | bedtools/bin/bedtools sort)\
                                       -a <(cat $*.sorted.bam.methyltest.sites.bed | bedtools/bin/bedtools sort) > $@


%.methylated_sites.distance_to_TSS.table: %.methylated_sites.distance_to_TSS.bed
	python $(SCRIPT_DIR)/calculate_methylation_by_distance.py --type ont -c $(CALL_THRESHOLD) -i $^ > $@

# Calculate methylation as a function of distance for the bisulfite data
NA12878.bisulfite.distance_to_TSS.bed: $(BISULFITE_BED) bedtools.version $(TSS_FILE)
	bedtools/bin/bedtools closest -D b -b <(zcat $(TSS_FILE) | bedtools/bin/bedtools sort)\
                                       -a <(cat $(BISULFITE_BED) | bedtools/bin/bedtools sort) > $@

%.bisulfite.distance_to_TSS.table: %.bisulfite.distance_to_TSS.bed
	python $(SCRIPT_DIR)/calculate_methylation_by_distance.py --type bisulfite -i $^ > $@


methylation_by_TSS_distance.pdf: NA12878.bisulfite.distance_to_TSS.table \
                                 NA12878.native.merged.methylated_sites.distance_to_TSS.table \
                                 NA12878.pcr.simpson.021616.methylated_sites.distance_to_TSS.table \
                                 NA12878.pcr_MSssI.simpson.021016.methylated_sites.distance_to_TSS.table
	Rscript $(SCRIPT_DIR)/methylation_plots.R TSS_distance_plot $^ $@


methylation_by_TSS_distance_by_chromosome.pdf: NA12878.bisulfite.distance_to_TSS.table \
                                               NA12878.native.merged.methylated_sites.distance_to_TSS.table
	Rscript $(SCRIPT_DIR)/methylation_plots.R TSS_distance_plot_by_chromosome $^ $@


#
# Site likelihood plot
#
site.likelihood.distribution.pdf: NA12878.pcr.simpson.021616.sorted.bam.methyltest.sites.tsv \
                     NA12878.pcr_MSssI.simpson.021016.sorted.bam.methyltest.sites.tsv \
                     NA12878.native.merged.sorted.bam.methyltest.sites.tsv
	Rscript $(SCRIPT_DIR)/methylation_plots.R site_likelihood_distribution $^ $@

#
# Accuracy results
#
accuracy.by_threshold.tsv \
accuracy.precision_recall.tsv \
accuracy.by_kmer.tsv: NA12878.pcr.simpson.021616.sorted.bam.methyltest.sites.bed \
                      NA12878.pcr_MSssI.simpson.021016.sorted.bam.methyltest.sites.bed
	python $(SCRIPT_DIR)/calculate_call_accuracy.py --unmethylated NA12878.pcr.simpson.021616.sorted.bam.methyltest.sites.bed \
                                                    --methylated NA12878.pcr_MSssI.simpson.021016.sorted.bam.methyltest.sites.bed

accuracy.roc.pdf: accuracy.precision_recall.tsv
	Rscript $(SCRIPT_DIR)/methylation_plots.R call_accuracy_roc $^ $@

accuracy.by_threshold.pdf: accuracy.by_threshold.tsv
	Rscript $(SCRIPT_DIR)/methylation_plots.R call_accuracy_by_threshold $^ $@

accuracy.by_kmer.pdf: accuracy.by_kmer.tsv
	Rscript $(SCRIPT_DIR)/methylation_plots.R call_accuracy_by_kmer $^ $@

