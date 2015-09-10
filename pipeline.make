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
ROOT_DIR:=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

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
	pip install --user biopython >> $@
	pip install --user h5py >> $@
	pip freeze >> $@

# Install bedtools
bedtools.version:
	git clone https://github.com/arq5x/bedtools2.git bedtools
	cd bedtools; make
	-cd bedtools; git log | head -1 > ../$@

# timestamp for saving different versions of plots
# http://stackoverflow.com/questions/12463770/insert-time-stamp-into-executable-name-using-makefile
NOW := $(shell date +'%y.%m.%d_%H:%M:%S')

.DEFAULT_GOAL := all
all: training_plots.pdf site_likelihood_plots.pdf read_classification_plot.pdf ProHum20kb_cpg_island_plot.pdf

##################################################
#
# Step 1. Prepare input data 
#
##################################################

# Convert a directory of FAST5 files to fasta using readtofasta.py
%.fasta: %.fast5
	python $(ROOT_DIR)/readtofasta.py $</ > $@

#
# Define variables for each data set
#
ECOLI_MSSI_DATA=M.SssI.e2925_ecoli.fasta
ECOLI_CONTROL_DATA=ERX708228.ecoli.fasta

LAMBDA_MSSI_DATA=M.SssI.lambda.fasta
LAMBDA_CONTROL_DATA=control.lambda.fasta

HUMAN_PROMEGA_DATA=ProHum20kb.fasta
HUMAN_GIAB_DATA=giab.NA24385.fasta

# For each data set that we use, define a variable containing its reference
$(ECOLI_MSSI_DATA)_REFERENCE=ecoli_k12.fasta
$(ECOLI_CONTROL_DATA)_REFERENCE=ecoli_k12.fasta

$(LAMBDA_MSSI_DATA)_REFERENCE=lambda.reference.fasta
$(LAMBDA_CONTROL_DATA)_REFERENCE=lambda.reference.fasta

$(HUMAN_PROMEGA_DATA)_REFERENCE=human_g1k_v37.fasta
$(HUMAN_GIAB_DATA)_REFERENCE=human_g1k_v37.fasta

#
# These variables control which datasets are used to train the model, test, etc
#
TRAINING_FASTA=$(ECOLI_MSSI_DATA)
TRAINING_CONTROL_FASTA=$(ECOLI_CONTROL_DATA)
TRAINING_REGION="gi|556503834|ref|NC_000913.3|:50000-700000"

TEST_FASTA=$(LAMBDA_MSSI_DATA)
TEST_CONTROL_FASTA=$(LAMBDA_CONTROL_DATA)

#
# These are derived convenience variables and should not be manually set
# 
# Derive dependent file names
TRAINING_REFERENCE=$($(TRAINING_FASTA)_REFERENCE)
TRAINING_BAM=$(TRAINING_FASTA:.fasta=.sorted.bam)
TRAINING_CONTROL_BAM=$(TRAINING_CONTROL_FASTA:.fasta=.sorted.bam)

TEST_REFERENCE=$($(TEST_FASTA)_REFERENCE)
TEST_BAM=$(TEST_FASTA:.fasta=.sorted.bam)
TEST_CONTROL_BAM=$(TEST_CONTROL_FASTA:.fasta=.sorted.bam)

##################################################
#
# Step 2. Align each data set to a reference
#
##################################################

# Index a reference genome with BWA
%.fasta.bwt: %.fasta bwa.version
	bwa/bwa index $<

# We use secondary expansion to construct the name of the variable
# containing the reference genome for this sample
.SECONDEXPANSION:
%.sorted.bam: %.fasta $$(%.fasta_REFERENCE) $$(%.fasta_REFERENCE).bwt bwa.version samtools.version
	make -f $(ROOT_DIR)/alignment.make OUTPUT=$@ \
                                       READS=$< \
                                       REFERENCE=$($*.fasta_REFERENCE) \
                                       THREADS=$(THREADS) \
                                       BWA=bwa/bwa \
                                       SAMTOOLS=samtools/samtools

%.bam.bai: %.bam
	samtools/samtools index $<

##################################################
#
# Step 3. Train the methylation model
#
##################################################

# Convert all CG dinucleotides of the reference genome to MG for training
$(TRAINING_REFERENCE).methylated: $(TRAINING_REFERENCE) pythonlibs.version
	python $(ROOT_DIR)/methylate_reference.py $< > $@

# Pretrain a model on unmethylated data to make the emissions better fit our HMM
r7.3_template_median68pA.model.pretrain.methyltrain \
r7.3_complement_median68pA_pop1.model.pretain.methyltrain \
r7.3_complement_median68pA_pop2.model.pretrain.methyltrain: $(TRAINING_BAM) $(TRAINING_BAM:.bam=.bam.bai) $(TRAINING_FASTA) $(TRAINING_REFERENCE).methylated initial_pretrain_models.fofn
	nanopolish/nanopolish methyltrain -t $(THREADS) \
                                      --train-unmethylated \
                                      --out-suffix ".pretrain.methyltrain" \
                                      -m initial_pretrain_models.fofn \
                                      -b $(TRAINING_CONTROL_BAM) \
                                      -r $(TRAINING_CONTROL_FASTA) \
                                      -g $(TRAINING_REFERENCE) \
                                      -w $(TRAINING_REGION) 

# Initialize methylation models from the base models
%.model.pretrain: %.model
	python $(ROOT_DIR)/methylate_model.py $< > $@

# Make a fofn of the initialized pretrain models
initial_pretrain_models.fofn: r7.3_template_median68pA.model.pretrain \
                              r7.3_complement_median68pA_pop1.model.pretrain \
                              r7.3_complement_median68pA_pop2.model.pretrain
	echo $^ | tr " " "\n" > $@

%.model.pretrain.initial_methyl: %.model.pretrain.methyltrain
	python $(ROOT_DIR)/methylate_model.py $< > $@

# Make a fofn of the initialized methylation from the pretrain models    
initial_methyl_models.fofn: r7.3_template_median68pA.model.pretrain.initial_methyl r7.3_complement_median68pA_pop1.model.pretrain.initial_methyl r7.3_complement_median68pA_pop2.model.pretrain.initial_methyl
	echo $^ | tr " " "\n" > $@

# Train the model with methylated 5-mers
r7.3_template_median68pA.model.methyltrain: $(TRAINING_BAM) $(TRAINING_BAM:.bam=.bam.bai) $(TRAINING_FASTA) $(TRAINING_REFERENCE).methylated initial_methyl_models.fofn
	nanopolish/nanopolish methyltrain -t $(THREADS) \
                                      -m initial_methyl_models.fofn \
                                      -b $(TRAINING_BAM) \
                                      -r $(TRAINING_FASTA) \
                                      -g $(TRAINING_REFERENCE).methylated \
                                      -w $(TRAINING_REGION)

# Run methyltrain on unmethylated data to get data to compare to
$(TRAINING_CONTROL_BAM).training.0.tsv $(TRAINING_CONTROL_BAM).methyltrain.0.tsv : $(TRAINING_CONTROL_BAM) $(TRAINING_CONTROL_BAM:.bam=.bam.bai) $(TRAINING_CONTROL_FASTA) $(TRAINING_REFERENCE) initial_methyl_models.fofn
	nanopolish/nanopolish methyltrain -t $(THREADS) \
                                      -m initial_methyl_models.fofn \
                                      --no-update-models \
                                      -b $(TRAINING_CONTROL_BAM) \
                                      -r $(TRAINING_CONTROL_FASTA) \
                                      -g $(TRAINING_REFERENCE) \
                                      -w $(TRAINING_REGION)

# Make a fofn of the trained methylation models 
trained_methyl_models.fofn: r7.3_template_median68pA.model.methyltrain r7.3_complement_median68pA_pop1.model.methyltrain r7.3_complement_median68pA_pop2.model.methyltrain
	echo $^ | tr " " "\n" > $@

# Make training plots
training_plots.pdf: r7.3_template_median68pA.model.methyltrain $(TRAINING_CONTROL_BAM).methyltrain.0.tsv
	Rscript $(ROOT_DIR)/methylation_plots.R training_plots
	cp $@ $@.$(NOW)

##################################################
#
# Step 4. Test the methylation model
#
##################################################
%.methyltest.sites.bed %.methyltest.reads.tsv: % %.bai trained_methyl_models.fofn
	$(eval TMP_BAM = $<)
	$(eval TMP_FASTA = $(TMP_BAM:.sorted.bam=.fasta))
	$(eval TMP_REF = $($(TMP_FASTA)_REFERENCE))
	nanopolish/nanopolish methyltest  -t 1 \
                                      -m trained_methyl_models.fofn \
                                      -b $(TMP_BAM) \
                                      -r $(TMP_FASTA) \
                                      -g $(TMP_REF)

# Convert a site BED file into a tsv file for R
%.methyltest.sites.tsv: %.methyltest.sites.bed
	cat $< | python $(ROOT_DIR)/annotated_bed_to_tsv.py > $@

site_likelihood_plots.pdf: $(TEST_BAM).methyltest.sites.tsv
	Rscript $(ROOT_DIR)/methylation_plots.R site_likelihood_plots $< $@
	cp $@ $@.$(NOW)

read_classification_plot.pdf: $(TEST_BAM).methyltest.reads.tsv $(TEST_CONTROL_BAM).methyltest.reads.tsv
	Rscript $(ROOT_DIR)/methylation_plots.R read_classification_plot $^ $@
	cp $@ $@.$(NOW)

##################################################
#
# Step 5. Human genome analysis
#
##################################################

# Download database of CpG islands from Irizarry's method
irizarry.cpg_islands.bed:
	wget http://rafalab.jhsph.edu/CGI/model-based-cpg-islands-hg19.txt
	cat model-based-cpg-islands-hg19.txt | python $(ROOT_DIR)/tsv_to_annotated_bed.py > $@

# Annotate the CpG islands with whether they are <= 2kb upstream of a gene
irizarry.cpg_islands.genes.bed: irizarry.cpg_islands.bed gencode_genes_2kb_upstream.bed bedtools.version
	bedtools/bin/bedtools map -o first -c 4 -a <(cat irizarry.cpg_islands.bed | bedtools/bin/bedtools sort) \
                                            -b <(cat gencode_genes_2kb_upstream.bed | bedtools/bin/bedtools sort) | \
                                            awk '{ print $$1 "\t" $$2 "\t" $$3 "\t" $$4 ";Gene=" $$5 }' > $@


# Download a bed file summarizing an NA12878 bisulfite experiment from ENCODE
ENCFF257GGV.bed:
	wget https://www.encodeproject.org/files/ENCFF257GGV/@@download/ENCFF257GGV.bed.gz
	gunzip ENCFF257GGV.bed.gz

# Calculate a summary data for each CpG island from the bisulfite data
NA12878.bisulfite_score.cpg_islands: irizarry.cpg_islands.genes.bed ENCFF257GGV.bed bedtools.version
	bedtools/bin/bedtools intersect -wb -b irizarry.cpg_islands.genes.bed -a ENCFF257GGV.bed | \
        python $(ROOT_DIR)/calculate_bisulfite_signal_for_cpg_islands.py > $@

# Calculate a summary score for each CpG island from the ONT reads
%.ont_score.cpg_islands: irizarry.cpg_islands.genes.bed %.sorted.bam.methyltest.sites.bed bedtools.version
	bedtools/bin/bedtools intersect -wb -b irizarry.cpg_islands.genes.bed \
                                        -a <(awk '{ print "chr" $$0 }' $*.sorted.bam.methyltest.sites.bed) | \
        python $(ROOT_DIR)/calculate_ont_signal_for_cpg_islands.py > $@

%_cpg_island_plot.pdf: NA12878.bisulfite_score.cpg_islands %.ont_score.cpg_islands
	Rscript $(ROOT_DIR)/methylation_plots.R human_cpg_island_plot $^ $@
	cp $@ $@.$(NOW)
	cp histogram.pdf $*.cpg_histogram.$(NOW).pdf

