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
all: training_plots_abcMG_event_mean.pdf ProHum20kb_cpg_island_plot.pdf
all-trained: M.SssI.ecoli_e2925.sqk006.sorted.bam.methyltrain.tsv loman.ecoli_k12.pcr.sqk006.sorted.bam.methyltrain.tsv

##################################################
#
# Step 1. Prepare input data 
#
##################################################

# Convert a directory of FAST5 files to fasta using readtofasta.py
#%.fasta: %.fast5
#	python $(ROOT_DIR)/readtofasta.py $</ > $@

#
# Define variables for each data set
#
ECOLI_K12_DATA=loman.ecoli_k12.sqk006.fasta
ECOLI_K12_PCR_DATA=loman.ecoli_k12.pcr.sqk006.fasta
ECOLI_E2925_MSSSI_DATA=M.SssI.ecoli_e2925.sqk006.fasta
ECOLI_E2925_NATIVE_DATA=native.ecoli_e2925.sqk006.fasta

LAMBDA_MSSSI_DATA=M.SssI.lambda.fasta
LAMBDA_CONTROL_DATA=control.lambda.fasta

#HUMAN_NA12878_DATA=093015.NA12878.fasta
#HUMAN_NA12878_DATA=oicr.run1.wash.NA12878.fasta
HUMAN_NA12878_DATA=NA12878.merged.fasta

# For each data set that we use, define a variable containing its reference
$(ECOLI_K12_DATA)_REFERENCE=ecoli_k12.fasta
$(ECOLI_K12_PCR_DATA)_REFERENCE=ecoli_k12.fasta
$(ECOLI_E2925_MSSSI_DATA)_REFERENCE=ecoli_k12.fasta
$(ECOLI_E2925_NATIVE_DATA)_REFERENCE=ecoli_k12.fasta

#$(LAMBDA_MSSI_DATA)_REFERENCE=lambda.reference.fasta
#$(LAMBDA_CONTROL_DATA)_REFERENCE=lambda.reference.fasta

$(HUMAN_NA12878_DATA)_REFERENCE=human_g1k_v37.fasta

#
# These variables control which datasets are used to train the model, test, etc
#
PCR_TRAINING_FASTA=$(ECOLI_K12_PCR_DATA)
DAM_TRAINING_FASTA=$(ECOLI_K12_DATA)
DCM_TRAINING_FASTA=$(ECOLI_K12_DATA)
MSSSI_TRAINING_FASTA=$(ECOLI_E2925_MSSSI_DATA)
MSSSI_CONTROL_FASTA=$(ECOLI_E2925_NATIVE_DATA)
TRAINING_REGION="gi|556503834|ref|NC_000913.3|:50000-3250000"

#TEST_FASTA=$(LAMBDA_MSSSI_DATA)
#TEST_CONTROL_FASTA=$(LAMBDA_CONTROL_DATA)

#
# These are derived convenience variables and should not be manually set
# 
# Derive dependent file names

# Reference file for training, must be the same for all training sets
TRAINING_REFERENCE=$($(PCR_TRAINING_FASTA)_REFERENCE)

# BAM file names
PCR_TRAINING_BAM=$(PCR_TRAINING_FASTA:.fasta=.sorted.bam)
DAM_TRAINING_BAM=$(DAM_TRAINING_FASTA:.fasta=.sorted.bam)
DCM_TRAINING_BAM=$(DCM_TRAINING_FASTA:.fasta=.sorted.bam)
MSSSI_TRAINING_BAM=$(MSSSI_TRAINING_FASTA:.fasta=.sorted.bam)
MSSSI_CONTROL_BAM=$(MSSSI_CONTROL_FASTA:.fasta=.sorted.bam)

TEST_REFERENCE=$($(TEST_FASTA)_REFERENCE)
TEST_BAM=$(TEST_FASTA:.fasta=.sorted.bam)
TEST_CONTROL_BAM=$(TEST_CONTROL_FASTA:.fasta=.sorted.bam)

# Convert a FAST5 file to FASTA using poretools
PORETOOLS=poretools

# Pass reads
%.pass.fasta: %.fast5/pass
	$(PORETOOLS) fasta --type 2D $< > $@

# Fail reads
%.fail.fasta: %.fast5/fail
	$(PORETOOLS) fasta --type 2D $< > $@

%.fasta: %.pass.fasta %.fail.fasta
	cat $^ > $@

# Special case for data sets with multiple runs that must be joined together
$(ECOLI_E2925_MSSSI_DATA): M.SssI.ecoli_e2925.sqk006.run1.fasta M.SssI.ecoli_e2925.sqk006.run2.fasta
	cat $^ > $@

$(HUMAN_NA12878_DATA): 093015.NA12878.fasta \
                       oicr.run1.wash.NA12878.fasta \
                       oicr.run2.NA12878.fasta
	cat $^ > $@

##################################################
#
# Step 2. Align each data set to a reference
#
##################################################

# Index a reference genome with BWA
%.fasta.bwt: %.fasta
	bwa/bwa index $<

# We use secondary expansion to construct the name of the variable
# containing the reference genome for this sample
.SECONDEXPANSION:
%.sorted.bam: %.fasta $$(%.fasta_REFERENCE) $$(%.fasta_REFERENCE).bwt bwa.version samtools.version
	bwa/bwa mem -t $(THREADS) -x ont2d $($*.fasta_REFERENCE) $< | \
        samtools/samtools view -Sb - | \
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

# Make a file-of-filenames from a model set
%.fofn: t.006.%.model c.p1.006.%.model c.p2.006.%.model
	echo $^ | tr " " "\n" > $@

#
# 3a. First, train using PCR-treated data. The
#     trained model should be similar to ONT's
#     model.
#

#
PCR_ECOLI_TRAINED=pcr.ecoli_k12.trained
$(PCR_ECOLI_TRAINED).fofn: $(PCR_TRAINING_BAM) $(PCR_TRAINING_BAM:.bam=.bam.bai) $(PCR_TRAINING_FASTA) $(TRAINING_REFERENCE) ont.fofn
	nanopolish/nanopolish methyltrain -t $(THREADS) \
                                      --train-kmers all \
                                      --out-fofn $@ \
                                      --out-suffix .$(PCR_ECOLI_TRAINED).model \
                                      -m ont.fofn \
                                      -b $(PCR_TRAINING_BAM) \
                                      -r $(PCR_TRAINING_FASTA) \
                                      -g $(TRAINING_REFERENCE) \
                                      $(TRAINING_REGION) 

# These files are built at the same time as output fofn. These rules make sure they get updated.
t.006.$(PCR_ECOLI_TRAINED).model: $(PCR_ECOLI_TRAINED).fofn
c.p1.006.$(PCR_ECOLI_TRAINED).model: $(PCR_ECOLI_TRAINED).fofn
c.p2.006.$(PCR_ECOLI_TRAINED).model: $(PCR_ECOLI_TRAINED).fofn
$(PCR_TRAINING_BAM).methyltrain.tsv: $(PCR_ECOLI_TRAINED).fofn

#
# 3b. Train the mCpG using the M.SssI treated ecoli data
#

# Expand the alphabet to include 5-mC k-mers
%.cpg_expanded.model: %.model
	python $(ROOT_DIR)/expand_model_alphabet.py --alphabet CpG $< > $@

# Convert all CG dinucleotides of the reference genome to MG
$(TRAINING_REFERENCE).cpg_methylated: $(TRAINING_REFERENCE) pythonlibs.version
	python $(ROOT_DIR)/methylate_reference.py --recognition CpG $< > $@

# Train the model
MSSSI_ECOLI_TRAINED=M.SssI.ecoli_e2925.trained
$(MSSSI_ECOLI_TRAINED).fofn: $(MSSSI_TRAINING_BAM) $(MSSSI_TRAINING_BAM:.bam=.bam.bai) $(MSSSI_TRAINING_FASTA) $(TRAINING_REFERENCE).cpg_methylated ont.cpg_expanded.fofn
	nanopolish/nanopolish methyltrain -t $(THREADS) \
                                      --train-kmers all \
                                      --out-fofn $@ \
                                      --out-suffix .$(MSSSI_ECOLI_TRAINED).model \
                                      -m ont.cpg_expanded.fofn \
                                      -b $(MSSSI_TRAINING_BAM) \
                                      -r $(MSSSI_TRAINING_FASTA) \
                                      -g $(TRAINING_REFERENCE).cpg_methylated \
                                      $(TRAINING_REGION)

t.006.$(MSSSI_ECOLI_TRAINED).model: $(MSSSI_ECOLI_TRAINED).fofn
c.p1.006.$(MSSSI_ECOLI_TRAINED).model: $(MSSSI_ECOLI_TRAINED).fofn
c.p2.006.$(MSSSI_ECOLI_TRAINED).model: $(MSSSI_ECOLI_TRAINED).fofn
$(MSSSI_TRAINING_BAM).methyltrain.tsv: $(MSSSI_ECOLI_TRAINED).fofn

# Perform the training with the unmethylated E2925 data as a negative control
MSSSI_CONTROL_TRAINED=native.ecoli_e2925.trained
$(MSSSI_CONTROL_TRAINED).fofn: $(MSSSI_CONTROL_BAM) $(MSSSI_CONTROL_BAM:.bam=.bam.bai) $(MSSSI_CONTROL_FASTA) $(TRAINING_REFERENCE).cpg_methylated ont.cpg_expanded.fofn
	nanopolish/nanopolish methyltrain -t $(THREADS) \
                                      --train-kmers all \
                                      --out-fofn $@ \
                                      --out-suffix .$(MSSSI_CONTROL_TRAINED).model \
                                      -m ont.cpg_expanded.fofn \
                                      -b $(MSSSI_CONTROL_BAM) \
                                      -r $(MSSSI_CONTROL_FASTA) \
                                      -g $(TRAINING_REFERENCE).cpg_methylated \
                                      $(TRAINING_REGION)

t.006.$(MSSSI_CONTROL_TRAINED).model: $(MSSSI_CONTROL_TRAINED).fofn
c.p1.006.$(MSSSI_CONTROL_TRAINED).model: $(MSSSI_CONTROL_TRAINED).fofn
c.p2.006.$(MSSSI_CONTROL_TRAINED).model: $(MSSSI_CONTROL_TRAINED).fofn
$(MSSSI_CONTROL_BAM).methyltrain.tsv: $(MSSSI_CONTROL_TRAINED).fofn


PCR_ECOLI_CPG_TRAINED=pcr.ecoli_k12.cpg.trained
$(PCR_ECOLI_CPG_TRAINED).fofn: $(PCR_TRAINING_BAM) $(PCR_TRAINING_BAM:.bam=.bam.bai) $(PCR_TRAINING_FASTA) $(TRAINING_REFERENCE).cpg_methylated ont.fofn
	nanopolish/nanopolish methyltrain -t $(THREADS) \
                                      --train-kmers all \
                                      --out-fofn $@ \
                                      --out-suffix .$(PCR_ECOLI_CPG_TRAINED).model \
                                      -m ont.cpg_expanded.fofn \
                                      -b $(PCR_TRAINING_BAM) \
                                      -r $(PCR_TRAINING_FASTA) \
                                      -g $(TRAINING_REFERENCE).cpg_methylated \
                                      $(TRAINING_REGION) 

# Perform the training with the unmethylated E2925 data as a negative control
MSSSI_CONTROL_K12_TRAINED=native.ecoli_e2925.k12.trained
$(MSSSI_CONTROL_K12_TRAINED).fofn: $(MSSSI_CONTROL_BAM) $(MSSSI_CONTROL_BAM:.bam=.bam.bai) $(MSSSI_CONTROL_FASTA) $(TRAINING_REFERENCE) ont.fofn
	nanopolish/nanopolish methyltrain -t $(THREADS) \
                                      --train-kmers all \
                                      --out-fofn $@ \
                                      --out-suffix .$(MSSSI_CONTROL_K12_TRAINED).model \
                                      -m ont.fofn \
                                      -b $(MSSSI_CONTROL_BAM) \
                                      -r $(MSSSI_CONTROL_FASTA) \
                                      -g $(TRAINING_REFERENCE) \
                                      $(TRAINING_REGION)

t.006.$(MSSSI_CONTROL_TRAINED).model: $(MSSSI_CONTROL_TRAINED).fofn
c.p1.006.$(MSSSI_CONTROL_TRAINED).model: $(MSSSI_CONTROL_TRAINED).fofn
c.p2.006.$(MSSSI_CONTROL_TRAINED).model: $(MSSSI_CONTROL_TRAINED).fofn
$(MSSSI_CONTROL_BAM).methyltrain.tsv: $(MSSSI_CONTROL_TRAINED).fofn


#
# 3c. Train the dam model
#

# Expand the alphabet of the pcr-trained model to include 5-mC k-mers
%.dam_expanded.model: %.model
	python $(ROOT_DIR)/expand_model_alphabet.py --alphabet dam $< > $@

# Convert all GATC sequences of the reference genome to GMTC
$(TRAINING_REFERENCE).dam.methylated: $(TRAINING_REFERENCE) pythonlibs.version
	python $(ROOT_DIR)/methylate_reference.py --recognition dam $< > $@

# Train the model with methylated 5-mers
DAM_ECOLI_TRAINED=dam.ecoli_k12.trained
$(DAM_ECOLI_TRAINED).fofn: $(DAM_TRAINING_BAM) $(DAM_TRAINING_BAM:.bam=.bam.bai) $(DAM_TRAINING_FASTA) $(TRAINING_REFERENCE).dam.methylated ont.dam_expanded.fofn
	nanopolish/nanopolish methyltrain -t $(THREADS) \
                                      --train-kmers all \
                                      --out-fofn $@ \
                                      --out-suffix .$(DAM_ECOLI_TRAINED).model \
                                      -m ont.dam_expanded.fofn \
                                      -b $(DAM_TRAINING_BAM) \
                                      -r $(DAM_TRAINING_FASTA) \
                                      -g $(TRAINING_REFERENCE).dam.methylated \
                                      $(TRAINING_REGION)

t.006.$(DAM_ECOLI_TRAINED).model: $(DAM_ECOLI_TRAINED).fofn
c.p1.006.$(DAM_ECOLI_TRAINED).model: $(DAM_ECOLI_TRAINED).fofn
c.p2.006.$(DAM_ECOLI_TRAINED).model: $(DAM_ECOLI_TRAINED).fofn
$(DAM_TRAINING_BAM).methyltrain.tsv: $(DAM_ECOLI_TRAINED).fofn

#
# 3d. Train the dcm model
#

#
%.dcm_expanded.model: %.model
	python $(ROOT_DIR)/expand_model_alphabet.py --alphabet dcm $< > $@

# Convert all dcm sequences of the reference genome to GMTC
$(TRAINING_REFERENCE).dcm.methylated: $(TRAINING_REFERENCE) pythonlibs.version
	python $(ROOT_DIR)/methylate_reference.py --recognition dcm $< > $@

# Train the model with methylated 5-mers
DCM_ECOLI_TRAINED=dcm.ecoli_k12.trained
$(DCM_ECOLI_TRAINED).fofn: $(DCM_TRAINING_BAM) $(DCM_TRAINING_BAM:.bam=.bam.bai) $(DCM_TRAINING_FASTA) $(TRAINING_REFERENCE).dcm.methylated ont.dcm_expanded.fofn
	nanopolish/nanopolish methyltrain -t $(THREADS) \
                                      --train-kmers all \
                                      --out-fofn $@ \
                                      --out-suffix .$(DCM_ECOLI_TRAINED).model \
                                      -m ont.dcm_expanded.fofn \
                                      -b $(DCM_TRAINING_BAM) \
                                      -r $(DCM_TRAINING_FASTA) \
                                      -g $(TRAINING_REFERENCE).dcm.methylated \
                                      $(TRAINING_REGION)

t.006.$(DCM_ECOLI_TRAINED).model: $(DCM_ECOLI_TRAINED).fofn
c.p1.006.$(DCM_ECOLI_TRAINED).model: $(DCM_ECOLI_TRAINED).fofn
c.p2.006.$(DCM_ECOLI_TRAINED).model: $(DCM_ECOLI_TRAINED).fofn
$(DCM_TRAINING_BAM).methyltrain.tsv: $(DCM_ECOLI_TRAINED).fofn

# Make training plots
training_plots_abcMG_event_mean.pdf: $(MSSSI_TRAINING_BAM).methyltrain.tsv $(PCR_TRAINING_BAM).methyltrain.tsv
	Rscript $(ROOT_DIR)/methylation_plots.R training_plots $^
	cp $@ $@.$(NOW).pdf

##################################################
#
# Step 4. Test the methylation model
#
##################################################
%.methyltest.sites.bed %.methyltest.reads.tsv %.methyltest.strand.tsv: % %.bai M.SssI.ecoli_e2925.trained.fofn
	$(eval TMP_BAM = $<)
	$(eval TMP_FASTA = $(TMP_BAM:.sorted.bam=.fasta))
	$(eval TMP_REF = $($(TMP_FASTA)_REFERENCE))
	nanopolish/nanopolish methyltest  -t $(THREADS) \
                                      -m M.SssI.ecoli_e2925.trained.fofn \
                                      -b $(TMP_BAM) \
                                      -r $(TMP_FASTA) \
                                      -g $(TMP_REF)

# Convert a site BED file into a tsv file for R
%.methyltest.sites.tsv: %.methyltest.sites.bed
	cat $< | python $(ROOT_DIR)/annotated_bed_to_tsv.py > $@

site_likelihood_plots.pdf: $(TEST_BAM).methyltest.sites.tsv
	Rscript $(ROOT_DIR)/methylation_plots.R site_likelihood_plots $< $@
	cp $@ $@.$(NOW).pdf

read_classification_plot.pdf: $(TEST_BAM).methyltest.reads.tsv $(TEST_CONTROL_BAM).methyltest.reads.tsv
	Rscript $(ROOT_DIR)/methylation_plots.R read_classification_plot $^ $@
	cp $@ $@.$(NOW).pdf

strand_classification_plot.pdf: $(TEST_BAM).methyltest.strand.tsv $(TEST_CONTROL_BAM).methyltest.strand.tsv
	Rscript $(ROOT_DIR)/methylation_plots.R strand_classification_plot $^ $@
	cp $@ $@.$(NOW).pdf


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
	cp $@ $@.$(NOW).pdf
	cp histogram.pdf $*.cpg_histogram.$(NOW).pdf
	cp histogram_chromosomes.pdf $*.cpg_histogram_chromosomes.$(NOW).pdf

%.site_comparison.tsv: ENCFF257GGV.bed %.sorted.bam.methyltest.sites.bed
	bedtools/bin/bedtools intersect -wb -a ENCFF257GGV.bed \
                                        -b <(awk '{ print "chr" $$0 }' $*.sorted.bam.methyltest.sites.bed) | \
                                        python ~/simpsonlab/users/jsimpson/code/methylation-analysis/calculate_bisulfite_vs_ont_signal.py > $@

%.site_comparison.pdf: %.site_comparison.tsv
	Rscript $(ROOT_DIR)/methylation_plots.R site_comparison_plot $^ $@
	cp $@ $@.$(NOW).pdf

##################################################
#
# Step 6. Global methylation analysis
#
##################################################
global_methylation.pdf: ProHum20kb.sorted.bam.methyltest.sites.tsv giab.NA24385.sorted.bam.methyltest.sites.tsv control.lambda.sorted.bam.methyltest.sites.tsv M.SssI.lambda.sorted.bam.methyltest.sites.tsv
	Rscript $(ROOT_DIR)/methylation_plots.R global_methylation $^ $@
	cp $@ $@.$(NOW).pdf
