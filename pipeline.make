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
all: training_plots_abcMG_event_mean.pdf

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

all-methylation-plots: NA12878.native.timp.093015.cpg_island_plot.pdf \
                       NA12878.native.simpson.101515.cpg_island_plot.pdf \
                       NA12878.native.simpson.103015.cpg_island_plot.pdf \
                       NA12878.native.merged.cpg_island_plot.pdf \
                       NA12878.pcr.simpson.021616.cpg_island_plot.pdf \
                       NA12878.pcr_MSssI.simpson.021016.cpg_island_plot.pdf

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
ECOLI_ER2925_PCR_RUN3_DATA=ecoli_er2925.pcr.timp.113015_recall.fasta
ECOLI_ER2925_PCR_MSSSI_RUN1_DATA=ecoli_er2925.pcr_MSssI.timp.113015.fasta
ECOLI_ER2925_PCR_MSSSI_RUN2_DATA=ecoli_er2925.pcr_MSssI.timp.021216.fasta

HUMAN_NA12878_NATIVE_RUN1_DATA=NA12878.native.timp.093015.fasta
HUMAN_NA12878_NATIVE_RUN2_DATA=NA12878.native.simpson.101515.fasta
HUMAN_NA12878_NATIVE_RUN3_DATA=NA12878.native.simpson.103015.fasta
HUMAN_NA12878_NATIVE_MERGED_DATA=NA12878.native.merged.fasta

HUMAN_NA12878_PCR_DATA=NA12878.pcr.simpson.021616.fasta
HUMAN_NA12878_PCR_MSSSI_DATA=NA12878.pcr_MSssI.simpson.021016.fasta

# For each data set that we use, define a variable containing its reference
$(ECOLI_K12_NATIVE_DATA)_REFERENCE=ecoli_k12.fasta
$(ECOLI_K12_PCR_RUN1_DATA)_REFERENCE=ecoli_k12.fasta

$(ECOLI_ER2925_NATIVE_RUN1_DATA)_REFERENCE=ecoli_k12.fasta
$(ECOLI_ER2925_NATIVE_RUN2_DATA)_REFERENCE=ecoli_k12.fasta

$(ECOLI_ER2925_MSSSI_DATA)_REFERENCE=ecoli_k12.fasta
$(ECOLI_ER2925_PCR_RUN1_DATA)_REFERENCE=ecoli_k12.fasta
$(ECOLI_ER2925_PCR_RUN2_DATA)_REFERENCE=ecoli_k12.fasta
$(ECOLI_ER2925_PCR_RUN3_DATA)_REFERENCE=ecoli_k12.fasta
$(ECOLI_ER2925_PCR_MSSSI_RUN1_DATA)_REFERENCE=ecoli_k12.fasta
$(ECOLI_ER2925_PCR_MSSSI_RUN2_DATA)_REFERENCE=ecoli_k12.fasta

$(HUMAN_NA12878_NATIVE_RUN1_DATA)_REFERENCE=human_g1k_v37.fasta
$(HUMAN_NA12878_NATIVE_RUN2_DATA)_REFERENCE=human_g1k_v37.fasta
$(HUMAN_NA12878_NATIVE_RUN3_DATA)_REFERENCE=human_g1k_v37.fasta
$(HUMAN_NA12878_NATIVE_MERGED_DATA)_REFERENCE=human_g1k_v37.fasta
$(HUMAN_NA12878_PCR_DATA)_REFERENCE=human_g1k_v37.fasta
$(HUMAN_NA12878_PCR_MSSSI_DATA)_REFERENCE=human_g1k_v37.fasta

# Reference file used for model training, must be the same for all training sets
TRAINING_REFERENCE=ecoli_k12.fasta
TRAINING_REGION="Chromosome:50000-3250000"
METHYLTRAIN_EXTRA_OPTS = 

TRAINED_MODEL_FOFN=ecoli_er2925.pcr_MSssI.timp.021216.alphabet_cpg.fofn

# The log-likelihood threshold required to make a call
CALL_THRESHOLD=2.5

# Convert a FAST5 file to FASTA using poretools
PORETOOLS=poretools

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
ecoli_k12.fasta:
	wget ftp://ftp.ensemblgenomes.org/pub/release-29/bacteria//fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.2.29.dna.genome.fa.gz
	gunzip Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.2.29.dna.genome.fa.gz
	mv Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.2.29.dna.genome.fa $@

human_g1k_v37.fasta:
	wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
	-gunzip human_g1k_v37.fasta.gz

#
# Rules for generating fasta files from FAST5
#

# Pass reads
%.pass.fasta: $(DATA_ROOT)/%.fast5/pass
	$(PORETOOLS) fasta --type 2D $< > $@

# Fail reads
%.fail.fasta: $(DATA_ROOT)/%.fast5/fail
	$(PORETOOLS) fasta --type 2D $< > $@

%.fasta: %.pass.fasta %.fail.fasta
	cat $^ > $@

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

# Make a file-of-filenames from a model set
%.fofn: t.006.%.model c.p1.006.%.model c.p2.006.%.model
	echo $^ | tr " " "\n" > $@

define generate-training-rules

$(eval DATASET=$(basename $1))
$(eval ALPHABET=$2)
$(eval PREFIX=$(DATASET).alphabet_$(ALPHABET))

# Training rule
$(PREFIX).fofn: $(DATASET).fasta $(DATASET).sorted.bam $(DATASET).sorted.bam.bai $(TRAINING_REFERENCE).alphabet_$(ALPHABET) ont.alphabet_$(ALPHABET).fofn
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

endef


#
# 3a. Train over a normal nucleotide (ACGT) alphabet
#

# Special case, the transformed reference for the nucleotide alphabet is just the reference
$(TRAINING_REFERENCE).alphabet_nucleotide: $(TRAINING_REFERENCE)
	ln -s $< $@

# Special case, the expanded model for the nucleotide alphabet is just the ONT model
%.alphabet_nucleotide.model: %.model
	ln -s $< $@
	
$(eval $(call generate-training-rules,$(ECOLI_K12_NATIVE_RUN1_DATA),nucleotide))
$(eval $(call generate-training-rules,$(ECOLI_K12_PCR_RUN1_DATA),nucleotide))
$(eval $(call generate-training-rules,$(ECOLI_ER2925_NATIVE_RUN1_DATA),nucleotide))
$(eval $(call generate-training-rules,$(ECOLI_ER2925_NATIVE_RUN2_DATA),nucleotide))
$(eval $(call generate-training-rules,$(ECOLI_ER2925_PCR_RUN1_DATA),nucleotide))
$(eval $(call generate-training-rules,$(ECOLI_ER2925_PCR_RUN2_DATA),nucleotide))
$(eval $(call generate-training-rules,$(ECOLI_ER2925_PCR_RUN3_DATA),nucleotide))
$(eval $(call generate-training-rules,$(ECOLI_ER2925_PCR_MSSSI_RUN1_DATA),nucleotide))
$(eval $(call generate-training-rules,$(ECOLI_ER2925_PCR_MSSSI_RUN2_DATA),nucleotide))

#
# 3b. Train over a CpG alphabet
#

# Expand the alphabet to include 5-mC k-mers
%.alphabet_cpg.model: %.model
	python $(ROOT_DIR)/expand_model_alphabet.py --alphabet cpg $< > $@

# Convert all CG dinucleotides of the reference genome to MG
$(TRAINING_REFERENCE).alphabet_cpg: $(TRAINING_REFERENCE) pythonlibs.version
	python $(ROOT_DIR)/methylate_reference.py --recognition cpg $< > $@

# train 
$(eval $(call generate-training-rules,$(ECOLI_ER2925_MSSSI_DATA),cpg))

# negative controls
$(eval $(call generate-training-rules,$(ECOLI_ER2925_NATIVE_RUN1_DATA),cpg))
$(eval $(call generate-training-rules,$(ECOLI_ER2925_NATIVE_RUN2_DATA),cpg))
$(eval $(call generate-training-rules,$(ECOLI_K12_PCR_RUN1_DATA),cpg))
$(eval $(call generate-training-rules,$(ECOLI_ER2925_PCR_RUN1_DATA),cpg))
$(eval $(call generate-training-rules,$(ECOLI_ER2925_PCR_RUN2_DATA),cpg))
$(eval $(call generate-training-rules,$(ECOLI_ER2925_PCR_MSSSI_RUN1_DATA),cpg))
$(eval $(call generate-training-rules,$(ECOLI_ER2925_PCR_MSSSI_RUN2_DATA),cpg))

#
# 3c. Train the dam model
#

# Expand the alphabet of the pcr-trained model to include 5-mC k-mers
%.alphabet_dam.model: %.model
	python $(ROOT_DIR)/expand_model_alphabet.py --alphabet dam $< > $@

# Convert all GATC sequences of the reference genome to GMTC
$(TRAINING_REFERENCE).alphabet_dam: $(TRAINING_REFERENCE) pythonlibs.version
	python $(ROOT_DIR)/methylate_reference.py --recognition dam $< > $@

$(eval $(call generate-training-rules,$(ECOLI_K12_NATIVE_DATA),dam))

#
# 3d. Train the dcm model
#

#
%.alphabet_dcm.model: %.model
	python $(ROOT_DIR)/expand_model_alphabet.py --alphabet dcm $< > $@

# Convert all dcm sequences of the reference genome
$(TRAINING_REFERENCE).alphabet_dcm: $(TRAINING_REFERENCE) pythonlibs.version
	python $(ROOT_DIR)/methylate_reference.py --recognition dcm $< > $@

$(eval $(call generate-training-rules,$(ECOLI_K12_NATIVE_DATA),dcm))

#
# 3e. Make training plots
#
training_plots_abcMG_event_mean.pdf: $(MSSSI_TRAINING_BAM).methyltrain.tsv $(PCR_TRAINING_BAM).methyltrain.tsv
	Rscript $(ROOT_DIR)/methylation_plots.R training_plots $^
	cp $@ $@.$(NOW).pdf

##################################################
#
# Step 4. Test the methylation model
#
##################################################
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
	wget http://web1.sph.emory.edu/users/hwu30/software/makeCGI/model-based-cpg-islands-hg19.txt
	cat model-based-cpg-islands-hg19.txt | python $(ROOT_DIR)/tsv_to_annotated_bed.py > $@

# Download gencode
gencode.v19.annotation.gtf:
	wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
	gunzip gencode.v19.annotation.gtf.gz

# Subset gencode to genes only
gencode.v19.annotation.genes.gtf: gencode.v19.annotation.gtf
	cat $^ | awk '$$0 ~ /#/ || ($$3 == "gene")' > $@

# Download gencode transcription start sites
gencode.v19.TSS.notlow.gff:
    wget http://genome.crg.es/~sdjebali/Gencode/version19/Fantom5_CAGE/Inputs/gencode.v19.TSS.notlow.gff

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
        python $(ROOT_DIR)/calculate_ont_signal_for_cpg_islands.py -c $(CALL_THRESHOLD) > $@

%.cpg_island_plot.pdf: NA12878.bisulfite_score.cpg_islands %.ont_score.cpg_islands
	Rscript $(ROOT_DIR)/methylation_plots.R human_cpg_island_plot $^ $@
	cp $@ $@.$(NOW).pdf
	cp histogram.pdf $*.cpg_histogram.$(NOW).pdf
	cp histogram_chromosomes.pdf $*.cpg_histogram_chromosomes.$(NOW).pdf

# Calculate methylation as a function of distance from a TSS for the ONT data
%.methylated_sites.distance_to_TSS.bed: %.sorted.bam.methyltest.sites.bed bedtools.version gencode.v19.TSS.notlow.gff
	bedtools/bin/bedtools closest -D b -b <(cat gencode.v19.TSS.notlow.gff | bedtools/bin/bedtools sort)\
                                       -a <(awk '{ print "chr" $$0 }' $*.sorted.bam.methyltest.sites.bed | bedtools/bin/bedtools sort) > $@


%.methylated_sites.distance_to_TSS.table: %.methylated_sites.distance_to_TSS.bed
	python ~/simpsonlab/users/jsimpson/code/methylation-analysis/calculate_methylation_by_distance.py --type ont -c $(CALL_THRESHOLD) -i $^ > $@

# Calculate methylation as a function of distance for the bisulfite data
bisulfite.distance_to_TSS.bed: ENCFF257GGV.bed bedtools.version gencode.v19.TSS.notlow.gff
	bedtools/bin/bedtools closest -D b -b <(cat gencode.v19.TSS.notlow.gff | bedtools/bin/bedtools sort)\
                                       -a ENCFF257GGV.bed > $@

bisulfite.distance_to_TSS.table:
	python ~/simpsonlab/users/jsimpson/code/methylation-analysis/calculate_methylation_by_distance.py --type bisulfite -i $^

##################################################
#
# Step 6. Global methylation analysis
#
##################################################
global_methylation.pdf: ProHum20kb.sorted.bam.methyltest.sites.tsv giab.NA24385.sorted.bam.methyltest.sites.tsv control.lambda.sorted.bam.methyltest.sites.tsv MSssI.lambda.sorted.bam.methyltest.sites.tsv
	Rscript $(ROOT_DIR)/methylation_plots.R global_methylation $^ $@
	cp $@ $@.$(NOW).pdf
