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
	pip install Cython >> $@
	pip install numpy >> $@
	-git clone https://github.com/arq5x/poretools.git
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
	cd nanopolish; git checkout 194cc3faaf; make
	-cd nanopolish; git log | head -1 > ../$@

#
PORETOOLS=poretools

#
# Resources, references, etc
#
HUMAN_REFERENCE=GRCh38.primary_assembly.genome.fa
ECOLI_REFERENCE=ecoli_k12.fasta
BISULFITE_BED=ENCFF279HCL.bed
DATA_ROOT=../data

#
# Output targets
#
.DEFAULT_GOAL := all
all: all-nucleotide all-cpg all-results-plots all-cancer-normal

all-software: pythonlibs.version bedtools.version nanopolish.version samtools.version bwa.version
all-download: all-software $(BISULFITE_BED) $(HUMAN_REFERENCE) $(ECOLI_REFERENCE) MCF10A.cyto.txt.gz MDAMB231.cyto.txt.gz

# Function from http://stackoverflow.com/questions/6145041/makefile-filter-out-strings-containing-a-character
# Removes all strings containing the argument from the input list
FILTER_OUT = $(foreach v,$(2),$(if $(findstring $(1),$(v)),,$(v)))
# Only keeps strings containing the argument
FILTER_IN = $(foreach v,$(2),$(if $(findstring $(1),$(v)),$(v),))

# Build a complete list of NA12878/Ecoli data sets
# For E.coli this is just all of the ecoli* runs in the data directory
# For NA12878 we have to manually add in the merged datasets as well
all_NA12878=$(notdir $(patsubst %.fast5,%.fasta,$(wildcard $(DATA_ROOT)/NA12878.*.fast5))) \
             NA12878.native.merged.fasta NA12878.native.r9.merged.fasta
all_ecoli=$(notdir $(patsubst %.fast5,%.fasta,$(wildcard $(DATA_ROOT)/ecoli*.fast5)))
all_cancer=mcf10a.rr1.timp.031116.fasta mcf10a.rr2.timp.031116.fasta mcf10a.merged.fasta \
           mdamb231.rr1.timp.031316.fasta mdamb231.rr2.timp.031316.fasta mdamb231.merged.fasta

# Rules to generate all trained models
all-nucleotide: $(patsubst %.fasta,%.alphabet_nucleotide.fofn,$(all_ecoli))
all-cpg: $(patsubst %.fasta,%.alphabet_cpg.fofn,$(all_ecoli))

# Rules to generate all output reports
all-island-plots: $(addprefix results/,$(patsubst %.fasta,%.cpg_island_plot.pdf,$(all_NA12878)))

all-training-plots: results/figure.emissions.pdf results/figure.shift_by_position.pdf

all-accuracy-plots: results/accuracy_roc.pdf results/figure.accuracy_by_threshold.pdf results/accuracy_by_kmer.pdf results/figure.site_likelihood_distribution.pdf

all-TSS-plots: results/figure.methylation_by_TSS_distance.pdf results/figure.methylation_by_TSS_distance_by_chromosome.pdf 

all-results-plots: all-accuracy-plots all-island-plots all-TSS-plots all-training-plots results/all.training.tables.tex

all-cancer-normal: results/mcf10a.bsnanocorr.pdf results/mdamb231.bsnanocorr.pdf results/cn.region.plot.pdf results/cn.strand.plot.pdf

##################################################
#
# Step 1. Prepare input data
#
##################################################

#
# Define variables for each data set
#
ECOLI_K12_PCR_RUN1_DATA=ecoli_k12.pcr.loman.250915.fasta

# R7 data sets
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

# Cancer pair
HUMAN_MCF10A_RR1_DATA=mcf10a.rr1.timp.031116.fasta
HUMAN_MCF10A_RR2_DATA=mcf10a.rr2.timp.031116.fasta
HUMAN_MCF10A_MERGED_DATA=mcf10a.merged.fasta

HUMAN_MDAMB231_RR1_DATA=mdamb231.rr1.timp.031316.fasta
HUMAN_MDAMB231_RR2_DATA=mdamb231.rr2.timp.031316.fasta
HUMAN_MDAMB231_MERGED_DATA=mdamb231.merged.fasta

# R9 data sets
ECOLI_ER2925_PCR_R9_RUN1_DATA=ecoli_er2925.pcr.r9.timp.061716.fasta
ECOLI_ER2925_PCR_MSSSI_R9_RUN1_DATA=ecoli_er2925.pcr_MSssI.r9.timp.061716.fasta

HUMAN_NA12878_NATIVE_R9_RUN1_DATA=NA12878.native.r9.simpson.062416.fasta
HUMAN_NA12878_NATIVE_R9_RUN2_DATA=NA12878.native.r9.timp.062216.fasta
HUMAN_NA12878_NATIVE_R9_RUN3_DATA=NA12878.native.r9.timp.062416.fasta
HUMAN_NA12878_NATIVE_R9_MERGED_DATA=NA12878.native.r9.merged.fasta

HUMAN_NA12878_PCR_R9_DATA=NA12878.pcr.r9.timp.081016.fasta
HUMAN_NA12878_PCR_R9_MSSSI_DATA=NA12878.pcr_MSssI.r9.timp.081016.fasta

# For each data set that we use, define a variable containing its reference

# E.coli samples
$(foreach file,$(all_ecoli),$(eval $(file)_REFERENCE=$(ECOLI_REFERENCE)))

# NA12878 samples
$(foreach file,$(all_NA12878),$(eval $(file)_REFERENCE=$(HUMAN_REFERENCE)))

# Cancer samples
$(foreach file,$(all_cancer),$(eval $(file)_REFERENCE=$(HUMAN_REFERENCE)))

#
# Reference and region file used for model training
#
TRAINING_REFERENCE=$(ECOLI_REFERENCE)
TRAINING_REGION="Chromosome:50000-3250000"

# for debugging
METHYLTRAIN_EXTRA_OPTS = 

# the name of the output model that we use for the human analysis
TRAINED_MODEL_FOFN=ecoli_er2925.pcr_MSssI.timp.021216.alphabet_cpg.fofn

# The log-likelihood threshold required to make a call
CALL_THRESHOLD=2.5

# The minimum mapping quality to use a read
MIN_MAPQ=20

#
# Download data
#

# Reference genomes
$(ECOLI_REFERENCE):
	curl -L ftp://ftp.ensemblgenomes.org/pub/release-29/bacteria//fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.2.29.dna.genome.fa.gz | gunzip - > $@

$(HUMAN_REFERENCE):
	curl -L ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/GRCh38.primary_assembly.genome.fa.gz | gunzip - > $@

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

$(HUMAN_MCF10A_MERGED_DATA): mcf10a.rr1.timp.031116.fasta \
                             mcf10a.rr2.timp.031116.fasta 
	cat $^ > $@

$(HUMAN_MDAMB231_MERGED_DATA): mdamb231.rr1.timp.031316.fasta \
                               mdamb231.rr2.timp.031316.fasta 
	cat $^ > $@

$(HUMAN_NA12878_NATIVE_R9_MERGED_DATA): $(HUMAN_NA12878_NATIVE_R9_RUN1_DATA) \
                                        $(HUMAN_NA12878_NATIVE_R9_RUN2_DATA) \
                                        $(HUMAN_NA12878_NATIVE_R9_RUN3_DATA)
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
        samtools/samtools view -q $(MIN_MAPQ) -Sb - |\
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

# initialize the R7 models
t.006.ont.model: $(SCRIPT_DIR)/models/r7.3_e6_70bps_6mer_template_median68pA.model
	$(SCRIPT_DIR)/initialize_model.sh $^ template t.006 SQK006 > $@

c.p1.006.ont.model: $(SCRIPT_DIR)/models/r7.3_e6_70bps_6mer_complement_median68pA_pop1.model
	$(SCRIPT_DIR)/initialize_model.sh $^ complement.pop1 c.p1.006 SQK006 > $@

c.p2.006.ont.model: $(SCRIPT_DIR)/models/r7.3_e6_70bps_6mer_complement_median68pA_pop2.model
	$(SCRIPT_DIR)/initialize_model.sh $^ complement.pop2 c.p2.006 SQK006 > $@

# initalize R9 models
MODEL_REV=rev2
t.007.ont.model: $(SCRIPT_DIR)/models/template_median68pA.r9.$(MODEL_REV).model
	$(SCRIPT_DIR)/initialize_model.sh $^ template t.007 SQK007 > $@

c.p1.007.ont.model: $(SCRIPT_DIR)/models/complement_median68pA_pop1.r9.$(MODEL_REV).model
	$(SCRIPT_DIR)/initialize_model.sh $^ complement.pop1 c.p1.007 SQK007 > $@

c.p2.007.ont.model: $(SCRIPT_DIR)/models/complement_median68pA_pop2.r9.$(MODEL_REV).model
	$(SCRIPT_DIR)/initialize_model.sh $^ complement.pop2 c.p2.007 SQK007 > $@

# Make a 5-mer model from the input 6-mer model (used by nanopolish to calculate scalings)
%.5mer.model: %.model
	python nanopolish/scripts/dropmodel.py --type base -i $^ > $@

# Make a file-of-filenames from a model set

# R7 models
%.R7.fofn: t.006.%.model c.p1.006.%.model c.p2.006.%.model
	echo $^ | tr " " "\n" > $@

# R9 models
%.R9.fofn: t.007.%.model c.p1.007.%.model c.p2.007.%.model t.007.%.5mer.model c.p1.007.%.5mer.model c.p2.007.%.5mer.model
	echo $^ | tr " " "\n" > $@

# this code generates a training rule for a data set/alphabet pair
define generate-training-rules

$(eval DATASET=$(basename $1))
$(eval ALPHABET=$2)
$(eval PORE_VERSION=$3)
$(eval PREFIX=$(DATASET).alphabet_$(ALPHABET))

# Training rule
$(PREFIX).fofn: $(DATASET).fasta $(DATASET).sorted.bam $(DATASET).sorted.bam.bai $(TRAINING_REFERENCE).alphabet_$(ALPHABET) ont.alphabet_$(ALPHABET).$(PORE_VERSION).fofn nanopolish.version
	nanopolish/nanopolish methyltrain -t $(THREADS) $(METHYLTRAIN_EXTRA_OPTS) \
                                      --train-kmers all \
                                      --out-fofn $$@ \
                                      --out-suffix .$(PREFIX).model \
                                      -m ont.alphabet_$(ALPHABET).$(PORE_VERSION).fofn \
                                      -b $(DATASET).sorted.bam \
                                      -r $(DATASET).fasta \
                                      -g $(TRAINING_REFERENCE).alphabet_$(ALPHABET) \
                                      --filter-policy $(PORE_VERSION) \
                                      $(TRAINING_REGION)

# Rules for additional files made by methyltrain
t.006.$(PREFIX).model: $(PREFIX).fofn
c.p1.006.$(PREFIX).model: $(PREFIX).fofn
c.p2.006.$(PREFIX).model: $(PREFIX).fofn
methyltrain.$(PREFIX).model.summary: $(PREFIX).fofn
methyltrain.$(PREFIX).model.round4.events.tsv: $(PREFIX).fofn

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
all_ecoli_R7=$(call FILTER_OUT,r9,$(all_ecoli))
all_ecoli_R9=$(call FILTER_IN,r9,$(all_ecoli))
$(foreach file,$(all_ecoli_R7),$(eval $(call generate-training-rules,$(file),nucleotide,R7)))
$(foreach file,$(all_ecoli_R9),$(eval $(call generate-training-rules,$(file),nucleotide,R9)))

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
$(foreach file,$(all_ecoli_R7),$(eval $(call generate-training-rules,$(file),cpg,R7)))
$(foreach file,$(all_ecoli_R9),$(eval $(call generate-training-rules,$(file),cpg,R9)))

#
# 3c. Make training plots
#

# Emissions figure, three panels showing baseline, strong methylation signal, weak/no methylation signal
PANEL_A_SET=ecoli_er2925.pcr.timp.021216.alphabet_nucleotide
PANEL_BC_SET1=ecoli_er2925.pcr.timp.021216.alphabet_cpg
PANEL_BC_SET2=ecoli_er2925.pcr_MSssI.timp.021216.alphabet_cpg

PANEL_A_KMER="AGGTAG"
PANEL_B_KMER="AGGTMG"
PANEL_C_KMER="TMGAGT"

results/figure.emissions.pdf: methyltrain.$(PANEL_A_SET).panelA.tsv \
                              methyltrain.$(PANEL_BC_SET1).panelB.tsv \
                              methyltrain.$(PANEL_BC_SET2).panelB.tsv \
                              methyltrain.$(PANEL_BC_SET1).panelC.tsv \
                              methyltrain.$(PANEL_BC_SET2).panelC.tsv
	Rscript $(SCRIPT_DIR)/methylation_plots.R make_emissions_figure $@ $^

%.panelA.tsv: %.model.round4.events.tsv
	cat $^ | awk 'NR == 1 || ($$1 == "t.006" && $$2 == $(PANEL_A_KMER))' > $@

%.panelB.tsv: %.model.round4.events.tsv
	cat $^ | awk 'NR == 1 || ($$1 == "t.006" && $$2 == $(PANEL_B_KMER))' > $@

%.panelC.tsv: %.model.round4.events.tsv
	cat $^ | awk 'NR == 1 || ($$1 == "t.006" && $$2 == $(PANEL_C_KMER))' > $@

# Mean-shift by position figure
%.alphabet_cpg.model.summary.delta: %.alphabet_cpg.model.summary
	python $(SCRIPT_DIR)/generate_kmer_deltas.py --summary $^ --ont-fofn ont.alphabet_cpg.fofn > $@

results/figure.shift_by_position.pdf: methyltrain.ecoli_er2925.pcr_MSssI.timp.021216.alphabet_cpg.model.summary.delta
	mkdir -p results
	Rscript $(SCRIPT_DIR)/methylation_plots.R make_mean_shift_by_position_figure $@ $^

#
# 3d. Make training tables from the summary files
#

# PCR data

# nucleotide
results/pcr.nucleotide.training.table.tex: methyltrain.$(ECOLI_ER2925_PCR_RUN1_DATA:.fasta=.alphabet_nucleotide.model.summary) \
                                           methyltrain.$(ECOLI_ER2925_PCR_RUN2_DATA:.fasta=.alphabet_nucleotide.model.summary) \
                                           methyltrain.$(ECOLI_K12_PCR_RUN1_DATA:.fasta=.alphabet_nucleotide.model.summary) \
                                           methyltrain.$(ECOLI_ER2925_PCR_R9_RUN1_DATA:.fasta=.alphabet_nucleotide.model.summary)
	mkdir -p results
	python $(SCRIPT_DIR)/generate_training_table.py --treatment pcr --alphabet nucleotide $^ > $@

# cpg
results/pcr.cpg.training.table.tex: methyltrain.$(ECOLI_ER2925_PCR_RUN1_DATA:.fasta=.alphabet_cpg.model.summary) \
                                           methyltrain.$(ECOLI_ER2925_PCR_RUN2_DATA:.fasta=.alphabet_cpg.model.summary) \
                                           methyltrain.$(ECOLI_K12_PCR_RUN1_DATA:.fasta=.alphabet_cpg.model.summary) \
                                           methyltrain.$(ECOLI_ER2925_PCR_R9_RUN1_DATA:.fasta=.alphabet_nucleotide.model.summary)
	mkdir -p results
	python $(SCRIPT_DIR)/generate_training_table.py --treatment pcr --alphabet cpg $^ > $@

# PCR + M.SssI

# nucleotide
results/pcr_MSssI.nucleotide.training.table.tex: methyltrain.$(ECOLI_ER2925_PCR_MSSSI_RUN1_DATA:.fasta=.alphabet_nucleotide.model.summary) \
                                                 methyltrain.$(ECOLI_ER2925_PCR_MSSSI_RUN2_DATA:.fasta=.alphabet_nucleotide.model.summary) \
                                                 methyltrain.$(ECOLI_ER2925_PCR_MSSSI_R9_RUN1_DATA:.fasta=.alphabet_nucleotide.model.summary)
	mkdir -p results
	python $(SCRIPT_DIR)/generate_training_table.py --treatment pcr_MSssI --alphabet nucleotide $^ > $@

# CpG
results/pcr_MSssI.cpg.training.table.tex: methyltrain.$(ECOLI_ER2925_PCR_MSSSI_RUN1_DATA:.fasta=.alphabet_cpg.model.summary) \
                                          methyltrain.$(ECOLI_ER2925_PCR_MSSSI_RUN2_DATA:.fasta=.alphabet_cpg.model.summary) \
                                          methyltrain.$(ECOLI_ER2925_PCR_MSSSI_R9_RUN1_DATA:.fasta=.alphabet_nucleotide.model.summary)
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
$(BISULFITE_BED):
	curl -L https://www.encodeproject.org/files/ENCFF279HCL/download/$(BISULFITE_BED).gz | gunzip - > $@

# Extract tsv for strand-based methylation plots
%.methyltest.phase.tsv: %.methyltest.sites.bed
	python $(SCRIPT_DIR)/phase_extract.py -c  $(CALL_THRESHOLD) -i $<

# Summarize methylation per cg
%.methyltest.cyto.txt: %.methyltest.phase.tsv
	python $(SCRIPT_DIR)/per_cg_methylation.py -i $<

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

results/%.cpg_island_plot.pdf: NA12878.bisulfite_score.cpg_islands %.ont_score.cpg_islands
	mkdir -p results
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

# Calculate methylation as a function of distance for the bisulfite data
%.bisulfite.distance_to_TSS.bed: %.bisulfite.bed bedtools.version $(TSS_FILE)
	bedtools/bin/bedtools closest -D b -b <(zcat $(TSS_FILE) | bedtools/bin/bedtools sort)\
                                       -a <(cat $< | bedtools/bin/bedtools sort) > $@

%.bisulfite.distance_to_TSS.table: %.bisulfite.distance_to_TSS.bed
	python $(SCRIPT_DIR)/calculate_methylation_by_distance.py --type bisulfite -i $^ > $@


results/figure.methylation_by_TSS_distance.pdf: NA12878.bisulfite.distance_to_TSS.table \
                                                NA12878.native.merged.methylated_sites.distance_to_TSS.table \
                                                NA12878.pcr.simpson.021616.methylated_sites.distance_to_TSS.table \
                                                NA12878.pcr_MSssI.simpson.021016.methylated_sites.distance_to_TSS.table \
                                                NA12878.pcr.r9.timp.081016.methylated_sites.distance_to_TSS.table \
                                                NA12878.native.r9.merged.methylated_sites.distance_to_TSS.table
	mkdir -p results
	Rscript $(SCRIPT_DIR)/methylation_plots.R TSS_distance_plot $^ $@


results/figure.methylation_by_TSS_distance_by_chromosome.pdf: NA12878.bisulfite.distance_to_TSS.table \
                                                              NA12878.native.merged.methylated_sites.distance_to_TSS.table \
                                                              NA12878.pcr.r9.timp.081016.methylated_sites.distance_to_TSS.table \
                                                              NA12878.pcr_MSssI.r9.timp.081016.methylated_sites.distance_to_TSS.table
	mkdir -p results
	Rscript $(SCRIPT_DIR)/methylation_plots.R TSS_distance_plot_by_chromosome $^ $@

#
# Site likelihood plot
#
results/figure.site_likelihood_distribution.pdf: NA12878.pcr.simpson.021616.sorted.bam.methyltest.sites.tsv \
                                                 NA12878.pcr_MSssI.simpson.021016.sorted.bam.methyltest.sites.tsv \
                                                 NA12878.native.merged.sorted.bam.methyltest.sites.tsv
	mkdir -p results/
	Rscript $(SCRIPT_DIR)/methylation_plots.R site_likelihood_distribution $^ $@

#
# Accuracy results
#
accuracy.by_threshold.tsv: NA12878.pcr.simpson.021616.sorted.bam.methyltest.sites.bed \
                           NA12878.pcr_MSssI.simpson.021016.sorted.bam.methyltest.sites.bed
	python $(SCRIPT_DIR)/calculate_call_accuracy.py --unmethylated NA12878.pcr.simpson.021616.sorted.bam.methyltest.sites.bed \
                                                    --methylated NA12878.pcr_MSssI.simpson.021016.sorted.bam.methyltest.sites.bed

accuracy.precision_recall.tsv: accuracy.by_threshold.tsv
accuracy.by_kmer.tsv: accuracy.by_threshold.tsv

results/accuracy_roc.pdf: accuracy.precision_recall.tsv
	mkdir -p results
	Rscript $(SCRIPT_DIR)/methylation_plots.R call_accuracy_roc $^ $@

results/figure.accuracy_by_threshold.pdf: accuracy.by_threshold.tsv
	mkdir -p results
	Rscript $(SCRIPT_DIR)/methylation_plots.R call_accuracy_by_threshold $^ $@

results/accuracy_by_kmer.pdf: accuracy.by_kmer.tsv
	mkdir -p results
	Rscript $(SCRIPT_DIR)/methylation_plots.R call_accuracy_by_kmer $^ $@

#
# Cancer/Normal results
#

# D/l illumina bs data

MCF10A.cyto.txt.gz:
	wget http://timplab.org/data/MCF10A.cyto.txt.gz

MDAMB231.cyto.txt.gz:
	wget http://timplab.org/data/MDAMB231.cyto.txt.gz


# Correlation plot

results/mcf10a.bsnanocorr.pdf: MCF10A.cyto.txt.gz mcf10a.merged.sorted.bam.methyltest.cyto.txt
	Rscript $(SCRIPT_DIR)/methylation_region_plot.R correlation $^ $@

results/mdamb231.bsnanocorr.pdf: MDAMB231.cyto.txt.gz mdamb231.merged.sorted.bam.methyltest.cyto.txt
	Rscript $(SCRIPT_DIR)/methylation_region_plot.R correlation $^ $@

# Region Plot

# Thresholds for filtering regions to plot
COV_THRESHOLD=5
SIZE_MN_THRESHOLD=3000
SIZE_MX_THRESHOLD=7000

filt.isl.bed.gz: $(CGI_FILE) \
                 MCF10A.cyto.txt.gz \
                 MDAMB231.cyto.txt.gz
	Rscript $(SCRIPT_DIR)/methylation_region_plot.R best_regions $^ $@

covd.regions.bed: bedtools.version \
                  mcf10a.merged.sorted.bam \
                  mdamb231.merged.sorted.bam \
                  filt.isl.bed.gz 
	bedtools/bin/bedtools genomecov -ibam mcf10a.merged.sorted.bam -bg | awk -v cov=${COV_THRESHOLD} '$$4 > cov' | bedtools/bin/bedtools merge -d 50 -i stdin >mcf10a.merged.sorted.bam.covd
	bedtools/bin/bedtools genomecov -ibam mdamb231.merged.sorted.bam -bg | awk -v cov=${COV_THRESHOLD} '$$4 > cov' | bedtools/bin/bedtools merge -d 50 -i stdin >mdamb231.merged.sorted.bam.covd
	bedtools/bin/bedtools intersect -a mcf10a.merged.sorted.bam.covd -b mdamb231.merged.sorted.bam.covd -u | \
	bedtools/bin/bedtools intersect -a stdin -b filt.isl.bed.gz -u | awk -v mn=${SIZE_MN_THRESHOLD} -v mx=${SIZE_MX_THRESHOLD} '($$3 - $$2) > mn && ($$3 - $$2) < mx' | \
	bedtools/bin/bedtools coverage -a stdin -b mcf10a.merged.sorted.bam mdamb231.merged.sorted.bam | sort -k 4rn >$@


results/cn.region.plot.pdf: covd.regions.bed \
	                    mcf10a.merged.sorted.bam.methyltest.cyto.txt \
                            mdamb231.merged.sorted.bam.methyltest.cyto.txt \
                            MCF10A.cyto.txt.gz \
                            MDAMB231.cyto.txt.gz
	mkdir -p results
	Rscript $(SCRIPT_DIR)/methylation_region_plot.R meth_diff_plot $^ $@

# Strand plot
results/cn.strand.plot.pdf: covd.regions.bed \
	                    mcf10a.merged.sorted.bam.methyltest.phase.tsv \
                            mdamb231.merged.sorted.bam.methyltest.phase.tsv
	mkdir -p results
	Rscript $(SCRIPT_DIR)/methylation_region_plot.R strand_plot $^ $@
