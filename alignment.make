HELL=/bin/bash -o pipefail

# do not leave failed files around
.DELETE_ON_ERROR:

# do not delete intermediate files
.SECONDARY:

# Run bwa
$(OUTPUT): $(READS) $(REFERENCE)
	$(BWA) mem -t $(THREADS) -x ont2d $(REFERENCE) $(READS) | $(SAMTOOLS) view -Sb - | $(SAMTOOLS) sort -f - $(OUTPUT)
