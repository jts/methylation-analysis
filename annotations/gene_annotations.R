library(GenomicRanges)

##Load gencode TSS, generated using awk script https://github.com/sdjebali/MakeGencodeTSS on v24 gencode 
gencode=read.delim(gzfile("gencode.v24.annotation_capped_sites_nr_with_confidence.gff.gz"), header=F, stringsAsFactors=F)

##Make array for BED
gencode.tss=cbind(gencode[,1], gencode[,4]-1, gencode[,5]-1, 'tss', 0, gencode[,7])

write.table(gencode.tss,gzfile("gencode.v24.tss.bed.gz"), row.names=F, col.names=F, sep="\t", quote=F)

gencode.tss.gr=GRanges(seqnames=gencode[,1], ranges=IRanges(start=gencode[,4],end=gencode[,5]), strand=gencode[,7])

##Find promoter region, using GenomicRanges Defaults (2kb upstream, 200bp downstream, strand aware)
gencode.promoter.gr=promoters(gencode.tss.gr)

gencode.promoter=cbind(as.character(seqnames(gencode.promoter.gr)), format(start(gencode.promoter.gr)-1, scientific=FALSE),
                       format(end(gencode.promoter.gr)-1, scientific=FALSE), "promoter", 0, as.character(strand(gencode.promoter.gr)))

write.table(gencode.promoter, gzfile("gencode.v24.promoter.bed.gz"), row.names=F, col.names=F, sep="\t", quote=F)


