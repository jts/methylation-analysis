library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

chroms=paste0("chr", c("1", "2", "3", "4", "5", "6", "7", "8","9","10",
                       "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
                       "21", "22", "X", "Y", "M"))
                       

##Load gencode TSS, generated using awk script https://github.com/sdjebali/MakeGencodeTSS on v24 gencode 
gencode=read.delim(gzfile("gencode.v24.annotation_capped_sites_nr_with_confidence.gff.gz"), header=F, stringsAsFactors=F)

##Make array for BED
gencode.tss=cbind(gencode[,1], gencode[,4]-1, gencode[,5]-1, 'tss', 0, gencode[,7])

write.table(gencode.tss,gzfile("gencode.v24.tss.bed.gz"), row.names=F, col.names=F, sep="\t", quote=F)

##Make GenomicRanges Object of TSS
gencode.tss.gr=GRanges(seqnames=gencode[,1], ranges=IRanges(start=gencode[,4],end=gencode[,5]), strand=gencode[,7])


##Find promoter region, using GenomicRanges Defaults (2kb upstream, 200bp downstream, strand aware)
gencode.promoter.gr=promoters(gencode.tss.gr)

##Set chrom lengths)
seqlengths(gencode.promoter.gr)[chroms]=seqlengths(Hsapiens)[chroms]

##Trim outof bounds promoters
gencode.promoter.gr=trim(gencode.promoter.gr)

gencode.promoter=cbind(as.character(seqnames(gencode.promoter.gr)), format(start(gencode.promoter.gr)-1, scientific=FALSE),
                       format(end(gencode.promoter.gr)-1, scientific=FALSE), "promoter", 0, as.character(strand(gencode.promoter.gr)))

write.table(gencode.promoter, gzfile("gencode.v24.promoter.bed.gz"), row.names=F, col.names=F, sep="\t", quote=F)


