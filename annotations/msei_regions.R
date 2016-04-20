library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(foreach)

##Get MseI Cuts
msei.cut="TTAA"

chroms=paste0("chr", c(1:22, "X", "Y", "M"))

frags=foreach(i=chroms, .combine='c') %do% {
    m<-matchPattern(msei.cut, Hsapiens[[i]])
    f=GRanges(seqnames=i, ranges=ranges(gaps(m)))
    return(f)
}

##Get frags that are between 3500 and 6000
sel.frags=frags[((width(frags)>3500)&(width(frags)<6000))]

sel.bedout=data.frame(chr=seqnames(sel.frags),
                      start=start(sel.frags)-1,
                      end=end(sel.frags),
                      names=".",
                      scores=".",
                      strands=".")
                      
write.table(sel.bedout, file=gzfile("msifrags.bed.gz"), quote=F, sep="\t", row.names=F, col.names=F)
