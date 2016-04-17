library(bsseq)
library(foreach)
library(ggplot2)
library(reshape2)

#chroms=paste0("chr", c(1:22, "X", "Y", "M"))

load_bed <- function(bedfile) {

    ##Load region of interest from bed file
    reg.bed=read.delim(gzfile(bedfile), sep="\t", header=F)
    reg.gr=GRanges(seqnames=reg.bed[,1], ranges=IRanges(start=reg.bed[,2]+1, end=reg.bed[,3]))

    return(reg.gr)
}


load_bsseq <- function(cytoreports) {
   
    samps=gsub(".cyto.txt.*", "", basename(cytoreports))
    
    methobj=foreach(i=1:length(cytoreports), .combine=combine) %do% {
        single=read.bismark(files=cytoreports[i], sampleNames=samps[i], rmZeroCov=T, fileType="cytosineReport")
        ##Bug workaround
        sampleNames(single)=samps[i]
        single
    }

}

load_strand <- function(pfiles) {

    samps=gsub(".phase.tsv.*", "", basename(pfiles))

    tot.raw=foreach(i=1:length(pfiles), .combine=rbind) %do% {
        
        dat.raw=read.delim(pfiles[i], header=T, stringsAsFactors=F)
        dat.raw$samp=samps[i]
        dat.raw$uid=paste(dat.raw$samp, dat.raw$readix, sep=".")

        dat.raw
    }
    
    return(tot.raw)
}

strand.overlap <- function(dat.raw, this.reg) {

    dat.raw=dat.raw[((dat.raw$start>=start(this.reg))&(dat.raw$start<=end(this.reg))),]

    return(dat.raw)
}



region.sort.filter <- function(cytoreports, reg, cthresh=5, regout="sorted.bed.gz") {
    ##Sorting and filtering based on coverage

    samp.cov=getCoverage(methobj, regions=reg, type="Cov", what="perRegionAverage")

    ##Sort on the first sample in order
    idx=order(-rowSums(samp.cov))
    idx=idx[(rowSums(samp.cov[idx,]>cthresh)==dim(methobj)[2])]
    idx=idx[!is.na(idx)]

    reg=reg[idx]

    reg.bedout=data.frame(chr=seqnames(reg),
                      start=start(reg)-1,
                      end=end(reg),
                      names=".",
                      scores=".",
                      strands=".")
                      
    write.table(reg.bedout, file=gzfile(regout), quote=F, sep="\t", row.names=F, col.names=F)
}


corr.plot <- function(methobj, cthresh=5, plotname="corr.pdf") {
    ##Add colors for level of coverage

    ##Assuming two samples only in bsseq
    methobj=methobj[,1:2]
    
    ##Looking for at least cthresh coverage in both samples
    enough.cov=rowSums(getCoverage(methobj, type="Cov", what="perBase")>cthresh)>=2

    ##Get meth values for CGs that have sufficient coverage in both samples
    enough.meth=as.data.frame(getMeth(methobj[enough.cov,], type="raw", what="perBase"))
    
    ##Get correlation value
    corval=cor(enough.meth[,1], enough.meth[,2])
    
    pdf(plotname)
    print(ggplot(enough.meth, aes_string(x=colnames(enough.meth)[1], y=colnames(enough.meth)[2]))+geom_point(alpha=.4)+theme_bw()+ggtitle(paste0("Correlation: ", format(corval, digits=2))))
    dev.off()

}

filt.cgs <- function(methobj, cthresh=5) {
##Filter all CGs that all samples aren't above a given threshold

    methobj=methobj[(rowSums(getCoverage(methobj, type="Cov", what="perBase")>cthresh))==dim(methobj)[2],]

    return(methobj)
}



plot.meth.reg <- function(methobj, this.reg) {
    ##Plot data for region
    
    coord=granges(methobj)[overlapsAny(granges(methobj), this.reg)]
    
    meth=getMeth(methobj, regions=this.reg, type="raw", what="perBase")[[1]]
    
    ##Check for empty region
    if (is.array(meth)) {
        cov=getCoverage(methobj, regions=this.reg, type="Cov", what="perBase")[[1]]
        
        mel.meth=melt(meth, value.name="meth")
        mel.cov=melt(cov, value.name="cov")
        mel.meth$coord=start(coord)
        mel.cov$coord=start(coord)
        
        print(ggplot(mel.meth, aes(x=coord, y=meth, color=Var2))+geom_smooth(se=F)+geom_point(alpha=.4)+theme_bw()+labs(title=paste0(i, ": ", as.character(seqnames(this.reg[1])))))
        
        print(ggplot(mel.cov, aes(x=coord, y=cov, color=Var2))+geom_smooth(se=F)+geom_point(alpha=.4)+theme_bw()+labs(title=paste0(i, ": ", as.character(seqnames(this.reg[1])))))
        
    }
}

plot.strand <- function(strand.loc, this.reg, called.only=T) {
### Plot strand data

    if (!called.only) {
        print(ggplot(strand.loc, aes(x=start, y=uid, color=factor(meth), group=uid))+geom_point(alpha=.5)+theme_bw()+
              facet_wrap(~samp, ncol=1, scale="free_y")+
              labs(title=paste0(i, ": ", as.character(seqnames(this.reg)))))
    }

    strand.loc=strand.loc[strand.loc$meth!=0,]
    
    print(ggplot(strand.loc, aes(x=start, y=uid, color=factor(meth), group=uid))+geom_point(alpha=.5)+theme_bw()+
          facet_wrap(~samp, ncol=1, scale="free_y")+
          labs(title=paste0(i, ": ", as.character(seqnames(this.reg)))))
    
    
}


##args=c("correlation", "/mithril/Data/NGS/Aligned/160405_rrbs/combined/MCF10A/MCF10A.cyto.txt.gz",
#       "~/jared/methylation-run/mcf10a.merged.sorted.bam.methyltest.cyto.txt","~/Dropbox/Temp/try.pdf")

##args=c("best_regions", "annotations/msifrags.bed.gz", "~/jared/methylation-run/mcf10a.merged.sorted.bam.methyltest.cyto.txt",
#       "~/jared/methylation-run/mdamb231.merged.sorted.bam.methyltest.cyto.txt", "msei_filt.bed.gz")

#args=c("meth_region_plot", "msei_filt.bed.gz", "~/jared/methylation-run/mcf10a.merged.sorted.bam.methyltest.cyto.txt",
#       "/mithril/Data/NGS/Aligned/160405_rrbs/combined/MCF10A/MCF10A.cyto.txt.gz",
#       "/mithril/Data/NGS/Aligned/160405_rrbs/combined/MDAMB231/MDAMB231.cyto.txt.gz", 
#       "~/jared/methylation-run/mdamb231.merged.sorted.bam.methyltest.cyto.txt", "~/Dropbox/Temp/reg_try2.pdf")

#args=c("strand_plot", "msei_filt.bed.gz", "~/jared/methylation-run/mcf10a.merged.sorted.bam.methyltest.phase.tsv",
#       "~/jared/methylation-run/mdamb231.merged.sorted.bam.methyltest.phase.tsv", "~/Dropbox/Temp/strand_try2.pdf")



args=commandArgs(TRUE)
command=args[1]



if(! interactive()) {
    if(command == "correlation") {        
        cytoreports=args[c(-1, -length(args))]
        methobj=load_bsseq(cytoreports)
        corr.plot(methobj, cthresh=5, plotname=args[length(args)])
    } else if(command == "best_regions") {

        reg.gr=load_bed(args[2]) 
        
        cytoreports=args[c(-1, -2, -length(args))]
        methobj=load_bsseq(cytoreports)
        region.sort.filter(cytoreports, reg.gr, cthresh=5, regout=args[length(args)])
        
    } else if(command == "meth_region_plot") {

        reg.gr=load_bed(args[2])

        cytoreports=args[c(-1, -2, -length(args))]
        methobj=load_bsseq(cytoreports)

        methobj=filt.cgs(methobj)
        
        pdf(args[length(args)])

        for (i in 1:length(reg.gr)) {
            plot.meth.reg(methobj, reg.gr[i]) 
        }

        dev.off()
 
    } else if(command == "strand_plot") {
        
        reg.gr=load_bed(args[2])

        pfiles=args[c(-1, -2, -length(args))]
        
        meth.strand=load_strand(pfiles)
                
        pdf(args[length(args)])

        for (i in 1:length(reg.gr)) {
            ##Filter region
            strand.loc=strand.overlap(meth.strand, reg.gr[i])
            ##plot filtered
            plot.strand(strand.loc, reg.gr[i], called.only=F)
        }

        dev.off()
              
    }
}


