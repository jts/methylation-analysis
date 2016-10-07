library(bsseq)
library(foreach)
library(ggplot2)
library(reshape2)

load_bed <- function(bedfile) {

    ##Load region of interest from bed file
    reg.bed=read.delim(gzfile(bedfile), sep="\t", header=F)
    reg.gr=GRanges(seqnames=reg.bed[,1], ranges=IRanges(start=reg.bed[,2]+1, end=reg.bed[,3]))

    return(reg.gr)
}


load_bsseq <- function(cytoreports) {
   
    samps=gsub("sorted.bam.methyltest.cyto.txt.*", "nanopore", basename(cytoreports))
    samps=gsub("cyto.txt.*", "bisulfite", samps)

    samps=toupper(samps)

    methobj=read.bismark(files=cytoreports, sampleNames=samps, rmZeroCov=T, fileType="cytosineReport")
    
    #methobj=foreach(i=1:length(cytoreports), .combine=combine) %do% {
    #    single=read.bismark(files=cytoreports[i], sampleNames=samps[i], rmZeroCov=T, fileType="cytosineReport")
    #    ##Bug workaround
    #    sampleNames(single)=samps[i]
    #    single
    #}

}

load_strand <- function(pfiles) {

    samps=gsub("sorted.bam.methyltest.phase.tsv.*", "nanopore", basename(pfiles))

    samps=toupper(samps)
    
    tot.raw=foreach(i=1:length(pfiles), .combine=rbind) %do% {
        
        dat.raw=read.delim(pfiles[i], header=T, stringsAsFactors=F)
        dat.raw$samp=samps[i]
        dat.raw$uid=paste(dat.raw$samp, dat.raw$readix, sep=";")

        dat.raw
    }
    
    return(tot.raw)
}

strand.overlap <- function(dat.raw, this.reg) {

    dat.raw=dat.raw[((dat.raw$start>=start(this.reg))&(dat.raw$start<=end(this.reg))),]

    return(dat.raw)
}


region.sort.filter <- function(methobj, reg, cthresh=10, diff.thresh=.2, regout="sorted.bed.gz") {
    ##Sorting and filtering based on coverage
    ##Assuming only two

    samp.cov=getCoverage(methobj, regions=reg, type="Cov", what="perRegionAverage")

    reg=reg[(rowSums(samp.cov>cthresh, na.rm=T)==dim(methobj)[2])]

    samp.meth=getMeth(methobj, regions=reg, type="raw", what="perRegion")

    diff.meth=samp.meth[,1]-samp.meth[,2]

    reg=reg[abs(diff.meth)>.2]

    reg.bedout=data.frame(chr=seqnames(reg),
                      start=start(reg)-1,
                      end=end(reg),
                      names=".",
                      scores=".",
                      strands=".")
                      
    write.table(reg.bedout, file=gzfile(regout), quote=F, sep="\t", row.names=F, col.names=F)
}

corr.plot <- function(methobj, cthresh=10, plotname="corr.pdf") {
    library(ggExtra)
    

    ##Assuming two samples only in bsseq
    methobj=methobj[,1:2]
    
    basecov=getCoverage(methobj, type="Cov", what="perBase")
    
    ##Looking for at least cthresh coverage in both samples
    enough.cov=rowSums(basecov>cthresh)>=2

    ##Get meth values for CGs that have sufficient coverage in both samples
    enough.meth=as.data.frame(getMeth(methobj[enough.cov,], type="raw", what="perBase"))

    ##Get coverage of lesser (probably nanopore) data
    enough.meth$mincov=apply(basecov[enough.cov,], 1, min)

    ##Get correlation value
    corval=cor(enough.meth[,1], enough.meth[,2])
                          
    cor.title=paste0(">", cthresh, "X: R=",format(corval, digits=2))

    ##Difference
    enough.meth$mdiff=enough.meth[,1]-enough.meth[,2]
    
    pdf(plotname, useDingbats=F)    
    vs.scat=ggplot(enough.meth, aes_string(x=colnames(enough.meth)[1], y=colnames(enough.meth)[2]))+geom_point(alpha=.4, show.legend=F)+theme_bw()+
          ggtitle(paste0("Correlation: ", cor.title))
    print(ggMarginal(vs.scat, type="histogram", margins="both"))

    print(ggplot(enough.meth, aes_string(x=colnames(enough.meth)[1], y="mdiff"))+geom_point(alpha=.4, show.legend=F)+geom_smooth(method="loess")+theme_bw()+
          ggtitle(paste0("Correlation: ", cor.title)))
          
    dev.off()
    
}

filt.cgs <- function(methobj, cthresh=5) {
    ##Filter all CGs that all samples aren't above a given threshold

        methobj=methobj[(rowSums(getCoverage(methobj, type="Cov", what="perBase")>cthresh))==dim(methobj)[2],]

        return(methobj)
}


plot.meth.reg <- function(methobj, this.reg, plotcov=F, plotdiff=F, ns=10, h=500) {
    ##Plot data for region

    idx=overlapsAny(granges(methobj), this.reg)

    this.obj=methobj[idx,]
    
    if (dim(this.obj)[1]>0) {
        this.obj=BSmooth(this.obj, ns=ns, h=h, verbose=F)
        
        coord=granges(this.obj)
        
        meth=data.frame(coord=start(coord), getMeth(this.obj, type="raw", what="perBase"))
        sm.meth=data.frame(coord=start(coord), getMeth(this.obj, type="smooth", what="perBase"))        
        
        mel.meth=melt(meth, value.name="meth", variable.name="sid", id.vars="coord")
        mel.smooth=melt(sm.meth, value.name="meth", variable.name="sid", id.vars="coord")
        
        ##Get samples and analysis type split out
        mel.meth$sample=sapply(strsplit(as.character(mel.meth$sid), split="\\."), head, 1)
        mel.meth$tech=sapply(strsplit(as.character(mel.meth$sid), split="\\."), tail, 1)

        mel.smooth$sample=sapply(strsplit(as.character(mel.smooth$sid), split="\\."), head, 1)
        mel.smooth$tech=sapply(strsplit(as.character(mel.smooth$sid), split="\\."), tail, 1)
        
        mel.meth=mel.meth[!is.na(mel.meth$meth),]
        mel.smooth=mel.smooth[!is.na(mel.smooth$meth),]
        
        print(ggplot(mel.meth, aes(x=coord, y=meth, color=sample, linetype=tech, shape=tech))+geom_point(alpha=.4)+
              geom_line(data=mel.smooth)+facet_wrap(~tech, ncol=1)+
              theme_bw()+labs(title=paste0("Methylation ", i, ": ", as.character(seqnames(this.reg[1]))))+
              xlab("Coordinate")+ylab("Percent Methylated")+
              guides(col=guide_legend(label.position="bottom", nrow=2))+theme(legend.position="bottom"))

        if (plotdiff) {
            
            diff.meth=data.frame(coord=start(coord), meth[,seq(from=2, to=dim(meth)[2], by=2)]-meth[,seq(from=3, to=dim(meth)[2], by=2)])
            diff.sm=data.frame(coord=start(coord), sm.meth[,seq(from=2, to=dim(sm.meth)[2], by=2)]-sm.meth[,seq(from=3, to=dim(sm.meth)[2], by=2)])

            colnames(diff.meth)=gsub("MCF10A", "MCF10A-MDAMB231", colnames(diff.meth))
            colnames(diff.sm)=gsub("MCF10A", "MCF10A-MDAMB231", colnames(diff.sm))
            
            mel.diff.meth=melt(diff.meth, value.name="diff.meth", variable.name="sample", id.vars="coord")
            mel.diff.smooth=melt(diff.sm, value.name="diff.meth", variable.name="sample", id.vars="coord")
            
            mel.diff.meth=mel.diff.meth[!is.na(mel.diff.meth$diff.meth),]
            mel.diff.smooth=mel.diff.smooth[!is.na(mel.diff.smooth$diff.meth),]
            
            print(ggplot(mel.diff.meth, aes(x=coord, y=diff.meth, color=sample))+geom_point(alpha=.4)+
                  geom_line(data=mel.diff.smooth)+
                  theme_bw()+labs(title=paste0("Methylation Difference ", i, ": ", as.character(seqnames(this.reg[1]))))+
                  guides(col=guide_legend(label.position="bottom", nrow=2))+theme(legend.position="bottom"))                        
        }
        
        
        if (plotcov) {
            ##Get coverage
            cov=data.frame(coord=start(coord), getCoverage(this.obj, type="Cov", what="perBase"))
            mel.cov=melt(cov, value.name="cov", variable.name="sid", id.vars="coord")

            mel.cov$sample=sapply(strsplit(as.character(mel.cov$sid), split="\\."), head, 1)
            mel.cov$tech=sapply(strsplit(as.character(mel.cov$sid), split="\\."), tail, 1)

            
            ##Get avg coverage for regions
            colMeans(cov[,2:dim(cov)[2]])
            covrep=paste(format(colMeans(cov[,2:dim(cov)[2]]), digits=2), collapse=";")
            
            print(ggplot(mel.cov, aes(x=coord, y=cov, color=sample, linetype=tech, shape=tech))+geom_point(alpha=.4)+
                  theme_bw()+labs(title=paste0("Coverage ", i, ": ", covrep))+
                  guides(col=guide_legend(label.position="bottom", nrow=2))+theme(legend.position="bottom")+
                  facet_wrap(~tech, ncol=1,scale="free_y"))
        }
        
        
    }
}



plot.strand <- function(strand.loc, this.reg, called.only=T) {
### Plot strand data; called.only shows the ones with a confident call only

    if (!called.only) {
        print(ggplot(strand.loc, aes(x=start, y=uid, color=factor(meth), group=uid))+geom_point(alpha=.5)+theme_bw()+
              facet_wrap(~samp, ncol=1, scale="free_y")+
              labs(title=paste0(i, ": ", as.character(seqnames(this.reg)))))
    }

    strand.loc=strand.loc[strand.loc$meth!=0,]

    
    print(ggplot(strand.loc, aes(x=start, y=uid, color=(meth==1), group=uid))+geom_point(alpha=.5)+theme_bw()+
          facet_wrap(~samp, ncol=1, scale="free_y")+guides(col=guide_legend(title="Is Methylated"))+
          xlim(start(this.reg), end(this.reg))+
          labs(title=paste0(i, ": ", as.character(seqnames(this.reg)))))
    
}

strand.loglik.hist <- function(strand.loc, highlight) {
    ##Plot distribution of loglik ratio per strand, highlighting specific strands

    strand.loc$col="blue"
    strand.loc$col[strand.loc$readix %in% highlight]="red"


    ##Remove uncalled CGs
    strand.loc=strand.loc[strand.loc$meth!=0,]
    
    print(ggplot(strand.loc, aes(x=loglikratio, color=col, group=readix))+theme_bw()+
          facet_wrap(~samp, ncol=1)+
          geom_density())

    print(ggplot(strand.loc, aes(x=loglikratio, color=col, group=readix))+theme_bw()+
          facet_wrap(~samp, ncol=1)+
          geom_freqpoly())

}
    


args=commandArgs(TRUE)
command=args[1]



if(! interactive()) {
    if(command == "correlation") {        
        cytoreports=args[c(-1, -length(args))]
        methobj=load_bsseq(cytoreports)
        corr.plot(methobj, cthresh=10, plotname=args[length(args)])
    } else if(command == "best_regions") {

        reg.gr=load_bed(args[2]) 

        cytoreports=args[c(-1, -2, -length(args))]
        methobj=load_bsseq(cytoreports)

        region.sort.filter(methobj, reg.gr, regout=args[length(args)])
        
    } else if(command == "meth_region_plot") {

        reg.gr=load_bed(args[2])

        cytoreports=args[c(-1, -2, -length(args))]
        methobj=load_bsseq(cytoreports)


        pdf(args[length(args)], height=8.5, width=11, useDingbats=F) 
        
        for (i in 1:length(reg.gr)) {
            plot.meth.reg(methobj, reg.gr[i], h=750, ns=40, plotcov=T)
        }

        dev.off()

    } else if(command == "meth_diff_plot") {

        reg.gr=load_bed(args[2])

        cytoreports=args[c(-1, -2, -length(args))]
        methobj=load_bsseq(cytoreports)

        pdf(args[length(args)], height=8.5, width=11, useDingbats=F)
        
        for (i in 1:length(reg.gr)) {
            plot.meth.reg(methobj, reg.gr[i], h=750, ns=40, plotcov=T, plotdiff=T)
        }

        dev.off()

        
        
    } else if(command == "strand_plot") {
        
        reg.gr=load_bed(args[2])

        pfiles=args[c(-1, -2, -length(args))]
        
        meth.strand=load_strand(pfiles)
                
        pdf(args[length(args)], height=8.5, width=11, useDingbats=F)

        for (i in 1:length(reg.gr)) {
            ##Filter region
            strand.loc=strand.overlap(meth.strand, reg.gr[i])
            ##plot filtered
            plot.strand(strand.loc, reg.gr[i])
        }

        dev.off()
              
    }
}

