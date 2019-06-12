require(edgeR)

analyzePair <- function(mouseGroup,dat,old=F,qcPlots=F) {

  # First, select only the counts from the samples that belong
  # to the specific groups (i,j). Filter for CDR that have
  # a length of >= 8 and <= 16.
  #

  key <- dat$key %>% filter(MouseGrp==mouseGroup)

  if(old) {
    dd <- dat$counts %>%
      select(CDR,key$SampleID) %>%
      data.frame(check.names=F) %>%
      column_to_rownames("CDR")
  } else {
      dd <- dat$counts %>%
      select(CDR,key$SampleID) %>%
      filter(nchar(CDR)>7 & nchar(CDR)<17) %>%
      data.frame(check.names=F) %>%
      column_to_rownames("CDR")
  }

  dd=dd[,key$SampleID]

  mouse = factor(key$Mouse)
  group = factor(key$Group)

  design=model.matrix(~mouse+group)

  # Then filter for CDR's that have more than 10 counts in
  # at least (minGroup - 1) samaples, where minGroup is the
  # size of the _smaller_ group.
  #

  if(old) {
    gt10=rowSums(dd)>10
    dg=dd[gt10,]
  } else {
    minGroupNum=min(table(group))
    ii=rowSums(dd>10)>=(minGroupNum-1)
    dg=dd[ii,]
  }

  # Then process the data using R/Bioconductions edgeR package. Normalize
  # using the standard edgeR normalization fundtion and use the exactTest
  # method for computing the p-value for the pairwise differential test.
  #

  d=DGEList(counts=dg,genes=data.frame(Symbol=rownames(dg)),group=group)
  d=calcNormFactors(d)

  scaleFactor=2^mean(log2(colSums(d$counts)))
  ds=(d$counts %*% diag(scaleFactor/(d$samples$norm.factors*d$samples$lib.size)))
  colnames(ds)=colnames(d$counts)

  # dn=cpm(d)
  # scaleFactor=2^(mean(log2(colSums(dd)))-mean(log2(colSums(dn))))
  # dn=dn*scaleFactor
  # ds=dn

  pseudo=min(ds[ds>0])
  cat("pseudo=",pseudo,"\n")

  d=estimateGLMCommonDisp(d,design)
  d=estimateGLMTrendedDisp(d,design)
  d=estimateGLMTagwiseDisp(d,design)

  if(qcPlots) {
    plotMDS(d,method="bcv",labels=cc(group,mouse),main=mouseGroup)
    plotBCV(d,main=mouseGroup)
  }

  fit=glmFit(d,design)
  group.coef=grep("group",colnames(design))
  lrt=glmLRT(fit,coef=group.coef)
  tbl=topTags(lrt,n=nrow(d))

  # For the purposes of computing both average expression in a natural scale
  # rescale the normalized counts to the original sizeof the dataset. This is
  # done by setting a scaling factor that is equal to the geometric mean of
  # of the total sample counts for the samples. Additionally to compute an
  # approximate of the fold change when one of the groups has zero counts compute
  # a pseudo count by taking the smallest non-zero value of the normalized counts.
  # The compute fold changes as (A+ps)/(B+ps) where ps is this pseudo count.

  # scaleFactor=2^mean(log2(colSums(d$counts)))
  # ds=(d$counts %*% diag(scaleFactor/(d$samples$norm.factors*d$samples$lib.size)))
  # colnames(ds)=colnames(d$counts)

  dds=as.data.frame(ds) %>%
    rownames_to_column("CDR") %>%
    as_tibble %>%
    gather(SampleID,Count,-1) %>%
    left_join(key) %>%
    select(CDR,Count,Group,Mouse) %>%
    spread(Group,Count)

  grps=levels(group)
  dds$Del=log2(dds[[grps[2]]]+pseudo)-log2(dds[[grps[1]]]+pseudo)
  nLogFC=dds %>% group_by(CDR) %>% summarize(logFC.Ps=mean(Del))

  ans=tibble(
    CDR=rownames(tbl$table),
    p.value=tbl$table$PValue,
    FDR=tbl$table$FDR,
    logFC=tbl$table$logFC
    ) %>%
    left_join(nLogFC) %>%
    left_join(data.frame(ds,check.names=F) %>% rownames_to_column("CDR"))

  ans

}

volcanoPlot <- function(ans,mouseGroup,titleTxt,classes) {

  plotFDR=ifelse(ans$FDR<1e-6,1e-6,ans$FDR)

  compTag=paste(classes[2],"/",classes[1])
  n.sig=filter(ans,FDR<0.05 & abs(logFC)>log2(2)) %>% nrow
  n.filt=nrow(ans)

  subTitle=paste(n.filt,"total CDR's",n.sig,paste("significant @ FDR <",round(qcut,2)))

  par(pty='s')
  plot(ans$logFC,plotFDR,log='y',axes=F,
    pch=20,col=8,
    xlab=expression('Log'[2] * ' Fold-Change'),
    ylab="FDR",
    main=paste(compTag," --- ",titleTxt,"\n",subTitle),
    xlim=c(-17,17),
    ylim=c(1e-6,1),
    asp=1)


  axis(1,4*(-4:4),4*(-4:4))
  axis(2,10^c(-6,-4,-2,0),labels=F)
  text(y=10^c(-6,-4,-2,0),par("usr")[1],
    labels=c(
        expression(10^-6),
        expression(10^-4),
        expression(10^-2),
        expression(10^0)
        ),
      srt=0,pos=2,xpd=T)
  box()

    jj=ans$FDR<qcut
    points(ans$logFC[jj],plotFDR[jj],pch=20,col="#880000")
  abline(h=0.05,col=1,lty=2,lwd=2)



}

analyzePooled <- function(g1,g2,dat,old=F,qcPlots=F) {

  # First, select only the counts from the samples that belong
  # to the specific groups (i,j). Filter for CDR that have
  # a length of >= 8 and <= 16.
  #

  key <- dat$key %>% filter(Group %in% c(g1,g2))

  if(old) {
    dd <- dat$counts %>%
      select(CDR,key$SampleID) %>%
      data.frame(check.names=F) %>%
      column_to_rownames("CDR")
  } else {
      dd <- dat$counts %>%
      select(CDR,key$SampleID) %>%
      filter(nchar(CDR)>7 & nchar(CDR)<17) %>%
      data.frame(check.names=F) %>%
      column_to_rownames("CDR")
  }

  dd=dd[,key$SampleID]

  mouse = factor(key$Mouse)
  group = factor(key$Group)

  design=model.matrix(~group)

  # Then filter for CDR's that have more than 10 counts in
  # at least (minGroup - 1) samaples, where minGroup is the
  # size of the _smaller_ group.
  #

  if(old) {
    gt10=rowSums(dd)>10
    dg=dd[gt10,]
  } else {
    minGroupNum=min(table(group))
    ii=rowSums(dd>10)>=(minGroupNum-1)
    dg=dd[ii,]
  }

  # Then process the data using R/Bioconductions edgeR package. Normalize
  # using the standard edgeR normalization fundtion and use the exactTest
  # method for computing the p-value for the pairwise differential test.
  #

  d=DGEList(counts=dg,genes=data.frame(Symbol=rownames(dg)),group=group)
  d=calcNormFactors(d)

  scaleFactor=2^mean(log2(colSums(d$counts)))
  ds=(d$counts %*% diag(scaleFactor/(d$samples$norm.factors*d$samples$lib.size)))
  colnames(ds)=colnames(d$counts)

  # dn=cpm(d)
  # scaleFactor=2^(mean(log2(colSums(dd)))-mean(log2(colSums(dn))))
  # dn=dn*scaleFactor
  # ds=dn

  pseudo=min(ds[ds>0])
  cat("pseudo=",pseudo,"\n")

  d=estimateGLMCommonDisp(d,design)
  d=estimateGLMTrendedDisp(d,design)
  d=estimateGLMTagwiseDisp(d,design)

  if(qcPlots) {
    plotMDS(d,method="bcv",labels=cc(group,mouse),main=c(g1,g2))
    plotBCV(d,main=c(g1,g2))
  }

  fit=glmFit(d,design)
  group.coef=grep("group",colnames(design))
  lrt=glmLRT(fit,coef=group.coef)
  tbl=topTags(lrt,n=nrow(d))

  # For the purposes of computing both average expression in a natural scale
  # rescale the normalized counts to the original sizeof the dataset. This is
  # done by setting a scaling factor that is equal to the geometric mean of
  # of the total sample counts for the samples. Additionally to compute an
  # approximate of the fold change when one of the groups has zero counts compute
  # a pseudo count by taking the smallest non-zero value of the normalized counts.
  # The compute fold changes as (A+ps)/(B+ps) where ps is this pseudo count.

  # scaleFactor=2^mean(log2(colSums(d$counts)))
  # ds=(d$counts %*% diag(scaleFactor/(d$samples$norm.factors*d$samples$lib.size)))
  # colnames(ds)=colnames(d$counts)

  # dds=as.data.frame(ds) %>%
  #   rownames_to_column("CDR") %>%
  #   as_tibble %>%
  #   gather(SampleID,Count,-1) %>%
  #   left_join(key) %>%
  #   select(CDR,Count,Group,Mouse) %>%
  #   spread(Group,Count)

  # grps=levels(group)
  # dds$Del=log2(dds[[grps[2]]]+pseudo)-log2(dds[[grps[1]]]+pseudo)
  # nLogFC=dds %>% group_by(CDR) %>% summarize(logFC.Ps=mean(Del))

  ans=tibble(
    CDR=rownames(tbl$table),
    p.value=tbl$table$PValue,
    FDR=tbl$table$FDR,
    logFC=tbl$table$logFC
    ) %>%
    left_join(data.frame(ds,check.names=F) %>% rownames_to_column("CDR"))

  ans

}
