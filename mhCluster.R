MH <- function(ani,bni) {
  ani=ani*1.0
  bni=bni*1.0
  aN=sum(ani)
  bN=sum(bni)
  da=sum(ani^2)/aN^2
  db=sum(bni^2)/bN^2
  rr=2*sum(ani*bni)/((da+db)*aN*bN)
  if(is.na(rr)) {
    return(0)
  } else {
    return(rr)
  }
}


MHij <- function(ds) {
  xx=matrix(NA,ncol=ncol(ds),nrow=ncol(ds))
  colnames(xx)=colnames(ds)
  rownames(xx)=colnames(ds)
  for(i in seq(ncol(xx))) {
    for(j in seq(ncol(xx))) {
      if(i<=j) {
        mhv=MH(ds[,i],ds[,j])
        xx[i,j]=mhv
        xx[j,i]=mhv
      }
      if(i==j) {
        xx[i,i]=1
      }
    }
  }
  return(xx)
}

sdist <- function(x){as.dist(1-x)}
hclust.avg<-function(x){hclust(x,method="average")}

doHeatMap.0 <- function(dd,label,labCol="",...) {
  rowDend=as.dendrogram(hclust.avg(sdist(MHij(dd))))

  heatmap.2(MHij(dd),main=label,trace='none',
            Rowv=rowDend,Colv=rev(rowDend),
            col=colorpanel(9,"yellow","#880000"),dendro="col",
            labCol=labCol,...)
}

doHeatMap <- function(dd,label,order=NULL,labCol="",MHFile=NULL,...) {
  ##rowDend=as.dendrogram(hclust.avg(sdist(MHij(dd))))
  if(is.null(order)) {
    order=seq(ncol(dd))
  }

  heatmap.2(MHij(dd[,order]),main=label,trace='none',
            Rowv=FALSE,Colv=FALSE,
            col=colorpanel(9,"yellow","#880000"),
            dendro="none",
            labCol=labCol[order],...)
  if(!is.null(MHFile)) {
    print(paste("saving file",MHFile))
    write.xls(round(MHij(dd[,order]),5),MHFile)
  }
}


plotMHCluster <- function(dat,mouseGrp) {

  dd=dat$counts %>% column_to_rownames("CDR") %>% select(matches(mouseGrp))
  dd=dd[rowSums(dd)>0,]

  group=dat$key %>% filter(MouseGrp==mouseGrp) %>% pull(Group) %>% factor
  mouse=dat$key %>% filter(MouseGrp==mouseGrp) %>% pull(Mouse) %>% factor

  BASE=paste0("Mouse-",mouseGrp)

  GROUPING="ORIG"

  ds=dd

  norm=exp(mean(log(colSums(ds))))
  dn=norm*sweep(ds,2,colSums(ds),"/")
  ds=dn
  dx=ds

  #halt("STOP L105")

  pdf(cc("mhCluster","iRep_190523",BASE,GROUPING,DATE(),".pdf"),height=8.5,width=8.5)

  cat("mouse =",mouseGrp,"...")

  par(cex.main=.8)

  for(count.cut in c(0,10,100)) {
    cat(count.cut,"\n")
    dat=ds
    dat=dat[rowSums(dat)>count.cut,]

    label=c(
      paste("Count Cut >",count.cut))

    #
    # Already ordered
    #
    pp=seq(ncol(dat))
    # tag=cc("mhCluster",BASE,"CUT",count.cut,DATE(),"MHTable.txt")
    # print(tag)
    doHeatMap(dat,label,pp,colnames(dat),margin=c(9/(11/8.5),9),MHFile=NULL)
    #doHeatMap(dat,label,pp,key[colnames(dat),"Group"],margin=c(9/(11/8.5),9),MHFile=tag)
  }

  dev.off()
  cat("\n")
}

plotMHClusterAll <- function(dat) {

  dd=dat$counts %>% column_to_rownames("CDR")
  dd=dd[rowSums(dd)>0,]

  BASE=paste0("ALL")

  GROUPING="ORIG"

  ds=dd

  norm=exp(mean(log(colSums(ds))))
  dn=norm*sweep(ds,2,colSums(ds),"/")
  ds=dn
  dx=ds

  #halt("STOP L105")

  pdf(cc("mhCluster","iRep_190523",BASE,GROUPING,DATE(),".pdf"),height=8.5,width=8.5)

  par(cex.main=.8)

  for(count.cut in c(0,10,100)) {
    cat(count.cut,"\n")
    dat=ds
    dat=dat[rowSums(dat)>count.cut,]

    label=c(
      paste("Count Cut >",count.cut))

    #
    # Already ordered
    #
    pp=seq(ncol(dat))
    # tag=cc("mhCluster",BASE,"CUT",count.cut,DATE(),"MHTable.txt")
    # print(tag)
    doHeatMap(dat,label,pp,colnames(dat),margin=c(9/(11/8.5),9),MHFile=NULL)
    #doHeatMap(dat,label,pp,key[colnames(dat),"Group"],margin=c(9/(11/8.5),9),MHFile=tag)
  }

  dev.off()
  cat("\n")
}

#####################################################
require(gplots)
require(stringr)
#####################################################
#halt("INCLUDE")
#####################################################
#####################################################


source("loadData.R")
dat=loadData()

mice=dat$key %>% distinct(MouseGrp) %>% pull

map(mice,~plotMHCluster(dat,.))

plotMHClusterAll(dat)


