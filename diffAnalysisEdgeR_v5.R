
source("loadData.R")
source("analyzePair.R")

require(openxlsx)

dat=loadData()
write_csv(dat$counts,"countTable.csv")

mouseGroups=dat$key %>% distinct(MouseGrp) %>% pull

qcut=0.05
fcCut=log2(4)

aSheets=list()

args<-commandArgs(trailing=T)
old<-as.logical(args[1])

pipeTag<-cc("CDRFilt_EdgeR_PAIRED_Method",ifelse(old,"Old","New"))

pdf(cc("volcanoPlots",pipeTag,".pdf"))

for(mg in mouseGroups) {

            cat("mouseGroup =",mg," ... ")
            sigTable=analyzePair(mg,dat,old=old)
            aSheets[[mg]]=sigTable %>% filter(FDR<qcut)
            titleTxt=paste(ifelse(old,"Old","New"),"Method")

            classes <- dat$key %>%
                  filter(MouseGrp==mg) %>%
                  arrange(Group) %>% mutate(Class=gsub("[-_]*\\d+$","",SampleID)) %>%
                  distinct(Class) %>%
                  pull(Class)
            volcanoPlot(sigTable,mg,titleTxt,classes)
            cat("\n")

        }

dev.off()

write.xlsx(aSheets,cc("diffAnalysis",pipeTag,".xlsx"))
