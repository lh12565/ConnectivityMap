library(glue)
library(magrittr)
library(dplyr)
#library(homologene)
library(parallel)
library(ConnectivityMap)
library(memoise)
library(ggplot2)
library(ggrepel)
data("rankMatrix")
data("instances")
source('cmap_function.R')
load("data/MSigDB_enrich.RData")  #MSigDB_enrich: Precalculate msgibdb enrichment to compute specificity score. This file is sourced from https://github.com/PavlidisLab/regenerationCMap/blob/master/analysis/00.DownloadPreprocess/00.DownloadPreprocess.R.
load("data/GPL96.RData")  #gpl96

degs<-readRDS("data/DEG.RDS")  #Please replace with your differentially expressed genes, ensuring the data contains columns for gene_name, log2FoldChange, pvalue, and padj

# save the randomized matrices for quick calculation
## You can skip this step and directly load the pre-computed results: load("data/memoRandomV.RData") load("data/memoKsCalc.RData")
data("instances")
data("rankMatrix")
d = 100000
randomV = function(length,d){
    replicate(d,sample(x=1:ncol(rankMatrix),
                       size = length,replace = FALSE) %>% sort)
}
memoRandomV = memoise(randomV)
memoKsCalc = memoise(ksCalc)

instanceLengths = table(instances$cmap_name)
allVRandoms = instanceLengths %>% unique %>% sapply(function(x){  #x represents the number of experimental replicates for a drug treatment across cell lines
    print(x)
    return(memoRandomV(x,d))
})
allKsCalc<-allVRandoms %>% sapply(function(x){
    memoKsCalc(x,nrow(instances))
})
use_data(memoRandomV,overwrite = TRUE)
use_data(memoKsCalc,overwrite = TRUE)


# connectivity score, specificity score, reliability score
out_cmap<-score_function(degs,log2FoldChange=0,deg_significance='pvalue',score_significance='p')  #You can adjust different parameters according to your needs
save(out_cmap,file="out_cmap.RData")
res<-out_cmap$instanceScores[ out_cmap$instanceScores$cmap_name%in%out_cmap$hits,]
head(res)

# rank
## Filtering and sorting, can be appropriately modified according to your own data
out_cmap$chemScores$cmap_name<-rownames(out_cmap$chemScores)
comb<-out_cmap$instanceScores
comb<-merge(comb,out_cmap$chemScores,by="cmap_name")
comb2<-comb[ comb$nonNull>0.5&comb$instanceCount>1&comb$reliable ==TRUE&comb$p<0.05,]
#comb2<-comb[ comb$nonNull>0.5&comb$instanceCount>1&comb$reliable ==TRUE&comb$FDR<0.05,]
ranks<-comb2[ order(comb2$ConScore,comb2$nonNull,comb2$enrichment,decreasing=T),]
save(ranks,file="ranks.RData")

# plot
obj<-ranks[,c(1,6)]
obj<-obj[obj$ConScore!=0,]
obj$Instances<-1:nrow(obj)
obj$label <- rep("", nrow(obj))
obj$label[c(1:5)]<-obj$cmap_name[c(1:5)]
#obj$label[c(1,4,7)]<-obj$cmap_name[c(1,4,7)]  #FDR
ggplot(obj,aes(x=Instances,y=ConScore))+geom_point()+theme_bw()+
  geom_text_repel(aes(label = label), size = 3,
            max.overlaps = 10000, key_glyph = draw_key_point)+
  theme(
    #axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()
  )
ggsave('dotplot_cmap.pdf',width=6,height=5)


