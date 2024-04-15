library(GEOquery)


filePaths = getGEOSuppFiles("GSE137166")

filePaths

rownames(filePaths[1])
dir()


dir(path = "./GSE137166")
setwd("./GSE137166")

system("tar -xvf GSE137166_RAW.tar")
system("gzip -d *.gz")

setwd("..")


gds <- getGEO('GSE137166', destdir=".")

gds <-gds[[1]]
  
gds

phenoData(gds)

phenoData(gds)$title


rownames(phenoData(gds))


targets_tab <-as.matrix(pData(phenoData(gds))[,c("title")])


rownames(targets_tab) <-rownames(pData(phenoData(gds)))

targets_tab

write.table(targets_tab, file="s2572046_targets_tab_main.txt",sep="\t",col.names=T,row.names=T)
