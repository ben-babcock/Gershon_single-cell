# save the rows of a DGE file to separate locations for fast access
# <cachedir>/expr/<region>/gene/<gene>.RDS
#
# Each row can be retrieved independently to overlay expression levels
# of a gene per cell on a t-SNE.
# Adapted from https://github.com/broadinstitute/dropviz/blob/master/prep-expr.R

library(plyr)
library(dplyr)
library(glue)

prep.dir <- '~/Dropbox (Wilhelmsen Lab)/Gershon MB Drop-Seq Shared/GeneExplorer/'

all.genes <- readRDS(glue("{prep.dir}/data/GT_Data.rds"))

expr.dir <- paste0(prep.dir,'/staged_data/GT')
xy.dir <- paste0(prep.dir,'/tsne')
N <- 100 # chunk size
expr <- all.genes
data.frame(GENE = 'a', FILE = 'b')

fh.dict <-  ldply(seq(1, nrow(expr), by=N), function(i) {
    fn <- glue("{expr.dir}/{i}.RDS")
    if (!file.exists(fn)) {
      sub.expr <- expr[i:min(nrow(expr),(i+N-1)),]
      saveRDS(as.matrix(sub.expr), fn)
      rm(sub.expr)
      gc()
    }
    fn <- strsplit(fn, "//")[[1]][2]
    data.frame(GENES=rownames(expr)[i:min(nrow(expr),(i+N-1))], FILE=fn)
  })

saveRDS(fh.dict, paste0(prep.dir,'/data/GT_expr.dict.rds'))

t <- readRDS(paste0(prep.dir,'/data/expr.dict.rds'))


