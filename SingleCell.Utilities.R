library(Seurat)
library(xlsx)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(sva)
library(viridis)
library(wesanderson)
library(svglite)
library(RSvgDevice)
library(caret)
library(data.table)
fonts <- list("Liberation Sans")

colors.use <- c(brewer.pal(8, 'Set1'), brewer.pal(7, 'Dark2'),  'blue4', 'orangered', 'orangered4')
more.colors <- c(brewer.pal(12, 'Paired'), brewer.pal(6, 'Accent'))
more.colors <- c(more.colors[1:10], more.colors[12], more.colors[15:20])

runComBat <- function(object, grouping.var){
  require(sva)
  genes.use <- names(which(Matrix::rowSums(object@raw.data) > 50))
  object@raw.data <- object@raw.data[genes.use, ]
  meta.data <- object@meta.data[colnames(object@raw.data), c('nUMI', 'nGene', 'Replicate')]
  edata <- as.matrix(object@raw.data)
  batch <- meta.data[, paste(grouping.var)]
  modcombat <- model.matrix(~1, data = meta.data)
  correct.data <- ComBat(dat = edata, 
                         batch = batch, 
                         mod = modcombat, 
                         prior.plots = FALSE,
                         par.prior = TRUE)
  return(correct.data)
}

compareTreatments <- function(object, group.test, group.name, levels.name){
  message('Beginning testing on: ', group.name, '\n\tLevel: ', group.test, '\n')
  object <- SetAllIdent(object, id = group.name)
  object <- SubsetData(object, ident.use = group.test)
  object <- SetAllIdent(object, id = paste(levels.name))
  levels.test <- levels(object@ident)
  dat <- FindMarkers(object, 
                     ident.1 = levels.test[2],
                     ident.2 = levels.test[1],
                     test.use = "bimod")
  colnames(dat) <- c('p_val', 
                     'avg_logFC', 
                     paste0('pct.', levels.test[2]), 
                     paste0('pct.', levels.test[1]), 
                     'p_val_adj')
  message('Writing OutFile\n')
  write.xlsx(x = dat, 
             file = paste0(group.name, '_', levels.name, '.xlsx'), 
             sheetName = paste(group.test), 
             append = T)
}

TestPCA <- function(object = NULL, 
                    genes.use = object@var.genes, 
                    use.imputed = F,
                    do.return = F,
                    do.rev = F,
                    custom.use = NULL){
  
  if(use.imputed){data.use <- object@imputed[genes.use, ]
  }else if(! is.null(custom.use)){
    data.use <- as.matrix(custom.use)
    }else{
    data.use <- object@scale.data[genes.use, ]
  }

  if(do.rev){
  pca.results <- svd(x = data.use)
  cell.embeddings <- pca.results$v
  gene.loadings <- pca.results$u %*% diag(pca.results$d)
  sdev <- pca.results$d/sqrt(max(1, nrow(data.use) - 1))
  }else{
  pca.results <- svd(x = t(data.use))
  cell.embeddings <- pca.results$u %*% diag(pca.results$d)
  gene.loadings <- pca.results$v
  sdev <- pca.results$d/sqrt(max(1, ncol(data.use) - 1))
  }
  PCVariance <- rbind(SD = sdev,
                      Proportion = (sdev^2)/sum(sdev^2), 
                      Cumulative = cumsum(sdev^2)/sum(sdev^2))
  m <- mean(PCVariance['Proportion', ])
  s <- sd(PCVariance['Proportion', ])
  PCVariance <- rbind(PCVariance, 
                      ZScore = (PCVariance['Proportion', ] - m)/s)
  print(PCVariance[, 1:20])
  if(do.return){
    return(PCVariance)
  }
}

MultiFeaturePlot <- function(Object = S2, 
                             gene.list = 'SmoM2-EYFP', 
                             colors = colors.use,
                             null.color = grey(0.7),
                             shape = 20,
                             alpha.colors = 0.9,
                             threshold = 0.0125,
                             title = NULL,
                             subtitle = NULL,
                             plot.intersection = F,
                             show.legend = T,
                             shape.legend = guide_legend('Is None'),
                             legend.title = 'Legend',
                             multi.label = 'Both',
                             pct.cutoff = 0.3,
                             pt.size = 0.5,
                             scale.color = F,
                             recover.raw = F){
  require(Seurat)

  #Get Plotting Coordinates
  frame <- GetDimReduction(Object, reduction.type = 'tsne', slot = "cell.embeddings")
  frame <- cbind(frame,
                 as.data.frame(matrix(data = 0, 
                                      nrow = nrow(frame),
                                      ncol = length(gene.list) + 3)))
  colnames(frame) <- c('tSNE_1', 'tSNE_2', gene.list, 'gene.sums', 'Plot.Status', 'is.null')

  #Should allow for multiple thresholds, not working*
  #*Fixed-ish. Requires multi-threshold input as named vector
  if(recover.raw){
    for(gene in gene.list){
      frame[colnames(Object@misc), gene] <- as.numeric(Object@misc[gene, ])
    }
  }else{
    for(gene in gene.list){
      frame[colnames(Object@data), gene] <- as.numeric(Object@data[gene, ])
      }
  }

  
  if(length(threshold) == 1){
    thresh <- threshold
    if(scale.color){
      frame[which(frame[, gene] >= thresh), 'Plot.Status'] <- paste0(gene)
      frame[which(frame[, 'Plot.Status'] == 0), 'Plot.Status'] <- 'None'
      frame[, 'Plot.Status'] <- factor(frame[, 'Plot.Status'], levels = c(gene.list, 'None'))
      frame1 <- frame[which(frame$Plot.Status == 'None'), ]
      frame2 <- frame[which(frame$Plot.Status != 'None'), ]
     d <- ggplot(mapping = aes(x = tSNE_1, y = tSNE_2)) +
        geom_point(data = frame1, colour = null.color, size = pt.size, shape = shape) +
        geom_point(data = frame2, aes_string(colour = gene), size = pt.size, shape = shape, show.legend = show.legend) +
        scale_colour_gradient(low = null.color, high = colors, aesthetics = 'colour', guide = 'colourbar') +
        ggtitle(title, subtitle = subtitle)
     return(d)
    }
  }
stopifnot(scale.color == F)

  
  #Plot Intersections: Cells that have multiple of genes vector
  if(plot.intersection){#Color by "Plot.status"
    for(gene in gene.list){
      if(length(threshold) > 1){
        thresh <- threshold[gene]
      }
      frame[which(frame[, gene] <= thresh), gene] <- paste0(0)
      frame[which(frame[, gene] > thresh), gene] <- paste0(1)
    }
    multis <- which(rowSums(apply(frame[, gene.list], 2, as.numeric)) >= (length(gene.list)*pct.cutoff))
    frame[multis, 'Plot.Status'] <- multi.label
    frame[which(frame[, 'Plot.Status'] != multi.label), 'Plot.Status'] <- 'None'
    frame[, 'Plot.Status'] <- factor(frame[, 'Plot.Status'], levels = c(multi.label, 'None'))
    frame1 <- frame[which(frame$Plot.Status == 'None'), ]
    frame2 <- frame[which(frame$Plot.Status != 'None'), ]
    
    d <- ggplot(mapping = aes(x = tSNE_1, y = tSNE_2)) +
      geom_point(data = frame1, color = null.color, size = pt.size, shape = shape) +
      geom_point(data = frame2, aes(color = Plot.Status), size = pt.size, shape = shape, show.legend = show.legend) +
      scale_colour_manual(values = colors) +
      ggtitle(title, subtitle = subtitle) +
      guides(colour = guide_legend(override.aes = list(size=10), title = legend.title))
    return(d)
    stopifnot(plot.intersection == F)
  }else{#Color by gene
    for(gene in gene.list){
      if(length(threshold) > 1){
        thresh <- threshold[gene]
      }else{
        thresh <- threshold}
      frame[which(frame[, gene] >= thresh), 'Plot.Status'] <- paste0(gene)
      frame[which(frame[, 'Plot.Status'] == 0), 'Plot.Status'] <- 'None'
    }
    frame[, 'Plot.Status'] <- factor(frame[, 'Plot.Status'], levels = c(gene.list, 'None'))
    frame1 <- frame[which(frame$Plot.Status == 'None'), ]
    frame2 <- frame[which(frame$Plot.Status != 'None'), ]
    #frame[which(frame[, 'Plot.Status'] == 'None'), 'is.null'] <- '1'
    
    d <- ggplot(mapping = aes(x = tSNE_1, y = tSNE_2)) +
      geom_point(data = frame1, color = null.color, size = pt.size, shape = shape) +
      geom_point(data = frame2, aes(color = Plot.Status), size = pt.size, shape = shape, show.legend = show.legend) +
      scale_colour_manual(values = colors) +
      ggtitle(title, subtitle = subtitle) +
      guides(colour = guide_legend(override.aes = list(size=10), title = legend.title))
    }
  return(d)
}

getPos <- function(mat, genes.use, c){
  tst <- sweep(x = mat[genes.use, c], MARGIN = 1, 
               STATS = threshold,
               FUN = '>')
  tst2 <- round(rowSums(tst)/length(x = c), 3)
  tst2
}

nPositiveCells <- function(object = NULL,
                           ident.test = NULL,
                           ident.compare = NULL,
                           replicate.name = 'Replicate',
                           threshold = 0, 
                           genes.use = NULL,
                           ident.name = NULL,
                           display.progress = T){
  require(Seurat)
  if(is.null(object) || is.null(ident.test)){
    stop('Please enter an object and primary identity')
  }
  if(is.null(genes.use)){
    genes.use <- rownames(object@data)
    message('Setting gene list to all genes')
  }
  if(is.null(ident.compare)){
    ident.compare <- levels(object@ident)[!(levels(object@ident) %in% ident.test)]
    message('Setting comparison identity to all other cells')
  }
  getPos <- function(mat, genes.use, c){
    tst <- sweep(x = mat[genes.use, c], MARGIN = 1, 
                 STATS = threshold,
                 FUN = '>')
    tst2 <- round(rowSums(tst)/length(x = c), 3)
    tst2
  }
  
  genes.use <- intersect(genes.use, rownames(object@data))
  cells <- colnames(object@data)
  replicate.levels <- unique(object@meta.data[, replicate.name])
  
  if (display.progress) {
    message("Finding Positive Cells")
    pb <- txtProgressBar(min = 0, max = length(replicate.levels), style = 3, file = stderr())
  }
  mtx.test <- matrix(data = genes.use, nrow = length(genes.use), ncol = 1)
  mtx.comp <- matrix(data = genes.use, nrow = length(genes.use), ncol = 1)
  mtx.out <- matrix(data = c('Gene', replicate.levels, 'Identity'), nrow = 1, ncol = length(replicate.levels) + 2)
  for(rep.name in replicate.levels){
    cells.test <- cells[which((object@meta.data[, ident.name] %in% ident.test) 
                              & (object@meta.data[, replicate.name] == rep.name))]
    cells.compare <- cells[which((object@meta.data[, ident.name] %in% ident.compare) 
                                 & (object@meta.data[, replicate.name] == rep.name))]
    mtx.test2 <- as.data.frame(getPos(object@data, genes.use = genes.use, c = cells.test))
    mtx.test <- cbind(mtx.test, mtx.test2[, 1])
    
    mtx.comp2 <- as.data.frame(getPos(object@data, genes.use = genes.use, c = cells.compare))
    mtx.comp <- cbind(mtx.comp, mtx.comp2[, 1])
    
    if (display.progress) {
      setTxtProgressBar(pb = pb, value = which(replicate.levels == rep.name))
    }
  }
  if (display.progress){
    close(con = pb)
  }
  mtx.test <- cbind(mtx.test, rep(ident.test, length(genes.use)))
  mtx.comp <- cbind(mtx.comp, rep('ident.compare', length(genes.use)))
  mtx.out <- rbind(mtx.out, mtx.test, mtx.comp)
  colnames(mtx.out) <- mtx.out[1,]
  mtx.out <- mtx.out[-1,]
  mtx.out <- as.data.frame(mtx.out)
  mtx.out[, 2:(ncol(mtx.out)-1)] <- apply(mtx.out[, 2:(ncol(mtx.out)-1)], 2, function(x) as.numeric(x))
  mtx.out[, 1] <- as.character(mtx.out[, 1])
  mtx.out[, ncol(mtx.out)] <- as.character(mtx.out[, ncol(mtx.out)])
  mtx.out$p.val <- rep(1, nrow(mtx.out))
  genes.use <- unique(mtx.out[,1])
  if (display.progress) {
    message("Testing Positive Values")
    pb <- txtProgressBar(min = 0, max = length(genes.use),
                         style = 3, file = stderr())
  }
  for(gene in genes.use){
    tryCatch({do.test <- mtx.out[which(mtx.out[, 1] == gene), ]
    p.val <- t.test(x = do.test[1, 2:(ncol(do.test)-2)], 
                    y = do.test[2, 2:(ncol(do.test)-2)])['p.value']
    mtx.out[which(mtx.out[, 1] == gene), 'p.val'] <- p.val
    if (display.progress) {
      setTxtProgressBar(pb = pb, value = which(genes.use == gene))
    }
    }, error=function(e){cat('ERROR: Gene', gene, 'In group', ident.test, conditionMessage(e), '\n')})
  }
  if (display.progress){
    close(con = pb)
  }
  mtx.out <- mtx.out[1:length(genes.use), ]
  #mtx.out <- mtx.out[which(mtx.out[, 'p.val'] <= 0.05), ]
  mtx.out$p.adj <- rep(1, nrow(mtx.out))
  p.to.adj <- which(mtx.out[, 'p.val'] <= 0.05)
  mtx.out[p.to.adj, 'p.adj'] <- p.adjust(p = mtx.out[p.to.adj, 'p.val'],
                                         n = nrow(object@data),
                                         method = "bonferroni")
  
  message('Sweeping Gene Sums')
  sum.mtx.test <- getPos(mat = object@data, genes.use = genes.use, c = cells[which(object@meta.data[, ident.name] %in% ident.test)])
  sum.mtx.ref <- getPos(mat = object@data, genes.use = genes.use, c = cells[which(object@meta.data[, ident.name] %in% ident.compare)])
  sum.mtx <- cbind(sum.mtx.test, sum.mtx.ref)
  sum.mtx <- cbind(sum.mtx, rownames(sum.mtx))
  colnames(sum.mtx) <- c(paste0('pct.', ident.test), 'pct.Others', 'Gene')
  sum.mtx <- sum.mtx[which(sum.mtx[, 'Gene'] %in% mtx.out[, 'Gene']), ]
  mtx.out <- merge(mtx.out, sum.mtx, by = 'Gene', sort = F)
  mtx.out <- mtx.out[order(mtx.out$p.adj, decreasing = F), ]
  rownames(mtx.out) <- seq(1:nrow(mtx.out))
  return(mtx.out[, c('Gene', paste0('pct.', ident.test), 'pct.Others', 'p.val', 'p.adj')])
}

aovTable <- function(object = NULL,
                     level.dep = 'ID',
                     level.ind = 'Treatment', 
                     rep.name = 'Replicate',
                     do.plot = F,
                     alternate.plot = T,
                     run.test = T,
                     prefix = 'C',
                     alternate.color = colors.use,
                     condition.first = "Vehicle", 
                     return.table = F){
  dat <- object@meta.data
  tests <- unique(dat[, paste(level.dep)])
  reps <- unique(dat[, paste(rep.name)])
  counts.mtx <- matrix(rep('0'), 
                       ncol = length(tests),
                       nrow = length(reps),
                       dimnames = list(reps, tests))
  getCounts <- function(b){
    r <- which(dat[, rep.name] == b)
    d <- sapply(tests,
                function(c) round(sum(dat[r, level.dep] == c)/length(r), 4))
    return(d)
  }
  counts.mtx <- t(sapply(reps, 
                         function(b) counts.mtx[b, ] = getCounts(b)))
  if(!is.null(prefix)){
    counts.mtx <- counts.mtx[,order(as.numeric(colnames(counts.mtx)))]
    colnames(counts.mtx) <- paste0(prefix, colnames(counts.mtx))
  }
  
  treatment <- sapply(reps, function(i) unique(dat[which(dat[, rep.name] == i), level.ind]))
  counts.mtx <- as.data.frame(cbind(counts.mtx, 'Treatment' = as.character(treatment)))
  counts.mtx[, 1:length(tests)] <- apply(counts.mtx[, 1:length(tests)], 2, function(x) as.numeric(as.character(x)))
  if(run.test){
    t <- sapply(X = colnames(counts.mtx)[1:length(tests)], 
                FUN = function(i) summary(aov(as.formula(paste0(i, ' ~ Treatment')), data = counts.mtx)))
   return(t) 
  }
    if(alternate.plot){
    counts.summary <- as.data.frame(counts.mtx[, 1:length(tests)])
    counts <- as.numeric(apply(counts.summary, 2, function(col) as.numeric(col)))
    counts <- cbind(counts, 'Level' = as.character(sapply(colnames(counts.summary), function(i) rep(i, length(reps)))))
    counts <- cbind(counts, 'Replicate' = rep(reps, length(tests)))
    counts <- cbind(counts, 'Condition' = rep(as.character(counts.mtx[, 'Treatment']), length(tests)))
    colnames(counts) <- c('Percent', colnames(counts)[2:4])
    counts <- as.data.frame(counts)
    counts$Percent <- as.numeric(as.character(counts$Percent))*100
    counts <- cbind(counts, 'sd' = rep(0, nrow(counts)))
    counts <- cbind(counts, 'mean' = rep(0, nrow(counts)))
    for(level in unique(counts$Level)){
      use <- which(counts$Level == level & counts$Condition == condition.first)
      counts[use, 'sd'] <- sd(counts[use, 'Percent'])
      counts[use, 'mean'] <- mean(counts[use, 'Percent'])
    }
    for(level in unique(counts$Level)){
      use <- which(counts$Level == level & counts$Condition != condition.first)
      counts[use, 'sd'] <- sd(counts[use, 'Percent'])
      counts[use, 'mean'] <- mean(counts[use, 'Percent'])
    }

    rownames(counts) <- paste0("r", 1:nrow(counts))
    counts.strip <- data.frame(matrix(ncol = ncol(counts), 
                                      nrow = length(unique(counts$Level))),
                               row.names = unique(counts$Level))
    colnames(counts.strip) <- colnames(counts)
    
    for(level in unique(counts$Level)){
      dat2 <- counts[which(counts$Condition == condition.first), ]
      ind <- rownames(dat2[which(dat2$Level == level), ])
      counts.strip[level, "Level"] <- level
      counts.strip[level, "Condition"] <- condition.first
      counts.strip[level, "mean"] <- as.numeric(dat2[ind, "mean"])[1]
      counts.strip[level, "sd"] <- as.numeric(dat2[ind, "sd"])[1]
    }
    counts.strip$Level <- factor(counts.strip$Level, levels = unique(counts.strip$Level))
    counts.strip$Condition <- factor(counts.strip$Condition, levels = unique(counts.strip$Condition))
    
    if(length(unique(counts$Condition)) > 1){
      condition.second <- as.character(counts[which(counts$Condition != condition.first), "Condition"])[1]
      counts.strip2 <- data.frame(matrix(ncol = ncol(counts), 
                                                         nrow = length(unique(counts$Level))),
                                                  row.names = unique(counts$Level))
      colnames(counts.strip) <- colnames(counts)
      for(level in unique(counts$Level)){
        dat2 <- counts[which(counts$Condition == condition.second), ]
        ind <- rownames(dat2[which(dat2$Level == level), ])
        counts.strip2[level, "Level"] <- level
        counts.strip2[level, "Condition"] <- condition.second 
        counts.strip2[level, "mean"] <- as.numeric(dat2[ind, "mean"])[1]
        counts.strip2[level, "sd"] <- as.numeric(dat2[ind, "sd"])[1]
      }
      counts$X <- paste0(counts$Level, "_", counts$Condition)
      counts.strip$X <- paste0(counts.strip$Level, "_", counts.strip$Condition)
      counts.strip2$X <- paste0(counts.strip2$Level, "_", counts.strip2$Condition)
      counts$C2 <- paste0(counts$Condition, "2")
      
      p <- ggplot(counts.strip, aes(x=X, ymin=mean-sd, ymax=mean+sd, color=Condition,
                                    middle = mean, min = mean, max = mean)) +
        geom_jitter(data = counts[which(counts$Condition == condition.first), ],
                    aes(y = Percent, color = C2), size = 2, show.legend = F) +
        geom_errorbar(width=.2, position=position_dodge(.8)) + 
        geom_boxplot(aes(y = mean), width = 0.2) +
        geom_jitter(data = counts[which(counts$Condition == condition.second), ],
                    aes(y = Percent, color = C2), size = 2, show.legend = F) +
        geom_errorbar(data = counts.strip2, width=.2, position=position_dodge(.8)) + 
        geom_boxplot(data = counts.strip2, aes(y = mean), width = 0.2) +
        scale_color_manual(values = rev(brewer.pal(4, 'Paired')))
        
        
      
      return(p)
    }else{
      counts.strip$Condition <- factor(counts.strip$Condition, levels = condition.first)
      p <- ggplot(counts.strip, aes(Level, mean, fill = Condition)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = alternate.color) +
        geom_jitter(data = counts, aes(Level, Percent, color = factor(Condition)), size = 2, show.legend = F) +
        scale_color_manual(values = c('black', grey(0, 0))) +
        geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                      position=position_dodge(.8))
      return(p)
    }

  }
  if(do.plot){
    counts.summary <- as.data.frame(counts.mtx[, 1:length(tests)])
    counts <- as.numeric(apply(counts.summary, 2, function(col) as.numeric(col)))
    counts <- cbind(counts, 'Level' = as.character(sapply(colnames(counts.summary), function(i) rep(i, length(reps)))))
    counts <- cbind(counts, 'Replicate' = rep(reps, length(tests)))
    counts <- cbind(counts, 'Condition' = rep(as.character(counts.mtx[, 'Treatment']), length(tests)))
    colnames(counts) <- c('Percent', colnames(counts)[2:4])
    counts <- as.data.frame(counts)
    counts$Percent <- as.numeric(as.character(counts$Percent))*100
    counts$Level <- factor(counts$Level, levels = c('C13', 'C7', 'C1', 'C2', 'C4', 'C3', "C0", 'Neurons', 'Astrocytes', 'Oligodendrocytes', 'Vascular_Fibroblasts', 'Microglia', 'Endothelial_Cells', "C6", "C5"))
    p <- ggplot(counts, aes(Level, Percent, fill = Condition)) +
      geom_boxplot(outlier.alpha = 0) +
      geom_jitter(aes(color = factor(Condition)), size = 2, show.legend = F)
    return(p)
  }else{
    sapply(X = colnames(counts.mtx)[1:length(tests)], 
           FUN = function(i) summary(aov(as.formula(paste0(i, ' ~ Treatment')), data = counts.mtx)))
    
  }
}










GetPopCounts <- function(object = NULL,
                         dep.group = 'Cluster',
                         ind.group = 'Replicate',
                         return.percent = F){
  getNumCells <- function(i, j){
    sum(object@meta.data[which(object@meta.data[, dep.group] == j), ind.group] == i)
  }
  pop.counts <- matrix(sample(0), ncol = length(unique(object@meta.data[, ind.group])),
                       nrow = length(unique(object@meta.data[, dep.group])),
                       dimnames = list(unique(object@meta.data[, dep.group]),
                                       unique(object@meta.data[, ind.group])))
  
  pop.counts <- sapply(colnames(pop.counts), 
                       FUN = function(i) sapply(rownames(pop.counts), 
                                                function(j) pop.counts[j, i] <- getNumCells(i, j)))
  if(! return.percent){
  return(pop.counts)
  }else{
    pop.sums <- colSums(pop.counts)
    pop.counts[1, ]/pop.sums
    pop.counts.percent <- t(apply(pop.counts, 1, function(x) x/pop.sums))*100
    return(round(pop.counts.percent, 2))
  }
}


getGeneSpaceCentroids <- function(object = NULL, genes.use = NULL, ident.use = NULL){
  require(Seurat)
  if(is.null(genes.use)){
    genes.use <- rownames(object@data)
  }
  object <- SetAllIdent(object, id = ident.use)
  gene.means <- sapply(levels(object@ident), FUN = function(i) 
                  apply(X = object@data[genes.use, WhichCells(object, i), drop = F], 
                        MARGIN = 1, FUN = function(x) log(x = mean(x = expm1(x = x)) + 1)))
  gene.means
}

plotGeneClusters <- function(object = Veh, 
                             genes.use = marker.genes, 
                             gene.clusters = genecl){
  require(Seurat)
  require(ape)
  if(is.null(genes.use)){
    genes.use <- rownames(object@data)
  }
  #Catch single-gene clusters and remove
  for(c in unique(gene.clusters)){
    if(length(which(gene.clusters == c)) <= 1){
      gene.clusters <- gene.clusters[which(gene.clusters != c)]
    }
  }
  
  gene.means <- sapply(levels(object@ident), FUN = function(i) 
    apply(X = object@data[genes.use, WhichCells(object, i), drop = F], 
          MARGIN = 1, FUN = function(x) log(x = mean(x = expm1(x = x)) + 1)))

  
  data.avg <- sapply(unique(gene.clusters), FUN = function(i) 
    apply(X = gene.means[which(gene.clusters == i), ], 
          MARGIN = 2, FUN = function(x) log(x = mean(x = expm1(x = x)) + 1)))
  colnames(data.avg) <- paste0('GeneCluster', unique(gene.clusters))
  data.dist <- dist(t(data.avg))
  gh <- hclust(d = data.dist)
  data.tree <- as.phylo(x = gh)
  plot.phylo(x = data.tree, direction = "downwards")
}



SimpleFeaturePlot <- function(Object = NULL, 
                              gene.use = NULL, 
                              color.use = 'red',
                              threshold = 0,
                              plot.intersection = F,
                              pt.size = 0.5,
                              alpha.colors = 0.9,
                              base.transparency = 0.7){
  require(Seurat)
  
      colorcode <- col2rgb(color.use)
      colorcode <- rgb(colorcode[1,], colorcode[2,], colorcode[3,], max=255, alpha = (255*alpha.colors))
    
    names(newcolors) <- c(grey(level = 0.7, alpha = base.transparency), colorcode)
    colors <- names(newcolors)
    colors <- colors[!is.na(colors)]
    
  #Get Plotting Coordinates
  frame <- GetDimReduction(Object, reduction.type = 'tsne', slot = "cell.embeddings")
  frame <- cbind(frame,
                 as.data.frame(matrix(data = 0, 
                                      nrow = nrow(frame),
                                      ncol = 3)))
  colnames(frame) <- c('tSNE_1', 'tSNE_2', gene.use, 'Plot.Status', 'is.null')
  
  positive_cells <- colnames(Object@data[, Object@data[gene.use, ] > threshold])
  frame[which(rownames(frame) %in% positive_cells), gene.use] <- '1'

  
  #Plot Intersections: Cells that have multiple of genes vector

      frame[which(frame[, gene.use] > 0), 'Plot.Status'] <- paste0(gene.use)
      frame[which(frame[, 'Plot.Status'] == 0), 'Plot.Status'] <- 'None'

    #frame[which(frame[, 'Plot.Status'] == 'None'), 'is.null'] <- '1'
    
    d <- ggplot(frame, aes(x = tSNE_1, y = tSNE_2)) +
      p +
      geom_point(data = frame[which(frame$Plot.Status != "None"), ], aes(x = tSNE_1, y = tSNE_2, color = Plot.Status), 
                 legend = T) +
      scale_colour_manual(values = c(colors.use[1:7],
                                     colors.use[8], 
                                     "green", 
                                     colors.use[9:10])) +
      ggtitle(title, subtitle = subtitle)
  return(d)
}

# Function for plotting colors side-by-side
pal <- function(col, border = "light gray", ...){
  n <- length(col)
  plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1),
       axes = FALSE, xlab = "", ylab = "", ...)
  rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)
}

PercentAbove <- function(x, threshold){
  return(length(x = x[x > threshold]) / length(x = x))
}

DotPlot2 <- function (object, genes.plot, cols.use = c("lightskyblue", "red"), 
                      col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
                      scale.by = "radius", scale.min = NA, scale.max = NA, group.by, 
                      plot.legend = FALSE, do.return = FALSE, x.lab.rot = FALSE, order = NA) 
{
  genes.plot <- rev(genes.plot)
  scale.func <- switch(EXPR = scale.by, size = scale_size, 
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  if (!missing(x = group.by)) {
    object <- SetAllIdent(object = object, id = group.by)
  }
  data.to.plot <- data.frame(FetchData(object = object, vars.all = genes.plot))
  colnames(x = data.to.plot) <- genes.plot
  data.to.plot$cell <- rownames(x = data.to.plot)
  data.to.plot$id <- object@ident
  if(!is.na(order[1])){
  data.to.plot$id <- factor(data.to.plot$id, levels = rev(order))
  }
  data.to.plot <- data.to.plot %>% gather(key = genes.plot, 
                                          value = expression, -c(cell, id))
  data.to.plot <- data.to.plot %>% group_by(id, genes.plot) %>% 
    summarize(avg.exp = mean(expm1(x = expression)), pct.exp = PercentAbove(x = expression, 
                                                                            threshold = 0))
  data.to.plot <- data.to.plot %>% ungroup() %>% group_by(genes.plot) %>% 
    mutate(avg.exp.scale = scale(x = avg.exp)) %>% mutate(avg.exp.scale = MinMax(data = avg.exp.scale, 
                                                                                 max = col.max, min = col.min))
  data.to.plot$genes.plot <- factor(x = data.to.plot$genes.plot, 
                                    levels = rev(x = genes.plot))
  data.to.plot$pct.exp[data.to.plot$pct.exp < dot.min] <- NA
  data.to.plot$pct.exp <- data.to.plot$pct.exp * 100
  p <- ggplot(data = data.to.plot, mapping = aes(x = genes.plot, 
                                                 y = id)) + geom_point(mapping = aes(size = pct.exp, 
                                                                                     color = avg.exp.scale)) + scale.func(range = c(0, dot.scale), 
                                                                                                                          limits = c(scale.min, scale.max)) + theme(axis.title.x = element_blank(), 
                                                                                                                                                                    axis.title.y = element_blank())
  if (length(x = cols.use) == 1) {
    p <- p + scale_color_distiller(palette = cols.use)
  }
  else {
    p <- p + scale_color_gradient(low = cols.use[1], high = cols.use[2])
  }
  if (!plot.legend) {
    p <- p + theme(legend.position = "none")
  }
  if (x.lab.rot) {
    p <- p + theme(axis.text.x = element_text(angle = 90, 
                                              vjust = 0.5))
  }
  suppressWarnings(print(p))
  if (do.return) {
    return(p)
  }
}

function (object, reduction.use = "pca", dim.1 = 1, dim.2 = 2, 
          cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE, 
          cols.use = NULL, group.by = "ident", pt.shape = NULL, do.hover = FALSE, 
          data.hover = "ident", do.identify = FALSE, do.label = FALSE, 
          label.size = 4, no.legend = FALSE, coord.fixed = FALSE, 
          no.axes = FALSE, dark.theme = FALSE, plot.order = NULL, 
          cells.highlight = NULL, cols.highlight = "red", sizes.highlight = 1, 
          plot.title = NULL, vector.friendly = FALSE, png.file = NULL, 
          png.arguments = c(10, 10, 100), na.value = "grey50", ...) 
{
  if (vector.friendly) {
    previous_call <- blank_call <- png_call <- match.call()
    blank_call$pt.size <- -1
    blank_call$do.return <- TRUE
    blank_call$vector.friendly <- FALSE
    png_call$no.axes <- TRUE
    png_call$no.legend <- TRUE
    png_call$do.return <- TRUE
    png_call$vector.friendly <- FALSE
    png_call$plot.title <- NULL
    blank_plot <- eval(blank_call, sys.frame(sys.parent()))
    png_plot <- eval(png_call, sys.frame(sys.parent()))
    png.file <- SetIfNull(x = png.file, default = paste0(tempfile(), 
                                                         ".png"))
    ggsave(filename = png.file, plot = png_plot, width = png.arguments[1], 
           height = png.arguments[2], dpi = png.arguments[3])
    to_return <- AugmentPlot(plot1 = blank_plot, imgFile = png.file)
    file.remove(png.file)
    if (do.return) {
      return(to_return)
    }
    else {
      print(to_return)
    }
  }
  embeddings.use <- GetDimReduction(object = object, reduction.type = reduction.use, 
                                    slot = "cell.embeddings")
  if (length(x = embeddings.use) == 0) {
    stop(paste(reduction.use, "has not been run for this object yet."))
  }
  cells.use <- SetIfNull(x = cells.use, default = colnames(x = object@data))
  dim.code <- GetDimReduction(object = object, reduction.type = reduction.use, 
                              slot = "key")
  dim.codes <- paste0(dim.code, c(dim.1, dim.2))
  data.plot <- as.data.frame(x = embeddings.use)
  cells.use <- intersect(x = cells.use, y = rownames(x = data.plot))
  data.plot <- data.plot[cells.use, dim.codes]
  ident.use <- as.factor(x = object@ident[cells.use])
  if (group.by != "ident") {
    ident.use <- as.factor(x = FetchData(object = object, 
                                         vars.all = group.by)[cells.use, 1])
  }
  data.plot$ident <- ident.use
  data.plot$x <- data.plot[, dim.codes[1]]
  data.plot$y <- data.plot[, dim.codes[2]]
  data.plot$pt.size <- pt.size
  if (!is.null(x = cells.highlight)) {
    if (is.character(x = cells.highlight)) {
      cells.highlight <- list(cells.highlight)
    }
    else if (is.data.frame(x = cells.highlight) || !is.list(x = cells.highlight)) {
      cells.highlight <- as.list(x = cells.highlight)
    }
    cells.highlight <- lapply(X = cells.highlight, FUN = function(cells) {
      cells.return <- if (is.character(x = cells)) {
        cells[cells %in% rownames(x = data.plot)]
      }
      else {
        cells <- as.numeric(x = cells)
        cells <- cells[cells <= nrow(x = data.plot)]
        rownames(x = data.plot)[cells]
      }
      return(cells.return)
    })
    cells.highlight <- Filter(f = length, x = cells.highlight)
    if (length(x = cells.highlight) > 0) {
      if (!no.legend) {
        no.legend <- is.null(x = names(x = cells.highlight))
      }
      names.highlight <- if (is.null(x = names(x = cells.highlight))) {
        paste0("Group_", 1L:length(x = cells.highlight))
      }
      else {
        names(x = cells.highlight)
      }
      sizes.highlight <- rep_len(x = sizes.highlight, 
                                 length.out = length(x = cells.highlight))
      cols.highlight <- rep_len(x = cols.highlight, length.out = length(x = cells.highlight))
      highlight <- rep_len(x = NA_character_, length.out = nrow(x = data.plot))
      if (is.null(x = cols.use)) {
        cols.use <- "black"
      }
      cols.use <- c(cols.use[1], cols.highlight)
      size <- rep_len(x = pt.size, length.out = nrow(x = data.plot))
      for (i in 1:length(x = cells.highlight)) {
        cells.check <- cells.highlight[[i]]
        index.check <- match(x = cells.check, rownames(x = data.plot))
        highlight[index.check] <- names.highlight[i]
        size[index.check] <- sizes.highlight[i]
      }
      plot.order <- sort(x = unique(x = highlight), na.last = TRUE)
      plot.order[is.na(x = plot.order)] <- "Unselected"
      highlight[is.na(x = highlight)] <- "Unselected"
      highlight <- as.factor(x = highlight)
      data.plot$ident <- highlight
      data.plot$pt.size <- size
      if (dark.theme) {
        cols.use[1] <- "white"
      }
    }
  }
  if (!is.null(x = plot.order)) {
    if (any(!plot.order %in% data.plot$ident)) {
      stop("invalid ident in plot.order")
    }
    plot.order <- rev(x = c(plot.order, setdiff(x = unique(x = data.plot$ident), 
                                                y = plot.order)))
    data.plot$ident <- factor(x = data.plot$ident, levels = plot.order)
    data.plot <- data.plot[order(data.plot$ident), ]
  }
  p <- ggplot(data = data.plot, mapping = aes(x = x, y = y)) + 
    geom_point(mapping = aes(colour = factor(x = ident), 
                             size = pt.size))
  if (!is.null(x = pt.shape)) {
    shape.val <- FetchData(object = object, vars.all = pt.shape)[cells.use, 
                                                                 1]
    if (is.numeric(shape.val)) {
      shape.val <- cut(x = shape.val, breaks = 5)
    }
    data.plot[, "pt.shape"] <- shape.val
    p <- ggplot(data = data.plot, mapping = aes(x = x, y = y)) + 
      geom_point(mapping = aes(colour = factor(x = ident), 
                               shape = factor(x = pt.shape), size = pt.size))
  }
  if (!is.null(x = cols.use)) {
    p <- p + scale_colour_manual(values = cols.use, na.value = na.value)
  }
  if (coord.fixed) {
    p <- p + coord_fixed()
  }
  p <- p + guides(size = FALSE)
  p2 <- p + xlab(label = dim.codes[[1]]) + ylab(label = dim.codes[[2]]) + 
    scale_size(range = c(min(data.plot$pt.size), max(data.plot$pt.size)))
  p3 <- p2 + SetXAxisGG() + SetYAxisGG() + SetLegendPointsGG(x = 6) + 
    SetLegendTextGG(x = 12) + no.legend.title + theme_bw() + 
    NoGrid()
  if (dark.theme) {
    p <- p + DarkTheme()
    p3 <- p3 + DarkTheme()
  }
  p3 <- p3 + theme(legend.title = element_blank())
  if (!is.null(plot.title)) {
    p3 <- p3 + ggtitle(plot.title) + theme(plot.title = element_text(hjust = 0.5))
  }
  if (do.label) {
    centers <- data.plot %>% dplyr::group_by(ident) %>% 
      summarize(x = median(x = x), y = median(x = y))
    p3 <- p3 + geom_point(data = centers, mapping = aes(x = x, 
                                                        y = y), size = 0, alpha = 0) + geom_text(data = centers, 
                                                                                                 mapping = aes(label = ident), size = label.size)
  }
  if (no.legend) {
    p3 <- p3 + theme(legend.position = "none")
  }
  if (no.axes) {
    p3 <- p3 + theme(axis.line = element_blank(), axis.text.x = element_blank(), 
                     axis.text.y = element_blank(), axis.ticks = element_blank(), 
                     axis.title.x = element_blank(), axis.title.y = element_blank(), 
                     panel.background = element_blank(), panel.border = element_blank(), 
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     plot.background = element_blank())
  }
  if (do.identify || do.hover) {
    if (do.bare) {
      plot.use <- p
    }
    else {
      plot.use <- p3
    }
    if (do.hover) {
      if (is.null(x = data.hover)) {
        features.info <- NULL
      }
      else {
        features.info <- FetchData(object = object, 
                                   vars.all = data.hover)
      }
      return(HoverLocator(plot = plot.use, data.plot = data.plot, 
                          features.info = features.info, dark.theme = dark.theme))
    }
    else if (do.identify) {
      return(FeatureLocator(plot = plot.use, data.plot = data.plot, 
                            dark.theme = dark.theme, ...))
    }
  }
  if (do.return) {
    if (do.bare) {
      return(p)
    }
    else {
      return(p3)
    }
  }
  if (do.bare) {
    print(p)
  }
  else {
    print(p3)
  }
}

MakeCountSummary <- function(object = NULL,
                             level.dep = 'Cluster',
                             level.ind = 'Treatment', 
                             rep.name = 'Replicate',
                             prefix = 'C'){
  dat <- object@meta.data
  tests <- unique(dat[, paste(level.dep)])
  reps <- unique(dat[, paste(rep.name)])
  counts.mtx <- matrix(rep('0'), 
                       ncol = length(tests),
                       nrow = length(reps),
                       dimnames = list(reps, tests))
  getCounts <- function(b){
    r <- which(dat[, rep.name] == b)
    d <- sapply(tests,
                function(c) round(sum(dat[r, level.dep] == c)/length(r), 4))
    return(d)
  }
  counts.mtx <- t(sapply(reps, 
                         function(b) counts.mtx[b, ] = getCounts(b)))
  
  counts.mtx <- counts.mtx[,order(as.numeric(colnames(counts.mtx)))]
  colnames(counts.mtx) <- paste0(prefix, colnames(counts.mtx))
  treatment <- sapply(reps, function(i) unique(dat[which(dat[, rep.name] == i), level.ind]))
  counts.mtx <- as.data.frame(cbind(counts.mtx, 'Treatment' = treatment))
  counts.mtx[, 1:length(tests)] <- apply(counts.mtx[, 1:length(tests)], 2, function(x) as.numeric(as.character(x)))
  sapply(X = colnames(counts.mtx)[1:length(tests)], 
         FUN = function(i) summary(aov(as.formula(paste0(i, ' ~ Treatment')), data = counts.mtx)))
  counts.summary <- as.data.frame(counts.mtx[, 1:length(tests)])
  counts <- as.numeric(apply(counts.summary, 2, function(col) as.numeric(col)))
  counts <- cbind(counts, 'Level' = as.character(sapply(colnames(counts.summary), function(i) rep(i, length(reps)))))
  counts <- cbind(counts, 'Replicate' = rep(reps, length(tests)))
  counts <- cbind(counts, 'Condition' = rep(as.character(counts.mtx[, 'Treatment']), length(tests)))
  colnames(counts) <- c('Percent', colnames(counts)[2:4])
  counts <- as.data.frame(counts)
  counts$Percent <- as.numeric(as.character(counts$Percent))*100
  counts$Level <- factor(counts$Level, levels = unique(counts$Level))
  return(counts)
}

MakeCountSummarySimple <- function(object = GT,
                             level.dep = 'Cluster',
                             rep.name = 'Replicate',
                             prefix = 'C'){
  dat <- object@meta.data
  tests <- unique(dat[, paste(level.dep)])
  reps <- unique(dat[, paste(rep.name)])
  counts.mtx <- matrix(rep('0'), 
                       ncol = length(tests),
                       nrow = length(reps),
                       dimnames = list(reps, tests))
  getCounts <- function(b){
    r <- which(dat[, rep.name] == b)
    d <- sapply(tests,
                function(c) round(sum(dat[r, level.dep] == c)/length(r), 4))
    return(d)
  }
  counts.mtx <- t(sapply(reps, 
                         function(b) counts.mtx[b, ] = getCounts(b)))
  
  counts.mtx <- counts.mtx[,order(as.numeric(colnames(counts.mtx)))]
  colnames(counts.mtx) <- paste0(prefix, colnames(counts.mtx))
  counts.mtx[, 1:length(tests)] <- apply(counts.mtx[, 1:length(tests)], 2, function(x) as.numeric(as.character(x)))
  counts.summary <- as.data.frame(counts.mtx[, 1:length(tests)])
  counts <- as.numeric(apply(counts.summary, 2, function(col) as.numeric(col)))
  counts <- cbind(counts, 'Level' = as.character(sapply(colnames(counts.summary), function(i) rep(i, length(reps)))))
  counts <- cbind(counts, 'Replicate' = rep(reps, length(tests)))
  colnames(counts) <- c('Percent', colnames(counts)[-1])
  counts <- as.data.frame(counts)
  counts$Percent <- as.numeric(as.character(counts$Percent))*100
  counts$Level <- factor(counts$Level, levels = unique(counts$Level))
  return(counts)
}

RunFrame <- function(Object = NULL, 
                     gene.use = NULL, 
                     threshold = 0){
  require(Seurat)
  
  #Get Plotting Coordinates
  frame <- GetDimReduction(Object, reduction.type = 'tsne', slot = "cell.embeddings")
  frame <- cbind(frame,
                 as.data.frame(matrix(data = 0, 
                                      nrow = nrow(frame),
                                      ncol = 3)))
  colnames(frame) <- c('tSNE_1', 'tSNE_2', gene.use, 'Plot.Status', 'is.null')
  
  positive_cells <- colnames(Object@data[, Object@data[gene.use, ] > threshold])
  frame[which(rownames(frame) %in% positive_cells), gene.use] <- '1'
  
  
  frame[which(frame[, gene.use] > 0), 'Plot.Status'] <- paste0(gene.use)
  frame[which(frame[, 'Plot.Status'] == 0), 'Plot.Status'] <- 'None'
  return(frame)
}

MakeQCPlots <- function(Object, n = 10, pc.max = 100){
  message(paste0("Running IRLBA/SVD PCA for ", pc.max, " components..."))
  Object <- RunPCA(Object, pcs.compute = pc.max, do.print = F)
  message(paste0("Running tSNE for ", n, " components..."))
  Object <- RunTSNE(Object, dims.use = 1:n)
  TSNEPlot(Object, colors.use = colors.use)
  message("Building QC plots...")
  
  plots <- list()
  
  plots[[1]] <- PCElbowPlot(Object, num.pc = pc.max)
  Object <- SetAllIdent(Object, "Replicate")
  
  
  plots[[2]] <- TSNEPlot(Object, colors.use = colors.use, pt.size = 0.3, 
                         plot.title = paste0("k=", n, 
                                             "; ", n, " PC; ",
                                             "All Cells; Replicate Effects"),
                         do.return = T)
  
  
  Object <- SetAllIdent(Object, "Sex")
  plots[[3]] <- TSNEPlot(Object, 
                         plot.title = paste0("k=", n, 
                                             "; ", n, " PC; ",
                                             "All Cells; Sex Effects"),
                         colors.use = viridis(3)[1:2], pt.size = 0.3, do.return = T)
  
  
  Object <- SetAllIdent(Object, "Date")
  plots[[4]] <- TSNEPlot(Object, 
                         plot.title = paste0("k=", n, 
                                             "; ", n, " PC; ",
                                             "All Cells; Date of Run Effects"),
                         colors.use = viridis(3)[1:3], pt.size = 0.3, do.return = T)
  
  
  Object <- SetAllIdent(Object, "ID")
  colors2 <- c(colors.use[1:4], colors.use[c(15, 12:14, 16:17)])
  plots[[5]] <- TSNEPlot(Object, colors.use = colors2, 
                         plot.title = paste0("k=", n, 
                                             "; ", n, " PC; ",
                                             "All Cells; Cell Class"),
                         plot.order = rev(c("Node A", "Node B", "Node C", "Node D",
                                            "Neurons", "Astrocytes", "Oligodendrocytes",
                                            "Microglia", "Vascular Fibroblasts",
                                            "Endothelial Cells")), 
                         pt.size = 0.3, do.return = T)
  
  
  genes <- c("Ccnd1", "Ccne2", "Ccna2", "Ccnb2")
  thresh <- c(2,1,2,2)
  names(thresh) <- genes
  colors2 <- brewer.pal(4,'Dark2')
  plots[[6]] <- MultiFeaturePlot(Object, gene.list = genes, colors = colors2, 
                                 threshold = thresh, pt.size = 0.3, 
                                 title = paste0("k=", n, 
                                                     "; ", n, " PC; ",
                                                     "All Cells; Cell Cycle"))
  file <- paste0("k", n, "_", n, "PC_AllCells.svg")
  devSVG(file = file, height = 11, width = 8.5)
  plot_grid(plots[[1]], plots[[2]],
            plots[[3]], plots[[4]],
            plots[[5]], plots[[6]],ncol = 2)
  dev.off() 
  message(pate0("QC plots saved to: ", file))
}

aovTable.norm <- function(object = NULL,
                     level.dep = 'ID',
                     level.ind = 'Treatment', 
                     rep.name = 'Replicate',
                     do.plot = F,
                     alternate.plot = T,
                     run.test = T,
                     prefix = 'C',
                     alternate.color = colors.use,
                     condition.first = "Vehicle", 
                     return.table = F){
  dat <- object@meta.data
  tests <- unique(dat[, paste(level.dep)])
  reps <- unique(dat[, paste(rep.name)])
  counts.mtx <- matrix(rep('0'), 
                       ncol = length(tests),
                       nrow = length(reps),
                       dimnames = list(reps, tests))
  getCounts <- function(b){
    r <- which(dat[, rep.name] == b)
    d <- sapply(tests,
                function(c) round(sum(dat[r, level.dep] == c)/length(r), 4))
    return(d)
  }
  counts.mtx <- t(sapply(reps, 
                         function(b) counts.mtx[b, ] = getCounts(b)))
  
  #counts.mtx <- counts.mtx[,order(as.numeric(colnames(counts.mtx)))]
  #colnames(counts.mtx) <- paste0(prefix, colnames(counts.mtx))
  treatment <- sapply(reps, function(i) unique(dat[which(dat[, rep.name] == i), level.ind]))
  counts.mtx <- as.data.frame(cbind(counts.mtx, 'Treatment' = treatment))
  counts.mtx[, 1:length(tests)] <- apply(counts.mtx[, 1:length(tests)], 2, function(x) as.numeric(as.character(x)))
  if(run.test){
    t <- sapply(X = colnames(counts.mtx)[1:length(tests)], 
                FUN = function(i) summary(aov(as.formula(paste0(i, ' ~ Treatment')), data = counts.mtx)))
    return(t) 
  }
  if(alternate.plot){
    counts.summary <- as.data.frame(counts.mtx[, 1:length(tests)])
    counts <- as.numeric(apply(counts.summary, 2, function(col) as.numeric(col)))
    counts <- cbind(counts, 'Level' = as.character(sapply(colnames(counts.summary), function(i) rep(i, length(reps)))))
    counts <- cbind(counts, 'Replicate' = rep(reps, length(tests)))
    counts <- cbind(counts, 'Condition' = rep(as.character(counts.mtx[, 'Treatment']), length(tests)))
    colnames(counts) <- c('Percent', colnames(counts)[2:4])
    counts <- as.data.frame(counts)
    counts$Percent <- as.numeric(as.character(counts$Percent))*100
    counts <- cbind(counts, 'sd' = rep(0, nrow(counts)))
    counts <- cbind(counts, 'mean' = rep(0, nrow(counts)))
    for(level in unique(counts$Level)){
      use <- which(counts$Level == level & counts$Condition == condition.first)
      counts[use, 'sd'] <- sd(counts[use, 'Percent'])
      counts[use, 'mean'] <- mean(counts[use, 'Percent'])
    }
    for(level in unique(counts$Level)){
      use <- which(counts$Level == level & counts$Condition != condition.first)
      counts[use, 'sd'] <- sd(counts[use, 'Percent'])
      counts[use, 'mean'] <- mean(counts[use, 'Percent'])
    }
    
    rownames(counts) <- paste0("r", 1:nrow(counts))
    counts.strip <- data.frame(matrix(ncol = ncol(counts), 
                                      nrow = length(unique(counts$Level))),
                               row.names = unique(counts$Level))
    colnames(counts.strip) <- colnames(counts)
    
    for(level in unique(counts$Level)){
      dat2 <- counts[which(counts$Condition == condition.first), ]
      ind <- rownames(dat2[which(dat2$Level == level), ])
      counts.strip[level, "Level"] <- level
      counts.strip[level, "Condition"] <- condition.first
      counts.strip[level, "mean"] <- as.numeric(dat2[ind, "mean"])[1]
      counts.strip[level, "sd"] <- as.numeric(dat2[ind, "sd"])[1]
    }
    counts.strip$Level <- factor(counts.strip$Level, levels = unique(counts.strip$Level))
    counts.strip$Condition <- factor(counts.strip$Condition, levels = unique(counts.strip$Condition))
    
    if(length(unique(counts$Condition)) > 1){
      condition.second <- as.character(counts[which(counts$Condition != condition.first), "Condition"])[1]
      counts.strip2 <- data.frame(matrix(ncol = ncol(counts), 
                                         nrow = length(unique(counts$Level))),
                                  row.names = unique(counts$Level))
      colnames(counts.strip) <- colnames(counts)
      for(level in unique(counts$Level)){
        dat2 <- counts[which(counts$Condition == condition.second), ]
        ind <- rownames(dat2[which(dat2$Level == level), ])
        counts.strip2[level, "Level"] <- level
        counts.strip2[level, "Condition"] <- condition.second 
        counts.strip2[level, "mean"] <- as.numeric(dat2[ind, "mean"])[1]
        counts.strip2[level, "sd"] <- as.numeric(dat2[ind, "sd"])[1]
      }
      counts$X <- paste0(counts$Level, "_", counts$Condition)
      counts.strip$X <- paste0(counts.strip$Level, "_", counts.strip$Condition)
      counts.strip2$X <- paste0(counts.strip2$Level, "_", counts.strip2$Condition)
      counts$C2 <- paste0(counts$Condition, "2")
      
      p <- ggplot(counts.strip, aes(x=X, ymin=mean-sd, ymax=mean+sd, color=Condition,
                                    middle = mean, min = mean, max = mean)) +
        geom_jitter(data = counts[which(counts$Condition == condition.first), ],
                    aes(y = Percent, color = C2), size = 2, show.legend = F) +
        geom_errorbar(width=.2, position=position_dodge(.8)) + 
        geom_boxplot(aes(y = mean), width = 0.2) +
        geom_jitter(data = counts[which(counts$Condition == condition.second), ],
                    aes(y = Percent, color = C2), size = 2, show.legend = F) +
        geom_errorbar(data = counts.strip2, width=.2, position=position_dodge(.8)) + 
        geom_boxplot(data = counts.strip2, aes(y = mean), width = 0.2) +
        scale_color_manual(values = rev(brewer.pal(4, 'Paired')))
      
      
      
      return(p)
    }else{
      counts.strip$Condition <- factor(counts.strip$Condition, levels = condition.first)
      p <- ggplot(counts.strip, aes(Level, mean, fill = Condition)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = alternate.color) +
        geom_jitter(data = counts, aes(Level, Percent, color = factor(Condition)), size = 2, show.legend = F) +
        scale_color_manual(values = c('black', grey(0, 0))) +
        geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                      position=position_dodge(.8))
      return(p)
    }
    
  }
  if(do.plot){
    counts.summary <- as.data.frame(counts.mtx[, 1:length(tests)])
    counts <- as.numeric(apply(counts.summary, 2, function(col) as.numeric(col)))
    counts <- cbind(counts, 'Level' = as.character(sapply(colnames(counts.summary), function(i) rep(i, length(reps)))))
    counts <- cbind(counts, 'Replicate' = rep(reps, length(tests)))
    counts <- cbind(counts, 'Condition' = rep(as.character(counts.mtx[, 'Treatment']), length(tests)))
    colnames(counts) <- c('Percent', colnames(counts)[2:4])
    counts <- as.data.frame(counts)
    counts$Percent <- as.numeric(as.character(counts$Percent))*100
    counts$Level <- factor(counts$Level, levels = c('C13', 'C7', 'C1', 'C2', 'C4', 'C3', "C0", 'Neurons', 'Astrocytes', 'Oligodendrocytes', 'Vascular_Fibroblasts', 'Microglia', 'Endothelial_Cells', "C6", "C5"))
    p <- ggplot(counts, aes(Level, Percent, fill = Condition)) +
      geom_boxplot(outlier.alpha = 0) +
      geom_jitter(aes(color = factor(Condition)), size = 2, show.legend = F)
    return(p)
  }else{
    sapply(X = colnames(counts.mtx)[1:length(tests)], 
           FUN = function(i) summary(aov(as.formula(paste0(i, ' ~ Treatment')), data = counts.mtx)))
    
  }
}








