### Custom functions by KJ ###

library(ggrepel)
library(cowplot)
library(patchwork)
library(ggplot2)
library(ggmosaic)
library(scales)
library(gplots)

options(warn=1)


###########################################################################

mycol <- c("burlywood4", "antiquewhite4", "antiquewhite3","burlywood2", "bisque1", "lightsalmon1",
           "khaki1", "gold", "orange", "#E38900", "chocolate", "orange4",  
           "mediumorchid3", "#A58AFF", "plum2", "orchid1", "hotpink1", "lightpink3", "lightpink1", "coral1", "indianred3",
           "#06A4FF","cornflowerblue", "slategray1", "cadetblue1", 
           "#99A800", "chartreuse3", "#00BC56", "#00C094","mediumturquoise","aquamarine", "darkseagreen1", "darkolivegreen1", "darkkhaki")
byr <- rev(RColorBrewer::brewer.pal(11,"RdYlBu")) %>% colorRampPalette                            
ryb <- RColorBrewer::brewer.pal(11,"RdYlBu") %>% colorRampPalette                    
wr <- colorRampPalette(rev(c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#FFFFBF", "white")))
wb  <- colorRampPalette(c("white", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4", "#313695"))
bgyor2 <- colorRampPalette(c("#A20543","#B01646","#D33C4E","#DD4B4B","#E65847","#F06943","#F67E4B","#FAA059","#FA9B58",
                            "#FDB567","#FDC272","#FDD885","#F8E38D","#EFEB91","#F1EA91","#E7F397","#E7F397","#D4ED9B",
                            "#C7E89E","#B2DFA1","#A0D7A4","#8CD0A4","#74CCA4","#65C0A4","#55AFAC","#469EB2","#388EB9",
                            "#387EB8","#456EB0","#4E61AA","#5A52A3") %>% rev())
spectral <- colorRampPalette(brewer.pal(11,"Spectral") %>% rev())

############################################################################


'%!in%' <- function(x,y)!('%in%'(x,y))


############################################################################


convertHumanGeneList <- function(x) {
  require("biomaRt")
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 <- getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x, mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows = T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}


############################################################################


PFlog1pPF = function(x){
  y = max(1, sum(x))
  y = log(1 + x/y)*1e9/y
  return(y)
}



############################################################################


RenameMetaCol <- function(seuObj,metaCol,namesOld,namesNew,levels=NULL) {
  
  tmp <- seuObj@meta.data[,metaCol] %>% as.character
  
  for (i in 1:length(namesOld)){
    tmp[tmp == namesOld[[i]]] <- namesNew[[i]]
  }
  
  if (is.null(levels)){
    seuObj@meta.data[,metaCol] <- factor(tmp) 
  } else {
    seuObj@meta.data[,metaCol] <- factor(tmp,levels=levels) 
  }
  Idents(seuObj) <- metaCol
  return(seuObj)
}


############################################################################


scale01 <- function(x){(x-min(x))/(max(x)-min(x))}


############################################################################


notUniqueGenes <- function(object) {
  rownames(object)[!(isUnique(rownames(object)))]
}


############################################################################


sumNotUniqueGenes <- function(macExp) {
  genList <- unique(notUniqueGenes(macExp))
  myList <- lapply(genList, function(gene) {
    list(name = gene, sums = colSums(macExp[rownames(macExp) == gene,]))
  })
  tmp <- do.call(rbind, lapply(myList, function(x) x$sums))
  rownames(tmp) <- genList
  macExp <- rbind(macExp[!(rownames(macExp) %in% genList), ], tmp)
  return(macExp)
}


############################################################################


wcsv <- function(data.frame,filename="my.csv",append= "FALSE",col_names="TRUE", num_threads = readr_threads()){
  data.frame <- cbind(rownames(data.frame),data.frame)
  write_csv(data.frame,filename,append)
}


############################################################################


addNewlines <- function(text,nChar=25,removeText="KEGG_") {
  # Remove "KEGG" prefix
  cleaned_text <- gsub(paste("^",removeText,sep = ""), "", text)
  
  # Replace underscores with spaces
  cleaned_text <- gsub("_", " ", cleaned_text)
  
  # Split the text into words
  words <- strsplit(cleaned_text, " ")[[1]]
  
  # Initialize variables
  current_length <- 0
  modified_text <- ""
  
  for (word in words) {
    word_length <- nchar(word)
    
    # If adding the current word exceeds the desired length, add a newline
    if (current_length + word_length > nChar) {
      modified_text <- paste0(modified_text, "\n")
      current_length <- 0
    }
    
    # Add the word and update the current length
    modified_text <- paste0(modified_text, word, " ")
    current_length <- current_length + word_length + 1
  }
  
  return(trimws(modified_text))
}


############################################################################


find.pK <- function(sweep.stats) # suppres creation of graphics
{
  "%ni%" <- Negate("%in%")
  if ("AUC" %ni% colnames(sweep.stats) == TRUE) {
    bc.mvn <- as.data.frame(matrix(0L,
                                   nrow = length(unique(sweep.stats$pK)),
                                   ncol = 5
    ))
    colnames(bc.mvn) <- c(
      "ParamID", "pK", "MeanBC", "VarBC",
      "BCmetric"
    )
    bc.mvn$pK <- unique(sweep.stats$pK)
    bc.mvn$ParamID <- 1:nrow(bc.mvn)
    x <- 0
    for (i in unique(bc.mvn$pK)) {
      x <- x + 1
      ind <- which(sweep.stats$pK == i)
      bc.mvn$MeanBC[x] <- mean(sweep.stats[ind, "BCreal"])
      bc.mvn$VarBC[x] <- sd(sweep.stats[ind, "BCreal"])^2
      bc.mvn$BCmetric[x] <- mean(sweep.stats[ind, "BCreal"]) / (sd(sweep.stats[
        ind,
        "BCreal"
      ])^2)
    }
    par(mar = rep(1, 4))
    return(bc.mvn)
  }
  if ("AUC" %in% colnames(sweep.stats) == TRUE) {
    bc.mvn <- as.data.frame(matrix(0L,
                                   nrow = length(unique(sweep.stats$pK)),
                                   ncol = 6
    ))
    colnames(bc.mvn) <- c(
      "ParamID", "pK", "MeanAUC", "MeanBC",
      "VarBC", "BCmetric"
    )
    bc.mvn$pK <- unique(sweep.stats$pK)
    bc.mvn$ParamID <- 1:nrow(bc.mvn)
    x <- 0
    for (i in unique(bc.mvn$pK)) {
      x <- x + 1
      ind <- which(sweep.stats$pK == i)
      bc.mvn$MeanAUC[x] <- mean(sweep.stats[ind, "AUC"])
      bc.mvn$MeanBC[x] <- mean(sweep.stats[ind, "BCreal"])
      bc.mvn$VarBC[x] <- sd(sweep.stats[ind, "BCreal"])^2
      bc.mvn$BCmetric[x] <- mean(sweep.stats[ind, "BCreal"]) / (sd(sweep.stats[
        ind,
        "BCreal"
      ])^2)
    }
    par(mar = rep(1, 4))
    x <- plot(
      x = bc.mvn$ParamID, y = bc.mvn$MeanAUC, pch = 18,
      col = "black", cex = 0.75, xlab = NA, ylab = NA
    )
    x <- lines(
      x = bc.mvn$ParamID, y = bc.mvn$MeanAUC, col = "black",
      lty = 2
    )
    par(new = TRUE)
    return(bc.mvn)
  }
}



############################################################################


library(viridis)
umapp <- function(macExp, i) {
  theme_feature <- theme_classic() +
    theme(
      axis.text = element_text(24), axis.ticks = element_blank(),
      axis.title = element_text(size = 20),
      legend.position = "bottom", legend.text = element_text(size = 18), legend.title = element_text(size = 20),
      plot.title = element_text(face = "italic", size = 60),
      legend.key.size = unit(1, "cm")
    )
  
  plot <- ggplot(macExp, aes(x = wnnUMAP_1, y = wnnUMAP_2)) +
    geom_jitter(size = 0.5, alpha = 0.5, aes_string(color = as.name(i))) + # czemu okregi a nie koła?
    labs(title = i, color = "log2(RPKM+1)") +
    scale_colour_gradientn(colours = viridis(50)) + # rev(inferno(50))
    # breaks=c(min,max))+
    # , labels=c("min", "max"))+
    xlab("wnnUMAP_1") +
    ylab("wnnUMAP_2") +
    theme_feature
  return(plot)
}



############################################################################


ensembleToHugo <- function(seuObj,assay){
  require("biomaRt")
  seuObj@active.assay <- assay
  mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl")
  lol <- getBM(filters = "ensembl_gene_id", 
               attributes = c("ensembl_gene_id", "external_gene_name"),
               values = rownames(seuObj@assays[["RNA"]]@data), mart = mart)
  ochseg <- subset(seuObj, features = seuObj$ensembl_gene_id )
  mygenes <- lol$ensembl_gene_id[match(lol$ensembl_gene_id, rownames(lol))]
  rownames(seuObj@assays[["RNA"]]@data) <- mygenes
}



############################################################################


testmac <- function(macExp, j) {
  for (i in j) {
    if (identical(colnames(macExp)[which(grepl(j, colnames(macExp)))], character(0))) {
      print(paste(i, " ", "is not in genes pool"))
    } else {
      return(colnames(macExp)[which(grepl(j, colnames(macExp)))])
    }
  }
}



############################################################################

testmacrev <- function(macExp, j) {
  for (i in j) {
    if (identical(rownames(macExp)[which(grepl(j, rownames(macExp)))], character(0))) {
      print(paste(i, " ", "is not in genes pool"))
    } else {
      return(rownames(macExp)[which(grepl(j, rownames(macExp)))])
    }
  }
}



############################################################################

testseu <- function(seuobj, j) {
  geneList <- list()
  n <- 0
  for (i in j) {
    n <- n + 1
    if (identical(rownames(seuobj)[which(grepl(i, rownames(seuobj)))], character(0))) {
      print(paste(i, " ", "is not in genes pool"))
    } else {
      print(rownames(seuobj)[which(grepl(i, rownames(seuobj)))])
      geneList[[n]] <- rownames(seuobj)[which(grepl(i, rownames(seuobj)))]
    }
  }
  return(unlist(geneList))
}



############################################################################


h2 <- function(i) {
  i[1:5, 1:5]
}


############################################################################


perc99 <- function(macExp) {
  for (i in colnames(macExp[, -c(1:4)])) {
    q99 <- quantile(macExp[which(macExp[, i] > 0), i], probs = 0.99)
    macExp[(which(macExp[, i] > q99)), i] <- q99
  }
  return(macExp)
}


############################################################################


greyClusters <- function(seuObj, idents, redu = "wnn.umap") {
  cellsInCluster <- WhichCells(seuObj, idents = idents)
  DimPlot(seuObj,
          reduction = redu,
          label = F, group.by = "ident",
          cells.highlight = list(cellsInCluster), cols.highlight = c("darkred"),
          cols = "grey", pt.size = 0.3
  ) + ggtitle(idents) + theme(legend.position = "none")
}



############################################################################


percTab <- function(alldata) {
  table(Idents(alldata), alldata$lab) %>%
    prop.table(margin = 2) %>%
    round(digits = 4)
}


############################################################################


mybarplot <- function(height, x = "Count", color = "p.adjust", showCategory = 8,
                      font.size = 12, title = "", label_format = 30, ...) { # barplot working on enrichKEGG and enrichPathway result objects (default stopped to working xd)
  object <- height
  colorBy <- match.arg(color, c("pvalue", "p.adjust", "qvalue"))
  if (x == "geneRatio" || x == "GeneRatio") {
    x <- "GeneRatio"
  } else if (x == "count" || x == "Count") {
    x <- "Count"
  }
  df <- fortify(object@result[1:showCategory, ], by = x)
  df$Description <- factor(df$Description, levels = rev(df$Description))
  if (colorBy %in% colnames(df)) {
    p <- ggplot(df, aes_string(
      x = x, y = "Description",
      fill = colorBy
    )) +
      theme_dose(font.size) +
      scale_fill_continuous(
        low = "red",
        high = "blue", name = color, guide = guide_colorbar(reverse = TRUE)
      )
  } else {
    p <- ggplot(df, aes_string(
      x = x, y = "Description",
      fill = "Description"
    )) +
      theme_dose(font.size) +
      theme(legend.position = "none")
  }
  p + geom_col() + ggtitle(title) + xlab(NULL) + ylab(NULL)
}



############################################################################


RenameGenesSeurat <- function(obj = ls.Seurat[[i]], newnames = HGNC.updated[[i]]$Suggested.Symbol) { # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  RNA <- obj@assays$RNA
  
  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]] <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]] <- newnames
    if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]] <- newnames
  } else {
    "Unequal gene sets: nrow(RNA) != nrow(newnames)"
  }
  obj@assays$RNA <- RNA
  return(obj)
}




############################################################################


h2m <- function(i) { # human-mouse
  lapply(i, function(j) {
    tmp <- unlist(strsplit(j, ""))
    tmp <- c(tmp[[1]], tolower(tmp[-1]))
    return(paste(tmp, collapse = ""))
  })
}
# as.character(expression(CD62L, CCR7, CXCR4))


#############################################################################


mosPlot <- function(labels,condition,legend=T,title='Percentage of cells in cluster originated from each sample (width of the column = number of cells in cluster)',
                    xlab="Cluster",ylab="Condition", eh=F){
  #### Mosaic plot
  cluster_names_all <- paste0(labels %>% table %>% sort %>% names)
  
  x <- data.frame(cluster_id = labels, condition = condition)
  
  if (eh == T){
    x$condition <- factor(x$condition, levels = c("control","CSI","CSI+aPD1"))
    return(ggplot(data = x) +
             geom_mosaic(aes(x = product(cluster_id), fill=condition), na.rm=F, show.legend = legend) +
             scale_fill_manual(values = rainbow(length(condition %>% table),s=0.8,v=0.8), breaks = c("control","CSI","CSI+aPD1")) +
             labs(x = xlab, y = ylab, title=title) +
             theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold")) +
             geom_hline(yintercept=0.5) +
             geom_text(aes(0,0.5,label = "50 %", vjust = -1)))
  }
  
  if (eh == F){
    #### ggplot geom_mosaic
    x$condition <- factor(x$condition, levels = rev(condition %>% table %>% names))
    mymosplot <- ggplot(data = x) +
      geom_mosaic(aes(x = product(cluster_id), fill=condition), na.rm=F, show.legend = legend) +
      scale_fill_manual(values = sample(mycol,(length(condition %>% table))), 
                        breaks = levels(x$condition),  guide = guide_legend(reverse = TRUE)) +
      labs(x = xlab, y = ylab, title=title) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold"),
            panel.background = element_rect(fill = "gray94", size = 2, linetype = "solid")) +
      geom_hline(yintercept=0) +
      geom_hline(yintercept=1) +
      geom_text(aes(0,-0.02,label = "0 %", vjust = -1)) +
      geom_text(aes(0,0.98,label = "100 %", vjust = -1))
    
    return(
      mymosplot)#theme(legend.position="left"))
  }
}


#############################################################################





#############################################################################


mosPlotNorm <- function(labels,resp,cond="resp"){
  
  if (cond %!in% c("resp","cond")){
    stop("wrong condition")
  }  
  
  cluster_names_all <- paste0(labels %>% table %>% sort %>% names)
  
  x <- data.frame(cluster_id = labels, cond = resp)
  
  ##normalize the response
  x <- table(x)
  
  
  if (cond == "cond"){
    control = sum(x[,"control"])
    noresp = sum(x[,"CSI"])
    resp = sum(x[,"CSI+aPD1"])
    
    x[,"control"]   <-  round(x[,"control"]   / (control/(resp + control + noresp)))
    x[,"CSI"]   <-  round(x[,"CSI"]   / (noresp/(resp + control + noresp)))
    x[,"CSI+aPD1"] <-  round(x[,"CSI+aPD1"] / (resp/(resp + control + noresp)))
    
    df <- as.data.frame(x)
    df = data.frame(x = rep(df$cluster_id, times = df$Freq), y = rep(df$cond, times = df$Freq))
    colnames(df) <- c("cluster","condition")
    
    return(mosPlot(df$cluster, df$condition, legend = T, eh=T, ylab="Condition",title='Proportion of cells in cluster originated from different condition after normalization')
    )
  } 
  
  
  
  if (cond == "resp") {
    control = sum(x[,"control"])
    noresp = sum(x[,"no resp"])
    resp = sum(x[,"resp"])
    
    
    x[,"control"]   <-  round(x[,"control"]   / (control/(resp + control + noresp)))
    x[,"no resp"]   <-  round(x[,"no resp"]   / (noresp/(resp + control + noresp)))
    x[,"resp"] <-  round(x[,"resp"] / (resp/(resp + control + noresp)))
    
    df <- as.data.frame(x)
    df = data.frame(x = rep(df$cluster_id, times = df$Freq), y = rep(df$cond, times = df$Freq))
    colnames(df) <- c("cluster","response")
    
    return(mosPlot(df$cluster, df$response, legend = T, ylab="Response status",title='Proportion of cells in cluster originated from different response status after normalization')
    )
  }
  
  
}


#############################################################################









############################################################################


fp <- function(object, features, dims = c(1, 2), cells = NULL, cols = if (blend) {
  c("lightgrey", "#ff0000", "#00ff00")
} else {
  c("#FDE725FF", "#440154FF")
}, mycol=rev(ryb(50)), pt.size=NULL, pt.s = 0.3, order = FALSE, min.cutoff = NA, max.cutoff = NA,
reduction = NULL, split.by = NULL, keep.scale = "feature",
shape.by = NULL, slot = "data", blend = FALSE, blend.threshold = 0.5,
label = FALSE, label.size = 7, label.color = "midnightblue", label.bg.color="white", bg.r = 0.1, repel = TRUE, max.iter=100,autosize=FALSE,
alpha=0.7, ncol = NULL, coord.fixed = FALSE, by.col = TRUE, sort.cell = NULL, jitter=0,
interactive = FALSE, combine = TRUE, raster = NULL, split.ncol=NULL, split.nrow=NULL) {
  library(patchwork)
  if (!is.null(x = sort.cell)) {
    warning("The sort.cell parameter is being deprecated. Please use the order ",
            "parameter instead for equivalent functionality.",
            call. = FALSE, immediate. = TRUE
    )
    if (isTRUE(x = sort.cell)) {
      order <- sort.cell
    }
  }
  if (interactive) {
    return(IFeaturePlot(
      object = object, feature = features[1],
      dims = dims, reduction = reduction, slot = slot
    ))
  }
  if (!(is.null(x = keep.scale)) && !(keep.scale %in% c(
    "feature",
    "all"
  ))) {
    stop("`keep.scale` must be set to either `feature`, `all`, or NULL")
  }
  no.right <- theme(
    axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(),
    axis.text.y.right = element_blank(), axis.title.y.right = element_text(
      face = "bold.italic",
      size = 14, margin = margin(r = 7)
    )
  )
  reduction <- reduction %||% DefaultDimReduc(object = object)
  if (length(x = dims) != 2 || !is.numeric(x = dims)) {
    stop("'dims' must be a two-length integer vector")
  }
  if (blend && length(x = features) != 2) {
    stop("Blending feature plots only works with two features")
  }
  if (blend) {
    default.colors <- eval(expr = formals(fun = FeaturePlot)$cols)
    cols <- switch(EXPR = as.character(x = length(x = cols)),
                   `0` = {
                     warning("No colors provided, using default colors",
                             call. = FALSE, immediate. = TRUE
                     )
                     default.colors
                   },
                   `1` = {
                     warning("Only one color provided, assuming specified is double-negative and augmenting with default colors",
                             call. = FALSE, immediate. = TRUE
                     )
                     c(cols, default.colors[2:3])
                   },
                   `2` = {
                     warning("Only two colors provided, assuming specified are for features and agumenting with '",
                             default.colors[1], "' for double-negatives",
                             call. = FALSE, immediate. = TRUE
                     )
                     c(default.colors[1], cols)
                   },
                   `3` = cols,
                   {
                     warning("More than three colors provided, using only first three",
                             call. = FALSE, immediate. = TRUE
                     )
                     cols[1:3]
                   }
    )
  }
  if (blend && length(x = cols) != 3) {
    stop("Blending feature plots only works with three colors; first one for negative cells")
  }
  dims <- paste0(Key(object = object[[reduction]]), dims)
  cells <- cells %||% colnames(x = object)
  data <- FetchData(object = object, vars = c(
    dims, "ident",
    features
  ), cells = cells, slot = slot)
  if (ncol(x = data) < 4) {
    stop("None of the requested features were found: ", paste(features,
                                                              collapse = ", "
    ), " in slot ", slot, call. = FALSE)
  } else if (!all(dims %in% colnames(x = data))) {
    stop("The dimensions requested were not found", call. = FALSE)
  }
  features <- colnames(x = data)[4:ncol(x = data)]
  min.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = min(data[
      ,
      feature
    ]), no = cutoff))
  }, cutoff = min.cutoff, feature = features)
  max.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = max(data[
      ,
      feature
    ]), no = cutoff))
  }, cutoff = max.cutoff, feature = features)
  check.lengths <- unique(x = vapply(X = list(
    features, min.cutoff,
    max.cutoff
  ), FUN = length, FUN.VALUE = numeric(length = 1)))
  if (length(x = check.lengths) != 1) {
    stop("There must be the same number of minimum and maximum cuttoffs as there are features")
  }
  brewer.gran <- ifelse(test = length(x = cols) == 1, yes = brewer.pal.info[cols, ]$maxcolors, no = length(x = cols))
  data[, 4:ncol(x = data)] <- sapply(
    X = 4:ncol(x = data),
    FUN = function(index) {
      data.feature <- as.vector(x = data[, index])
      min.use <- Seurat:::SetQuantile(cutoff = min.cutoff[index -
                                                            3], data.feature)
      max.use <- Seurat:::SetQuantile(cutoff = max.cutoff[index -
                                                            3], data.feature)
      data.feature[data.feature < min.use] <- min.use
      data.feature[data.feature > max.use] <- max.use
      if (brewer.gran == 2) {
        return(data.feature)
      }
      data.cut <- if (all(data.feature == 0)) {
        0
      } else {
        as.numeric(x = as.factor(x = cut(
          x = as.numeric(x = data.feature),
          breaks = brewer.gran
        )))
      }
      return(data.cut)
    }
  )
  colnames(x = data)[4:ncol(x = data)] <- features
  rownames(x = data) <- cells
  data$split <- if (is.null(x = split.by)) {
    RandomName()
  } else {
    switch(EXPR = split.by,
           ident = Idents(object = object)[cells,
                                           drop = TRUE
           ],
           object[[split.by, drop = TRUE]][cells,
                                           drop = TRUE
           ]
    )
  }
  if (!is.factor(x = data$split)) {
    data$split <- factor(x = data$split)
  }
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  plots <- vector(mode = "list", length = ifelse(test = blend,
                                                 yes = 4, no = length(x = features) * length(x = levels(x = data$split))
  ))
  xlims <- c(floor(x = min(data[, dims[1]])), ceiling(x = max(data[
    ,
    dims[1]
  ])))
  ylims <- c(floor(min(data[, dims[2]])), ceiling(x = max(data[
    ,
    dims[2]
  ])))
  if (blend) {
    ncol <- 4
    color.matrix <- Seurat:::BlendMatrix(
      two.colors = cols[2:3], col.threshold = blend.threshold,
      negative.color = cols[1]
    )
    cols <- cols[2:3]
    colors <- list(
      color.matrix[, 1], color.matrix[1, ],
      as.vector(x = color.matrix)
    )
  }
  for (i in 1:length(x = levels(x = data$split))) {
    ident <- levels(x = data$split)[i]
    data.plot <- data[as.character(x = data$split) == ident, ,
                      drop = FALSE
    ]
    if (blend) {
      features <- features[1:2]
      no.expression <- features[colMeans(x = data.plot[
        ,
        features
      ]) == 0]
      if (length(x = no.expression) != 0) {
        stop("The following features have no value: ",
             paste(no.expression, collapse = ", "),
             call. = FALSE
        )
      }
      data.plot <- cbind(
        data.plot[, c(dims, "ident")],
        Seurat:::BlendExpression(data = data.plot[, features[1:2]])
      )
      features <- colnames(x = data.plot)[4:ncol(x = data.plot)]
    }
    for (j in 1:length(x = features)) {
      feature <- features[j]
      if (blend) {
        cols.use <- as.numeric(x = as.character(x = data.plot[
          ,
          feature
        ])) + 1
        cols.use <- colors[[j]][sort(x = unique(x = cols.use))]
      } else {
        cols.use <- NULL
      }
      data.single <- data.plot[, c(
        dims, "ident", feature,
        shape.by
      )]
      plot <- mySingleDimPlot(
        data = data.single, dims = dims, alpha=alpha,
        col.by = feature, order = order, pt.size = pt.size,
        cols = cols.use, shape.by = shape.by, label = FALSE,
        raster = raster, my.pt.size = pt.s, autosize = autosize,
        jitter = jitter
      ) + scale_x_continuous(limits = xlims) +
        scale_y_continuous(limits = ylims) + theme_cowplot() +
        Seurat:::CenterTitle()
      if (label) {
        plot <- fplabelClustersWithOutilnedText(plot = plot, id = "ident", repel = repel, size = label.size, color = label.color, bg.color=label.bg.color, bg.r=bg.r,max.iter=max.iter)
      }
      if (length(x = levels(x = data$split)) > 1) {
        plot <- plot + theme(panel.border = element_rect(
          fill = NA,
          colour = "black"
        ))
        plot <- plot + if (i == 1) {
          labs(title = feature)
        } else {
          labs(title = NULL)
        }
        if (j == length(x = features) && !blend) {
          suppressMessages(expr = plot <- plot + scale_y_continuous(
            sec.axis = dup_axis(name = ident),
            limits = ylims
          ) + no.right)
        }
        if (j != 1) {
          plot <- plot + theme(
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(), axis.text.y = element_blank(),
            axis.title.y.left = element_blank()
          )
        }
        if (i != length(x = levels(x = data$split))) {
          plot <- plot + theme(
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank(), axis.text.x = element_blank(),
            axis.title.x = element_blank()
          )
        }
      } else {
        plot <- plot + labs(title = feature)
      }
      if (!blend) {
        plot <- plot + guides(color = NULL)
        cols.grad <- cols
        if (length(x = cols) == 1) {
          plot <- plot + scale_color_brewer(palette = cols)
        } else if (length(x = cols) > 1) {
          unique.feature.exp <- unique(data.plot[, feature])
          if (length(unique.feature.exp) == 1) {
            warning(
              "All cells have the same value (",
              unique.feature.exp, ") of ", feature, "."
            )
            if (unique.feature.exp == 0) {
              cols.grad <- cols[1]
            } else {
              cols.grad <- cols
            }
          }
          plot <- suppressMessages(expr = plot + scale_color_gradientn(
            colors = cols.grad,
            guide = "colorbar"
          ))
        }
      }
      if (!(is.null(x = keep.scale)) && keep.scale == "feature" &&
          !blend) {
        max.feature.value <- max(data.single[, feature])
        min.feature.value <- min(data.single[, feature])
        plot <- suppressMessages(plot & scale_color_gradientn(
          colors = mycol,
          limits = c(min.feature.value, max.feature.value)
        ))
      }
      if (coord.fixed) {
        plot <- plot + coord_fixed()
      }
      plot <- plot
      plots[[(length(x = features) * (i - 1)) + j]] <- plot
    }
  }
  if (blend) {
    blend.legend <- Seurat:::BlendMap(color.matrix = color.matrix)
    for (ii in 1:length(x = levels(x = data$split))) {
      suppressMessages(expr = plots <- append(
        x = plots,
        values = list(blend.legend + scale_y_continuous(
          sec.axis = dup_axis(name = ifelse(test = length(x = levels(x = data$split)) >
                                              1, yes = levels(x = data$split)[ii], no = "")),
          expand = c(0, 0)
        ) + labs(
          x = features[1], y = features[2],
          title = if (ii == 1) {
            paste("Color threshold:", blend.threshold)
          } else {
            NULL
          }
        ) + no.right), after = 4 * ii - 1
      ))
    }
  }
  plots <- Filter(f = Negate(f = is.null), x = plots)
  if (is.null(x = ncol)) {
    ncol <- 2
    if (length(x = features) == 1) {
      ncol <- 1
    }
    if (length(x = features) > 6) {
      ncol <- 3
    }
    if (length(x = features) > 9) {
      ncol <- 4
    }
  }
  ncol <- ifelse(test = is.null(x = split.by) || blend, yes = ncol,
                 no = length(x = features))
  legend <- if (blend) {
    "none"
  } else {
    split.by %iff% "none"
  }
  if (combine) {
    if (by.col && !is.null(x = split.by) && !blend) {
      plots <- lapply(X = plots, FUN = function(x) {
        return(suppressMessages(expr = x + theme_cowplot() +
                                  ggtitle("") + scale_y_continuous(
                                    sec.axis = dup_axis(name = ""),
                                    limits = ylims
                                  ) + no.right))
      })
      nsplits <- length(x = levels(x = data$split))
      idx <- 1
      for (i in (length(x = features) * (nsplits - 1) +
                 1):(length(x = features) * nsplits)) {
        plots[[i]] <- suppressMessages(expr = plots[[i]] +
                                         scale_y_continuous(
                                           sec.axis = dup_axis(name = features[[idx]]),
                                           limits = ylims
                                         ) + no.right)
        idx <- idx + 1
      }
      idx <- 1
      for (i in which(x = 1:length(x = plots) %% length(x = features) ==
                      1)) {
        plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]]) +
          theme(plot.title = element_text(hjust = 0.5))
        idx <- idx + 1
      }
      idx <- 1
      if (length(x = features) == 1) {
        for (i in 1:length(x = plots)) {
          plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]]) +
            theme(plot.title = element_text(hjust = 0.5))
          idx <- idx + 1
        }
        ncol <- 1
        nrow <- nsplits
      } else {
        if (is.null(split.ncol)){
          nrow <- split.by %iff% length(x = levels(x = data$split))
        }
        else {
          nrow <- split.ncol
          ncol <- split.nrow
        }
      }
      plots <- plots[c(do.call(what = rbind, args = split(
        x = 1:length(x = plots),
        f = ceiling(x = seq_along(along.with = 1:length(x = plots)) / length(x = features))
      )))]
      plots <- wrap_plots(plots, ncol = nrow, nrow = ncol)
      if (!is.null(x = legend) && legend == "none") {
        plots <- plots & NoLegend()
      }
    } else {
      plots <- wrap_plots(plots, ncol = ncol, nrow = split.by %iff%
                            length(x = levels(x = data$split)))
    }
    if (!is.null(x = legend) && legend == "none") {
      plots <- plots & NoLegend()
    }
    if (!(is.null(x = keep.scale)) && keep.scale == "all" &&
        !blend) {
      max.feature.value <- max(data.plot[, features])
      min.feature.value <- min(data.plot[, features])
      plots <- suppressMessages(plots & scale_color_gradientn(
        colors = mycol,
        limits = c(min.feature.value, max.feature.value)
      ))
    }
  }
  return(plots)
}







############################################################################


fp2 <- function(object, features, dims = c(1, 2), cells = NULL, cols = if (blend) {
  c("lightgrey", "#ff0000", "#00ff00")
} else {
  c("#FDE725FF", "#440154FF")
}, pt.size=NULL, pt.s = 1, order = FALSE, min.cutoff = NA, max.cutoff = NA,
reduction = NULL, split.by = NULL, keep.scale = "feature",
shape.by = NULL, slot = "data", blend = FALSE, blend.threshold = 0.5,
label = FALSE, label.size = 7, label.color = "midnightblue", label.bg.color="white", bg.r = 0.1, repel = TRUE, max.iter=1,autosize=FALSE,
alpha=0.7, ncol = NULL, coord.fixed = FALSE, by.col = TRUE, sort.cell = NULL,
interactive = FALSE, combine = TRUE, raster = NULL, split.ncol=NULL, split.nrow=NULL) {
  library(patchwork)
  if (!is.null(x = sort.cell)) {
    warning("The sort.cell parameter is being deprecated. Please use the order ",
            "parameter instead for equivalent functionality.",
            call. = FALSE, immediate. = TRUE
    )
    if (isTRUE(x = sort.cell)) {
      order <- sort.cell
    }
  }
  if (interactive) {
    return(IFeaturePlot(
      object = object, feature = features[1],
      dims = dims, reduction = reduction, slot = slot
    ))
  }
  if (!(is.null(x = keep.scale)) && !(keep.scale %in% c(
    "feature",
    "all"
  ))) {
    stop("`keep.scale` must be set to either `feature`, `all`, or NULL")
  }
  no.right <- theme(
    axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(),
    axis.text.y.right = element_blank(), axis.title.y.right = element_text(
      face = "bold.italic",
      size = 14, margin = margin(r = 7)
    )
  )
  reduction <- reduction %||% DefaultDimReduc(object = object)
  if (length(x = dims) != 2 || !is.numeric(x = dims)) {
    stop("'dims' must be a two-length integer vector")
  }
  if (blend && length(x = features) != 2) {
    stop("Blending feature plots only works with two features")
  }
  if (blend) {
    default.colors <- eval(expr = formals(fun = FeaturePlot)$cols)
    cols <- switch(EXPR = as.character(x = length(x = cols)),
                   `0` = {
                     warning("No colors provided, using default colors",
                             call. = FALSE, immediate. = TRUE
                     )
                     default.colors
                   },
                   `1` = {
                     warning("Only one color provided, assuming specified is double-negative and augmenting with default colors",
                             call. = FALSE, immediate. = TRUE
                     )
                     c(cols, default.colors[2:3])
                   },
                   `2` = {
                     warning("Only two colors provided, assuming specified are for features and agumenting with '",
                             default.colors[1], "' for double-negatives",
                             call. = FALSE, immediate. = TRUE
                     )
                     c(default.colors[1], cols)
                   },
                   `3` = cols,
                   {
                     warning("More than three colors provided, using only first three",
                             call. = FALSE, immediate. = TRUE
                     )
                     cols[1:3]
                   }
    )
  }
  if (blend && length(x = cols) != 3) {
    stop("Blending feature plots only works with three colors; first one for negative cells")
  }
  dims <- paste0(Key(object = object[[reduction]]), dims)
  cells <- cells %||% colnames(x = object)
  data <- FetchData(object = object, vars = c(
    dims, "ident",
    features
  ), cells = cells, slot = slot)
  if (ncol(x = data) < 4) {
    stop("None of the requested features were found: ", paste(features,
                                                              collapse = ", "
    ), " in slot ", slot, call. = FALSE)
  } else if (!all(dims %in% colnames(x = data))) {
    stop("The dimensions requested were not found", call. = FALSE)
  }
  features <- colnames(x = data)[4:ncol(x = data)]
  min.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = min(data[
      ,
      feature
    ]), no = cutoff))
  }, cutoff = min.cutoff, feature = features)
  max.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = max(data[
      ,
      feature
    ]), no = cutoff))
  }, cutoff = max.cutoff, feature = features)
  check.lengths <- unique(x = vapply(X = list(
    features, min.cutoff,
    max.cutoff
  ), FUN = length, FUN.VALUE = numeric(length = 1)))
  if (length(x = check.lengths) != 1) {
    stop("There must be the same number of minimum and maximum cuttoffs as there are features")
  }
  brewer.gran <- ifelse(test = length(x = cols) == 1, yes = brewer.pal.info[cols, ]$maxcolors, no = length(x = cols))
  data[, 4:ncol(x = data)] <- sapply(
    X = 4:ncol(x = data),
    FUN = function(index) {
      data.feature <- as.vector(x = data[, index])
      min.use <- Seurat:::SetQuantile(cutoff = min.cutoff[index -
                                                            3], data.feature)
      max.use <- Seurat:::SetQuantile(cutoff = max.cutoff[index -
                                                            3], data.feature)
      data.feature[data.feature < min.use] <- min.use
      data.feature[data.feature > max.use] <- max.use
      if (brewer.gran == 2) {
        return(data.feature)
      }
      data.cut <- if (all(data.feature == 0)) {
        0
      } else {
        as.numeric(x = as.factor(x = cut(
          x = as.numeric(x = data.feature),
          breaks = brewer.gran
        )))
      }
      return(data.cut)
    }
  )
  colnames(x = data)[4:ncol(x = data)] <- features
  rownames(x = data) <- cells
  data$split <- if (is.null(x = split.by)) {
    RandomName()
  } else {
    switch(EXPR = split.by,
           ident = Idents(object = object)[cells,
                                           drop = TRUE
           ],
           object[[split.by, drop = TRUE]][cells,
                                           drop = TRUE
           ]
    )
  }
  if (!is.factor(x = data$split)) {
    data$split <- factor(x = data$split)
  }
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  plots <- vector(mode = "list", length = ifelse(test = blend,
                                                 yes = 4, no = length(x = features) * length(x = levels(x = data$split))
  ))
  xlims <- c(floor(x = min(data[, dims[1]])), ceiling(x = max(data[
    ,
    dims[1]
  ])))
  ylims <- c(floor(min(data[, dims[2]])), ceiling(x = max(data[
    ,
    dims[2]
  ])))
  if (blend) {
    ncol <- 4
    color.matrix <- Seurat:::BlendMatrix(
      two.colors = cols[2:3], col.threshold = blend.threshold,
      negative.color = cols[1]
    )
    cols <- cols[2:3]
    colors <- list(
      color.matrix[, 1], color.matrix[1, ],
      as.vector(x = color.matrix)
    )
  }
  for (i in 1:length(x = levels(x = data$split))) {
    ident <- levels(x = data$split)[i]
    data.plot <- data[as.character(x = data$split) == ident, ,
                      drop = FALSE
    ]
    if (blend) {
      features <- features[1:2]
      no.expression <- features[colMeans(x = data.plot[
        ,
        features
      ]) == 0]
      if (length(x = no.expression) != 0) {
        stop("The following features have no value: ",
             paste(no.expression, collapse = ", "),
             call. = FALSE
        )
      }
      data.plot <- cbind(
        data.plot[, c(dims, "ident")],
        Seurat:::BlendExpression(data = data.plot[, features[1:2]])
      )
      features <- colnames(x = data.plot)[4:ncol(x = data.plot)]
    }
    for (j in 1:length(x = features)) {
      feature <- features[j]
      if (blend) {
        cols.use <- as.numeric(x = as.character(x = data.plot[
          ,
          feature
        ])) + 1
        cols.use <- colors[[j]][sort(x = unique(x = cols.use))]
      } else {
        cols.use <- NULL
      }
      data.single <- data.plot[, c(
        dims, "ident", feature,
        shape.by
      )]
      plot <- mySingleDimPlot(
        data = data.single, dims = dims, alpha=alpha,
        col.by = feature, order = order, pt.size = pt.size,
        cols = cols.use, shape.by = shape.by, label = FALSE,
        raster = raster, my.pt.size = pt.s, autosize = autosize,
        jitter = jitter
      ) + scale_x_continuous(limits = xlims) +
        scale_y_continuous(limits = ylims) + theme_cowplot() +
        Seurat:::CenterTitle()
      if (label) {
        plot <- fplabelClustersWithOutilnedText(plot = plot, id = "ident", repel = repel, size = label.size, color = label.color, bg.color=label.bg.color, bg.r=bg.r,max.iter=max.iter)
      }
      if (length(x = levels(x = data$split)) > 1) {
        plot <- plot + theme(panel.border = element_rect(
          fill = NA,
          colour = "black"
        ))
        plot <- plot + if (i == 1) {
          labs(title = feature)
        } else {
          labs(title = NULL)
        }
        if (j == length(x = features) && !blend) {
          suppressMessages(expr = plot <- plot + scale_y_continuous(
            sec.axis = dup_axis(name = ident),
            limits = ylims
          ) + no.right)
        }
        if (j != 1) {
          plot <- plot + theme(
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(), axis.text.y = element_blank(),
            axis.title.y.left = element_blank()
          )
        }
        if (i != length(x = levels(x = data$split))) {
          plot <- plot + theme(
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank(), axis.text.x = element_blank(),
            axis.title.x = element_blank()
          )
        }
      } else {
        plot <- plot + labs(title = feature)
      }
      if (!blend) {
        plot <- plot + guides(color = NULL)
        cols.grad <- cols
        if (length(x = cols) == 1) {
          plot <- plot + scale_color_brewer(palette = cols)
        } else if (length(x = cols) > 1) {
          unique.feature.exp <- unique(data.plot[, feature])
          if (length(unique.feature.exp) == 1) {
            warning(
              "All cells have the same value (",
              unique.feature.exp, ") of ", feature, "."
            )
            if (unique.feature.exp == 0) {
              cols.grad <- cols[1]
            } else {
              cols.grad <- cols
            }
          }
          plot <- suppressMessages(expr = plot + scale_color_gradientn(
            colors = cols.grad,
            guide = "colorbar"
          ))
        }
      }
      if (!(is.null(x = keep.scale)) && keep.scale == "feature" &&
          !blend) {
        max.feature.value <- max(data.single[, feature])
        min.feature.value <- min(data.single[, feature])
        plot <- suppressMessages(plot & scale_color_gradientn(
          colors = mycol,
          limits = c(min.feature.value, max.feature.value)
        ))
      }
      if (coord.fixed) {
        plot <- plot + coord_fixed()
      }
      plot <- plot
      plots[[(length(x = features) * (i - 1)) + j]] <- plot
    }
  }
  if (blend) {
    blend.legend <- Seurat:::BlendMap(color.matrix = color.matrix)
    for (ii in 1:length(x = levels(x = data$split))) {
      suppressMessages(expr = plots <- append(
        x = plots,
        values = list(blend.legend + scale_y_continuous(
          sec.axis = dup_axis(name = ifelse(test = length(x = levels(x = data$split)) >
                                              1, yes = levels(x = data$split)[ii], no = "")),
          expand = c(0, 0)
        ) + labs(
          x = features[1], y = features[2],
          title = if (ii == 1) {
            paste("Color threshold:", blend.threshold)
          } else {
            NULL
          }
        ) + no.right), after = 4 * ii - 1
      ))
    }
  }
  plots <- Filter(f = Negate(f = is.null), x = plots)
  if (is.null(x = ncol)) {
    ncol <- 2
    if (length(x = features) == 1) {
      ncol <- 1
    }
    if (length(x = features) > 6) {
      ncol <- 3
    }
    if (length(x = features) > 9) {
      ncol <- 4
    }
  }
  ncol <- ifelse(test = is.null(x = split.by) || blend, yes = ncol,
                 no = length(x = features))
  legend <- if (blend) {
    "none"
  } else {
    split.by %iff% "none"
  }
  if (combine) {
    if (by.col && !is.null(x = split.by) && !blend) {
      plots <- lapply(X = plots, FUN = function(x) {
        return(suppressMessages(expr = x + theme_cowplot() +
                                  ggtitle("") + scale_y_continuous(
                                    sec.axis = dup_axis(name = ""),
                                    limits = ylims
                                  ) + no.right))
      })
      nsplits <- length(x = levels(x = data$split))
      idx <- 1
      for (i in (length(x = features) * (nsplits - 1) +
                 1):(length(x = features) * nsplits)) {
        plots[[i]] <- suppressMessages(expr = plots[[i]] +
                                         scale_y_continuous(
                                           sec.axis = dup_axis(name = features[[idx]]),
                                           limits = ylims
                                         ) + no.right)
        idx <- idx + 1
      }
      idx <- 1
      for (i in which(x = 1:length(x = plots) %% length(x = features) ==
                      1)) {
        plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]]) +
          theme(plot.title = element_text(hjust = 0.5))
        idx <- idx + 1
      }
      idx <- 1
      if (length(x = features) == 1) {
        for (i in 1:length(x = plots)) {
          plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]]) +
            theme(plot.title = element_text(hjust = 0.5))
          idx <- idx + 1
        }
        ncol <- 1
        nrow <- nsplits
      } else {
        if (is.null(split.ncol)){
          nrow <- split.by %iff% length(x = levels(x = data$split))
        }
        else {
          nrow <- split.ncol
          ncol <- split.nrow
        }
      }
      plots <- plots[c(do.call(what = rbind, args = split(
        x = 1:length(x = plots),
        f = ceiling(x = seq_along(along.with = 1:length(x = plots)) / length(x = features))
      )))]
      plots <- wrap_plots(plots, ncol = nrow, nrow = ncol)
      if (!is.null(x = legend) && legend == "none") {
        plots <- plots & NoLegend()
      }
    } else {
      plots <- wrap_plots(plots, ncol = ncol, nrow = split.by %iff%
                            length(x = levels(x = data$split)))
    }
    if (!is.null(x = legend) && legend == "none") {
      plots <- plots & NoLegend()
    }
    if (!(is.null(x = keep.scale)) && keep.scale == "all" &&
        !blend) {
      max.feature.value <- max(data.plot[, features])
      min.feature.value <- min(data.plot[, features])
      plots <- suppressMessages(plots & scale_color_gradientn(
        colors = mycol,
        limits = c(min.feature.value, max.feature.value)
      ))
    }
  }
  return(plots)
}









############################################################################

featurePlotMultiple <- function(object.list, features.to.plot,
                                feature.labels = NULL,
                                filename,
                                assay = NULL,
                                ncol = 2,
                                minCutOff = NA,
                                maxCutOff = NA,
                                colorPalette = brewer.pal(11, "RdYlBu")[11:1]) {
  if (is.null(feature.labels)) {
    feature.labels <- features.to.plot
  }
  if (!is.null(assay)) {
    for (i in 1:length(object.list)) {
      DefaultAssay(object.list[[i]]) <- assay
    }
  }
  
  p.list <- lapply(object.list, function(object.to.plot) {
    lapply(seq_along(features.to.plot), function(i) {
      if (features.to.plot[i] %in% rownames(object.to.plot)) {
        FeaturePlot(object.to.plot,
                    features = features.to.plot[i],
                    pt.size = 0.9, min.cutoff = minCutOff, max.cutoff = maxCutOff
        ) +
          scale_colour_gradientn(colours = colorPalette) +
          ggtitle(feature.labels[i])
      } else {
        ggplot() +
          theme_void()
      }
    })
  })
  
  png(file = filename, width = 3840 / 2, height = 2400 / 2)
  print(CombinePlots(plots = unlist(p.list, recursive = F), ncol = ncol))
  dev.off()
}




############################################################################





splitDataFrame <- function(df, listToSplitBy) {
  n <- 0
  marklist <- list()
  for (i in listToSplitBy) {
    n <- n + 1
    marklist[[n]] <- markers[df$cluster == i, ]
  }
  return(marklist)
}









############################################################################


fplabelClustersWithOutilnedText <- function(plot, id, clusters = NULL, labels = NULL, split.by = NULL, repel = TRUE, box = FALSE, geom = "GeomPoint", position = "median", color = "white", bg.color = "midnightblue", bg.r=0.1,
                                            max.iter=1,
                                            ...) {
  xynames <- unlist(
    x = Seurat:::GetXYAesthetics(plot = plot, geom = geom),
    use.names = TRUE
  )
  if (!id %in% colnames(x = plot$data)) {
    stop("Cannot find variable ", id, " in plotting data")
  }
  if (!is.null(x = split.by) && !split.by %in% colnames(x = plot$data)) {
    warning("Cannot find splitting variable ", id, " in plotting data")
    split.by <- NULL
  }
  data <- plot$data[, c(xynames, id, split.by)]
  possible.clusters <- as.character(x = na.omit(object = unique(x = data[
    ,
    id
  ])))
  groups <- clusters %||% as.character(x = na.omit(object = unique(x = data[
    ,
    id
  ])))
  if (any(!groups %in% possible.clusters)) {
    stop("The following clusters were not found: ", paste(groups[!groups %in%
                                                                   possible.clusters], collapse = ","))
  }
  pb <- ggplot_build(plot = plot)
  if (geom == "GeomSpatial") {
    xrange.save <- layer_scales(plot = plot)$x$range$range
    yrange.save <- layer_scales(plot = plot)$y$range$range
    data[, xynames["y"]] <- max(data[, xynames["y"]]) - data[
      ,
      xynames["y"]
    ] + min(data[, xynames["y"]])
    if (!pb$plot$plot_env$crop) {
      y.transform <- c(0, nrow(x = pb$plot$plot_env$image)) -
        pb$layout$panel_params[[1]]$y.range
      data[, xynames["y"]] <- data[, xynames["y"]] + sum(y.transform)
    }
  }
  data <- cbind(data, color = pb$data[[1]][[1]])
  labels.loc <- lapply(X = groups, FUN = function(group) {
    data.use <- data[data[, id] == group, , drop = FALSE]
    data.medians <- if (!is.null(x = split.by)) {
      do.call(what = "rbind", args = lapply(X = unique(x = data.use[
        ,
        split.by
      ]), FUN = function(split) {
        medians <- apply(
          X = data.use[data.use[, split.by] ==
                         split, xynames, drop = FALSE], MARGIN = 2,
          FUN = median, na.rm = TRUE
        )
        medians <- as.data.frame(x = t(x = medians))
        medians[, split.by] <- split
        return(medians)
      }))
    } else {
      as.data.frame(x = t(x = apply(X = data.use[, xynames,
                                                 drop = FALSE
      ], MARGIN = 2, FUN = median, na.rm = TRUE)))
    }
    data.medians[, id] <- group
    data.medians$color <- data.use$color[1]
    return(data.medians)
  })
  if (position == "nearest") {
    labels.loc <- lapply(X = labels.loc, FUN = function(x) {
      group.data <- data[as.character(x = data[, id]) ==
                           as.character(x[3]), ]
      nearest.point <- nn2(data = group.data[, 1:2], query = as.matrix(x = x[c(
        1,
        2
      )]), k = 1)$nn.idx
      x[1:2] <- group.data[nearest.point, 1:2]
      return(x)
    })
  }
  labels.loc <- do.call(what = "rbind", args = labels.loc)
  labels.loc[, id] <- factor(x = labels.loc[, id], levels = levels(data[
    ,
    id
  ]))
  labels <- labels %||% groups
  if (length(x = unique(x = labels.loc[, id])) != length(x = labels)) {
    stop(
      "Length of labels (", length(x = labels), ") must be equal to the number of clusters being labeled (",
      length(x = labels.loc), ")."
    )
  }
  names(x = labels) <- groups
  for (group in groups) {
    labels.loc[labels.loc[, id] == group, id] <- labels[group]
  }
  if (box) {
    geom.use <- ifelse(test = repel, yes = geom_label_repel,
                       no = geom_label
    )
    plot <- plot + geom.use(
      data = labels.loc, mapping = aes_string(
        x = xynames["x"],
        y = xynames["y"], label = id, fill = id
      ), show.legend = FALSE,
      ...
    ) + scale_fill_manual(values = labels.loc$color[order(labels.loc[
      ,
      id
    ])])
  } else {
    geom.use <- ifelse(test = repel, yes = geom_text_repel, no = geom_text)
    plot <- plot + geom.use(
      data = labels.loc, mapping = aes_string(
        x = xynames["x"],
        y = xynames["y"], label = id
      ), max.iter=max.iter,
      color = color, bg.color = bg.color, bg.r = bg.r, show.legend = FALSE,
      ...
    )
  }
  if (geom == "GeomSpatial") {
    plot <- suppressMessages(expr = plot + coord_fixed(
      xlim = xrange.save,
      ylim = yrange.save
    ))
  }
  return(plot)
}








############################################################################


dplabelClustersWithOutilnedText <- function(plot, id, clusters = NULL, labels = NULL, split.by = NULL, repel = TRUE, box = FALSE, geom = "GeomPoint", position = "median", color = "white", bg.r = 0.1, max.iter=1, mylegend=mylegend,
                                            ...) {
  xynames <- unlist(
    x = Seurat:::GetXYAesthetics(plot = plot, geom = geom),
    use.names = TRUE
  )
  if (!id %in% colnames(x = plot$data)) {
    stop("Cannot find variable ", id, " in plotting data")
  }
  if (!is.null(x = split.by) && !split.by %in% colnames(x = plot$data)) {
    warning("Cannot find splitting variable ", id, " in plotting data")
    split.by <- NULL
  }
  data <- plot$data[, c(xynames, id, split.by)]
  possible.clusters <- as.character(x = na.omit(object = unique(x = data[
    ,
    id
  ])))
  groups <- clusters %||% as.character(x = na.omit(object = unique(x = data[
    ,
    id
  ])))
  if (any(!groups %in% possible.clusters)) {
    stop("The following clusters were not found: ", paste(groups[!groups %in%
                                                                   possible.clusters], collapse = ","))
  }
  pb <- ggplot_build(plot = plot)
  if (geom == "GeomSpatial") {
    xrange.save <- layer_scales(plot = plot)$x$range$range
    yrange.save <- layer_scales(plot = plot)$y$range$range
    data[, xynames["y"]] <- max(data[, xynames["y"]]) - data[
      ,
      xynames["y"]
    ] + min(data[, xynames["y"]])
    if (!pb$plot$plot_env$crop) {
      y.transform <- c(0, nrow(x = pb$plot$plot_env$image)) -
        pb$layout$panel_params[[1]]$y.range
      data[, xynames["y"]] <- data[, xynames["y"]] + sum(y.transform)
    }
  }
  data <- cbind(data, color = pb$data[[1]][[1]])
  labels.loc <- lapply(X = groups, FUN = function(group) {
    data.use <- data[data[, id] == group, , drop = FALSE]
    data.medians <- if (!is.null(x = split.by)) {
      do.call(what = "rbind", args = lapply(X = unique(x = data.use[
        ,
        split.by
      ]), FUN = function(split) {
        medians <- apply(
          X = data.use[data.use[, split.by] ==
                         split, xynames, drop = FALSE], MARGIN = 2,
          FUN = median, na.rm = TRUE
        )
        medians <- as.data.frame(x = t(x = medians))
        medians[, split.by] <- split
        return(medians)
      }))
    } else {
      as.data.frame(x = t(x = apply(X = data.use[, xynames,
                                                 drop = FALSE
      ], MARGIN = 2, FUN = median, na.rm = TRUE)))
    }
    data.medians[, id] <- group
    data.medians$color <- data.use$color[1]
    return(data.medians)
  })
  if (position == "nearest") {
    labels.loc <- lapply(X = labels.loc, FUN = function(x) {
      group.data <- data[as.character(x = data[, id]) ==
                           as.character(x[3]), ]
      nearest.point <- nn2(data = group.data[, 1:2], query = as.matrix(x = x[c(
        1,
        2
      )]), k = 1)$nn.idx
      x[1:2] <- group.data[nearest.point, 1:2]
      return(x)
    })
  }
  labels.loc <- do.call(what = "rbind", args = labels.loc)
  labels.loc[, id] <- factor(x = labels.loc[, id], levels = levels(data[
    ,
    id
  ]))
  labels <- labels %||% groups
  if (length(x = unique(x = labels.loc[, id])) != length(x = labels)) {
    stop(
      "Length of labels (", length(x = labels), ") must be equal to the number of clusters being labeled (",
      length(x = labels.loc), ")."
    )
  }
  names(x = labels) <- groups
  for (group in groups) {
    labels.loc[labels.loc[, id] == group, id] <- labels[group]
  }
  if (box) {
    geom.use <- ifelse(test = repel, yes = geom_label_repel,
                       no = geom_label
    )
    plot <- plot + geom.use(
      data = labels.loc, mapping = aes_string(
        x = xynames["x"],
        y = xynames["y"], label = id, fill = id
      ), show.legend = FALSE,
      ...
    ) + scale_fill_manual(values = labels.loc$color[order(labels.loc[
      ,
      id,
    ])])
  } else {
    geom.use <- ifelse(test = repel, yes = geom_text_repel, no = geom_text)
    labels.loc <- labels.loc[order(labels.loc[,3]),]
    plot <- plot + geom.use(
      data = labels.loc, mapping = aes_string(
        x = xynames["x"],
        y = xynames["y"], label = id
      ), max.iter=max.iter,
      color = color, bg.color = labels.loc$color[order(labels.loc[,id,])], bg.r = bg.r, show.legend = FALSE, ...
    )
    
    
  }
  if (geom == "GeomSpatial") {
    plot <- suppressMessages(expr = plot + coord_fixed(
      xlim = xrange.save,
      ylim = yrange.save
    ))
  }
  if (mylegend == T){
    return(plot)
  }
  else {
    return(plot + NoLegend())
  }
}





############################################################################







############################################################################


dplabelClustersWithOutilnedText2 <- function(plot, id, clusters = NULL, labels = NULL, split.by = NULL, repel = TRUE, box = FALSE, geom = "GeomPoint", position = "median", color = "white", bg.r = 0.1, max.iter=1,
                                             ...) {
  xynames <- unlist(
    x = Seurat:::GetXYAesthetics(plot = plot, geom = geom),
    use.names = TRUE
  )
  if (!id %in% colnames(x = plot$data)) {
    stop("Cannot find variable ", id, " in plotting data")
  }
  if (!is.null(x = split.by) && !split.by %in% colnames(x = plot$data)) {
    warning("Cannot find splitting variable ", id, " in plotting data")
    split.by <- NULL
  }
  data <- plot$data[, c(xynames, id, split.by)]
  possible.clusters <- as.character(x = na.omit(object = unique(x = data[
    ,
    id
  ])))
  groups <- clusters %||% as.character(x = na.omit(object = unique(x = data[
    ,
    id
  ])))
  if (any(!groups %in% possible.clusters)) {
    stop("The following clusters were not found: ", paste(groups[!groups %in%
                                                                   possible.clusters], collapse = ","))
  }
  pb <- ggplot_build(plot = plot)
  if (geom == "GeomSpatial") {
    xrange.save <- layer_scales(plot = plot)$x$range$range
    yrange.save <- layer_scales(plot = plot)$y$range$range
    data[, xynames["y"]] <- max(data[, xynames["y"]]) - data[
      ,
      xynames["y"]
    ] + min(data[, xynames["y"]])
    if (!pb$plot$plot_env$crop) {
      y.transform <- c(0, nrow(x = pb$plot$plot_env$image)) -
        pb$layout$panel_params[[1]]$y.range
      data[, xynames["y"]] <- data[, xynames["y"]] + sum(y.transform)
    }
  }
  data <- cbind(data, color = pb$data[[1]][[1]])
  labels.loc <- lapply(X = groups, FUN = function(group) {
    data.use <- data[data[, id] == group, , drop = FALSE]
    data.medians <- if (!is.null(x = split.by)) {
      do.call(what = "rbind", args = lapply(X = unique(x = data.use[
        ,
        split.by
      ]), FUN = function(split) {
        medians <- apply(
          X = data.use[data.use[, split.by] ==
                         split, xynames, drop = FALSE], MARGIN = 2,
          FUN = median, na.rm = TRUE
        )
        medians <- as.data.frame(x = t(x = medians))
        medians[, split.by] <- split
        return(medians)
      }))
    } else {
      as.data.frame(x = t(x = apply(X = data.use[, xynames,
                                                 drop = FALSE
      ], MARGIN = 2, FUN = median, na.rm = TRUE)))
    }
    data.medians[, id] <- group
    data.medians$color <- data.use$color[1]
    return(data.medians)
  })
  if (position == "nearest") {
    labels.loc <- lapply(X = labels.loc, FUN = function(x) {
      group.data <- data[as.character(x = data[, id]) ==
                           as.character(x[3]), ]
      nearest.point <- nn2(data = group.data[, 1:2], query = as.matrix(x = x[c(
        1,
        2
      )]), k = 1)$nn.idx
      x[1:2] <- group.data[nearest.point, 1:2]
      return(x)
    })
  }
  labels.loc <- do.call(what = "rbind", args = labels.loc)
  labels.loc[, id] <- factor(x = labels.loc[, id], levels = levels(data[
    ,
    id
  ]))
  labels <- labels %||% groups
  if (length(x = unique(x = labels.loc[, id])) != length(x = labels)) {
    stop(
      "Length of labels (", length(x = labels), ") must be equal to the number of clusters being labeled (",
      length(x = labels.loc), ")."
    )
  }
  names(x = labels) <- groups
  for (group in groups) {
    labels.loc[labels.loc[, id] == group, id] <- labels[group]
  }
  if (box) {
    geom.use <- ifelse(test = repel, yes = geom_label_repel,
                       no = geom_label
    )
    plot <- plot + geom.use(
      data = labels.loc, mapping = aes_string(
        x = xynames["x"],
        y = xynames["y"], label = id, fill = id
      ), show.legend = FALSE,
      ...
    ) + scale_fill_manual(values = labels.loc$color[order(labels.loc[
      ,
      id,
    ])])
  } else {
    geom.use <- ifelse(test = repel, yes = geom_text_repel, no = geom_text)
    labels.loc <- labels.loc[order(labels.loc$ident),]
    plot <- plot + geom.use(
      data = labels.loc, mapping = aes_string(
        x = xynames["x"],
        y = xynames["y"], label = id
      ), max.iter=max.iter,
      color = color, bg.color = labels.loc$color[order(labels.loc[,id,])], bg.r = bg.r, show.legend = FALSE, ...
    )
  }
  if (geom == "GeomSpatial") {
    plot <- suppressMessages(expr = plot + coord_fixed(
      xlim = xrange.save,
      ylim = yrange.save
    ))
  }
  return(plot)
}





############################################################################








############################################################################



dp <- function(object, dims = c(1, 2), cells = NULL, cols = NULL, pt.s = 0.3, reduction = NULL, gb = NULL, 
               split.by = NULL, shape.by = NULL, order = NULL, shuffle = TRUE, seed = 1, label = TRUE, 
               label.size = 5, label.color = "black", bg.r=0.1, label.box = FALSE, max.iter=1, repel = TRUE, 
               alpha = 0.5, cells.highlight = NULL, cols.highlight = "#DE2D26", sizes.highlight = 1, 
               autosize = FALSE, na.value = "grey50", ncol = NULL, mylegend=F,combine = TRUE, raster = NULL, 
               raster.dpi = c(512, 512), jitter=0) {
  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }
  group.by <- gb
  remove(gb)
  pt.size <- pt.s
  remove(pt.s)
  
  reduction <- reduction %||% DefaultDimReduc(object = object)
  cells <- cells %||% colnames(x = object)
  data <- Embeddings(object = object[[reduction]])[cells, dims]
  data <- as.data.frame(x = data)
  dims <- paste0(Key(object = object[[reduction]]), dims)
  object[["ident"]] <- Idents(object = object)
  orig.groups <- group.by
  group.by <- group.by %||% "ident"
  data <- cbind(data, object[[group.by]][cells, , drop = FALSE])
  group.by <- colnames(x = data)[3:ncol(x = data)]
  for (group in group.by) {
    if (!is.factor(x = data[, group])) {
      data[, group] <- factor(x = data[, group])
    }
  }
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  if (!is.null(x = split.by)) {
    data[, split.by] <- object[[split.by, drop = TRUE]]
  }
  if (isTRUE(x = shuffle)) {
    set.seed(seed = seed)
    data <- data[sample(x = 1:nrow(x = data)), ]
  }
  plots <- lapply(X = group.by, FUN = function(x) {
    x=group.by
    plot <- mySingleDimPlot(
      data = data[, c(
        dims, x, split.by,
        shape.by
      )], dims = dims, col.by = x, cols = cols,
      my.pt.size = pt.size, shape.by = shape.by, order = order,
      label = FALSE, alpha=alpha,cells.highlight = cells.highlight,
      cols.highlight = cols.highlight, sizes.highlight = sizes.highlight, autosize = autosize,
      na.value = na.value, raster = raster, raster.dpi = raster.dpi, jitter=jitter
    )
    if (label) {
      plot <- dplabelClustersWithOutilnedText(plot = plot, id = x, repel = repel, size = label.size, split.by = split.by, box = label.box, color = label.color, alpha=1,max.iter=max.iter, mylegend=mylegend)
    }
    if (!is.null(x = split.by)) {
      plot <- plot + Seurat:::FacetTheme() + facet_wrap(
        facets = vars(!!sym(x = split.by)),
        ncol = if (length(x = group.by) > 1 || is.null(x = ncol)) {
          length(x = unique(x = data[, split.by]))
        } else {
          ncol
        }
      )
    }
    plot <- if (is.null(x = orig.groups)) {
      plot + labs(title = NULL)
    } else {
      plot + Seurat:::CenterTitle()
    }
  })
  if (!is.null(x = split.by)) {
    ncol <- 1
  }
  if (combine) {
    plots <- wrap_plots(plots, ncol = orig.groups %iff% ncol)
  }
  if (mylegend == TRUE){
    return(plots)
  }
  else {
    return(plots + NoLegend())
  }
}





############################################################################










############################################################################



dp2 <- function(object, dims = c(1, 2), cells = NULL, cols = NULL, pt.size = 1, reduction = NULL, group.by = NULL, split.by = NULL, shape.by = NULL, order = NULL, shuffle = FALSE, seed = 1, label = TRUE, label.size = 7, label.color = "black", bg.r=0.1, label.box = FALSE, max.iter=1, repel = TRUE, alpha = 0.5, cells.highlight = NULL, cols.highlight = "#DE2D26", sizes.highlight = 1, autosize = FALSE, na.value = "grey50", ncol = NULL, combine = TRUE, raster = NULL, raster.dpi = c(512, 512)) {
  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }
  reduction <- reduction %||% DefaultDimReduc(object = object)
  cells <- cells %||% colnames(x = object)
  data <- Embeddings(object = object[[reduction]])[cells, dims]
  data <- as.data.frame(x = data)
  dims <- paste0(Key(object = object[[reduction]]), dims)
  object[["ident"]] <- Idents(object = object)
  orig.groups <- group.by
  group.by <- group.by %||% "ident"
  data <- cbind(data, object[[group.by]][cells, , drop = FALSE])
  group.by <- colnames(x = data)[3:ncol(x = data)]
  for (group in group.by) {
    if (!is.factor(x = data[, group])) {
      data[, group] <- factor(x = data[, group])
    }
  }
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  if (!is.null(x = split.by)) {
    data[, split.by] <- object[[split.by, drop = TRUE]]
  }
  if (isTRUE(x = shuffle)) {
    set.seed(seed = seed)
    data <- data[sample(x = 1:nrow(x = data)), ]
  }
  plots <- lapply(X = group.by, FUN = function(x) {
    x="ident"
    plot <- mySingleDimPlot(
      data = data[, c(
        dims, x, split.by,
        shape.by
      )], dims = dims, col.by = x, cols = cols,
      my.pt.size = pt.size, shape.by = shape.by, order = order,
      label = FALSE, alpha=alpha,cells.highlight = cells.highlight,
      cols.highlight = cols.highlight, sizes.highlight = sizes.highlight, autosize = autosize,
      na.value = na.value, raster = raster, raster.dpi = raster.dpi
    )
    if (label) {
      plot <- dplabelClustersWithOutilnedText2(plot = plot, id = x, repel = repel, size = label.size, split.by = split.by, box = label.box, color = label.color, alpha=1,max.iter=max.iter)
    }
    if (!is.null(x = split.by)) {
      plot <- plot + Seurat:::FacetTheme() + facet_wrap(
        facets = vars(!!sym(x = split.by)),
        ncol = if (length(x = group.by) > 1 || is.null(x = ncol)) {
          length(x = unique(x = data[, split.by]))
        } else {
          ncol
        }
      )
    }
    plot <- if (is.null(x = orig.groups)) {
      plot + labs(title = NULL)
    } else {
      plot + Seurat:::CenterTitle()
    }
  })
  if (!is.null(x = split.by)) {
    ncol <- 1
  }
  if (combine) {
    plots <- wrap_plots(plots, ncol = orig.groups %iff% ncol)
  }
  return(plots)
}





############################################################################








############################################################################


mySingleDimPlot <- function(data, dims, col.by = NULL, cols = NULL, pt.size = NULL, my.pt.size = my.pt.size,
                            shape.by = NULL, alpha.by = NULL, order = NULL, label = FALSE,
                            repel = FALSE, label.size = label.size, cells.highlight = NULL, cols.highlight = "#DE2D26", alpha = 1,
                            sizes.highlight = 1, na.value = "grey50", raster = NULL, jitter=jitter,
                            raster.dpi = NULL, autosize=FALSE) {
  if (autosize==TRUE){
    pt.size <- pt.size %||% AutoPointSize(data = data, raster = raster)
  }
  else {
    pt.size <- my.pt.size}
  if ((nrow(x = data) > 1e+05) & !isFALSE(raster)) {
    message(
      "Rasterizing points since number of points exceeds 100,000.",
      "\nTo disable this behavior set `raster=FALSE`"
    )
  }
  raster <- raster %||% (nrow(x = data) > 1e+05)
  if (!is.null(x = raster.dpi)) {
    if (!is.numeric(x = raster.dpi) || length(x = raster.dpi) !=
        2) {
      stop("'raster.dpi' must be a two-length numeric vector")
    }
  }
  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }
  if (!is.data.frame(x = data)) {
    data <- as.data.frame(x = data)
  }
  if (is.character(x = dims) && !all(dims %in% colnames(x = data))) {
    stop("Cannot find dimensions to plot in data")
  } else if (is.numeric(x = dims)) {
    dims <- colnames(x = data)[dims]
  }
  if (!is.null(x = cells.highlight)) {
    highlight.info <- Seurat:::SetHighlight(
      cells.highlight = cells.highlight,
      cells.all = rownames(x = data), sizes.highlight = pt.size,
      pt.size = pt.size, cols.highlight = cols.highlight, col.base = cols[1] %||%
        "#C3C3C3"
    )
    order <- highlight.info$plot.order
    data$highlight <- highlight.info$highlight
    col.by <- "highlight"
    pt.size <- highlight.info$size
    cols <- highlight.info$color
  }
  if (!is.null(x = order) && !is.null(x = col.by)) {
    if (typeof(x = order) == "logical") {
      if (order) {
        data <- data[order(
          !is.na(x = data[, col.by]),
          data[, col.by]
        ), ]
      }
    } else {
      order <- rev(x = c(order, setdiff(x = unique(x = data[
        ,
        col.by
      ]), y = order)))
      data[, col.by] <- factor(x = data[, col.by], levels = order)
      new.order <- order(x = data[, col.by])
      data <- data[new.order, ]
      if (length(x = pt.size) == length(x = new.order)) {
        pt.size <- pt.size[new.order]
      }
    }
  }
  if (!is.null(x = col.by) && !col.by %in% colnames(x = data)) {
    warning("Cannot find ", col.by, " in plotting data, not coloring plot")
    col.by <- NULL
  } else {
    col.index <- match(x = col.by, table = colnames(x = data))
    if (grepl(pattern = "^\\d", x = col.by)) {
      col.by <- paste0("x", col.by)
    } else if (grepl(pattern = "-", x = col.by)) {
      col.by <- gsub(
        pattern = "-", replacement = ".",
        x = col.by
      )
    }
    colnames(x = data)[col.index] <- col.by
  }
  if (!is.null(x = shape.by) && !shape.by %in% colnames(x = data)) {
    warning("Cannot find ", shape.by, " in plotting data, not shaping plot")
  }
  if (!is.null(x = alpha.by) && !alpha.by %in% colnames(x = data)) {
    warning("Cannot find alpha variable ", alpha.by, " in data, setting to NULL",
            call. = FALSE, immediate. = TRUE
    )
    alpha.by <- NULL
  }
  plot <- ggplot(data = data)
  plot <- if (isTRUE(x = raster)) {
    plot + geom_scattermore(mapping = aes_string(
      x = dims[1],
      y = dims[2], color = paste0("`", col.by, "`"), shape = shape.by,
      alpha = alpha.by
    ), pointsize = pt.size, pixels = raster.dpi)
  } else {
    plot + geom_point(mapping = aes_string(
      x = dims[1], y = dims[2],
      color = paste0("`", col.by, "`"), shape = shape.by
      
    ), size = pt.size, alpha = alpha, position=position_jitter(h=jitter,w=jitter)) 
    }
  plot <- plot + guides(color = guide_legend(override.aes = list(size = 3))) +
    labs(color = NULL, title = col.by, ) + Seurat:::CenterTitle()
  if (label && !is.null(x = col.by)) {
    plot <- LabelClusters(
      plot = plot, id = col.by, repel = repel,
      size = label.size
    )
  }
  if (!is.null(x = cols)) {
    if (length(x = cols) == 1 && (is.numeric(x = cols) ||
                                  cols %in% rownames(x = brewer.pal.info))) {
      scale <- scale_color_brewer(palette = cols, na.value = na.value)
    } else if (length(x = cols) == 1 && (cols %in% c(
      "alphabet",
      "alphabet2", "glasbey", "polychrome", "stepped"
    ))) {
      colors <- DiscretePalette(length(unique(data[[col.by]])),
                                palette = cols
      )
      scale <- scale_color_manual(values = colors, na.value = na.value)
    } else {
      scale <- scale_color_manual(values = cols, na.value = na.value)
    }
    plot <- plot + scale
  }
  plot <- plot + theme_cowplot()
  return(plot)
}




###############################################################################


netVisual_heatmap <- function (object, comparison = c(1, 2), measure = c("count",
                                                                         "weight"), signaling = NULL, slot.name = c("netP", "net"),
                               color.use = NULL, color.heatmap = c("#2166ac", "#b2182b"),
                               title.name = NULL, width = NULL, height = NULL, font.size = 8,
                               font.size.title = 10, cluster.rows = FALSE, cluster.cols = FALSE,
                               sources.use = NULL, targets.use = NULL, remove.isolate = FALSE,
                               row.show = NULL, col.show = NULL, col.rot = 45)
{
  if (!is.null(measure)) {
    measure <- match.arg(measure)
  }
  slot.name <- match.arg(slot.name)
  if (is.list(object@net[[1]])) {
    message("Do heatmap based on a merged object \n")
    obj1 <- object@net[[comparison[1]]][[measure]]
    obj2 <- object@net[[comparison[2]]][[measure]]
    net.diff <- obj2 - obj1
    if (measure == "count") {
      if (is.null(title.name)) {
        title.name = "Differential number of interactions"
      }
    }
    else if (measure == "weight") {
      if (is.null(title.name)) {
        title.name = "Differential interaction strength"
      }
    }
    legend.name = "Relative values"
  }
  else {
    message("Do heatmap based on a single object \n")
    if (!is.null(signaling)) {
      net.diff <- slot(object, slot.name)$prob[, , signaling]
      if (is.null(title.name)) {
        title.name = paste0(signaling, " signaling network")
      }
      legend.name <- "Communication Prob."
    }
    else if (!is.null(measure)) {
      net.diff <- object@net[[measure]]
      if (measure == "count") {
        if (is.null(title.name)) {
          title.name = "Number of interactions"
        }
      }
      else if (measure == "weight") {
        if (is.null(title.name)) {
          title.name = "Interaction strength"
        }
      }
      legend.name <- title.name
    }
  }
  net <- net.diff
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source", "target")
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) {
        sources.use <- rownames(net.diff)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) {
        targets.use <- rownames(net.diff)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- rownames(net.diff)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]],
                                          df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }
  mat <- net
  if (is.null(color.use)) {
    color.use <- scPalette(ncol(mat))
  }
  names(color.use) <- colnames(mat)
  if (!is.null(row.show)) {
    mat <- mat[row.show, ]
  }
  if (!is.null(col.show)) {
    mat <- mat[, col.show]
    color.use <- color.use[col.show]
  }
  if (min(mat) < 0) {
    color.heatmap.use = colorRamp3(c(min(mat), 0, max(mat)),
                                   c(color.heatmap[1], "#f7f7f7", color.heatmap[2]))
    colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*",
                                                                      "\\1", min(mat, na.rm = T))) + 1), 0, round(max(mat,
                                                                                                                      na.rm = T), digits = nchar(sub(".*\\.(0*).*", "\\1",
                                                                                                                                                     max(mat, na.rm = T))) + 1))
  }
  else {
    if (length(color.heatmap) == 3) {
      color.heatmap.use = colorRamp3(c(0, min(mat), max(mat)),
                                     color.heatmap)
    }
    else if (length(color.heatmap) == 2) {
      color.heatmap.use = colorRamp3(c(min(mat), max(mat)),
                                     color.heatmap)
    }
    else if (length(color.heatmap) == 1) {
      color.heatmap.use = (grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9,
                                                                                 name = color.heatmap))))(100)
    }
    colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*",
                                                                      "\\1", min(mat, na.rm = T))) + 1), round(max(mat,
                                                                                                                   na.rm = T), digits = nchar(sub(".*\\.(0*).*", "\\1",
                                                                                                                                                  max(mat, na.rm = T))) + 1))
  }
  df <- data.frame(group = colnames(mat))
  rownames(df) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use),
                                      which = "column", show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  row_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use),
                                      which = "row", show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  ha1 = rowAnnotation(Strength = anno_barplot(rowSums(abs(mat)),
                                              border = FALSE, gp = gpar(fill = color.use, col = color.use)),
                      show_annotation_name = FALSE)
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(abs(mat)),
                                                  border = FALSE, gp = gpar(fill = color.use, col = color.use)),
                          show_annotation_name = FALSE)
  if (sum(abs(mat) > 0) == 1) {
    color.heatmap.use = c("white", color.heatmap.use)
  }
  else {
    mat[mat == 0] <- NA
  }
  ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white",
                name = legend.name, bottom_annotation = col_annotation,
                left_annotation = row_annotation, top_annotation = ha2,
                right_annotation = ha1, cluster_rows = cluster.rows,
                cluster_columns = cluster.rows, row_names_side = "left",
                row_names_rot = 0, row_names_gp = gpar(fontsize = font.size),
                column_names_gp = gpar(fontsize = font.size), column_title = title.name,
                column_title_gp = gpar(fontsize = font.size.title), column_names_rot = col.rot,
                row_title = "Sources (Sender)", row_title_gp = gpar(fontsize = font.size.title),
                row_title_rot = 90, heatmap_legend_param = list(title_gp = gpar(fontsize = 8,
                                                                                fontface = "plain"), title_position = "leftcenter-rot",
                                                                border = NA, legend_height = unit(20, "mm"), labels_gp = gpar(fontsize = 8),
                                                                grid_width = unit(2, "mm")))
  return(ht1)
}


####################################################





####################################################

netVisual_heatmap <- function (object, comparison = c(1, 2), measure = c("count",
                                                                         "weight"), signaling = NULL, slot.name = c("netP", "net"),
                               color.use = NULL, color.heatmap = c("#2166ac", "#b2182b"),
                               title.name = NULL, width = NULL, height = NULL, font.size = 8,
                               font.size.title = 10, cluster.rows = FALSE, cluster.cols = FALSE,
                               sources.use = NULL, targets.use = NULL, remove.isolate = FALSE,
                               row.show = NULL, col.show = NULL)
{
  if (!is.null(measure)) {
    measure <- match.arg(measure)
  }
  slot.name <- match.arg(slot.name)
  if (is.list(object@net[[1]])) {
    message("Do heatmap based on a merged object \n")
    obj1 <- object@net[[comparison[1]]][[measure]]
    obj2 <- object@net[[comparison[2]]][[measure]]
    net.diff <- obj2 - obj1
    if (measure == "count") {
      if (is.null(title.name)) {
        title.name = "Differential number of interactions"
      }
    }
    else if (measure == "weight") {
      if (is.null(title.name)) {
        title.name = "Differential interaction strength"
      }
    }
    legend.name = "Relative values"
  }
  else {
    message("Do heatmap based on a single object \n")
    if (!is.null(signaling)) {
      net.diff <- slot(object, slot.name)$prob[, , signaling]
      if (is.null(title.name)) {
        title.name = paste0(signaling, " signaling network")
      }
      legend.name <- "Communication Prob."
    }
    else if (!is.null(measure)) {
      net.diff <- object@net[[measure]]
      if (measure == "count") {
        if (is.null(title.name)) {
          title.name = "Number of interactions"
        }
      }
      else if (measure == "weight") {
        if (is.null(title.name)) {
          title.name = "Interaction strength"
        }
      }
      legend.name <- title.name
    }
  }
  net <- net.diff
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source", "target")
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) {
        sources.use <- rownames(net.diff)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) {
        targets.use <- rownames(net.diff)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- rownames(net.diff)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]],
                                          df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }
  mat <- net
  if (is.null(color.use)) {
    color.use <- scPalette(ncol(mat))
  }
  names(color.use) <- colnames(mat)
  if (!is.null(row.show)) {
    mat <- mat[row.show, ]
  }
  if (!is.null(col.show)) {
    mat <- mat[, col.show]
    color.use <- color.use[col.show]
  }
  # if (min(mat) < 0) {
  #      color.heatmap.use = colorRamp3(c(min(mat), 0, max(mat)),
  #                                     c(color.heatmap[1], "#f7f7f7", color.heatmap[2]))
  #      colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*",
  #                                                                        "\\1", min(mat, na.rm = T))) + 1), 0, round(max(mat,
  #                                                                                                                        na.rm = T), digits = nchar(sub(".*\\.(0*).*", "\\1",
  #                                                                                                                                                       max(mat, na.rm = T))) + 1))
  # }
  # else {
  #      if (length(color.heatmap) == 3) {
  #           color.heatmap.use = colorRamp3(c(0, min(mat), max(mat)),
  #                                          color.heatmap)
  #      }
  #      else if (length(color.heatmap) == 2) {
  #           color.heatmap.use = colorRamp3(c(min(mat), max(mat)),
  #                                          color.heatmap)
  #      }
  #      else if (length(color.heatmap) == 1) {
  #           color.heatmap.use = (grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9,
  #                                                                                      name = color.heatmap))))(100)
  #      }
  #      colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*",
  #                                                                        "\\1", min(mat, na.rm = T))) + 1), round(max(mat,
  #                                                                                                                     na.rm = T), digits = nchar(sub(".*\\.(0*).*", "\\1",
  #                                                                                                                                                    max(mat, na.rm = T))) + 1))
  # }
  color.heatmap.fun <- colorRampPalette(c("white","#b2182b"))
  color.heatmap.use <- color.heatmap.fun(100)
  df <- data.frame(group = colnames(mat))
  rownames(df) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use),
                                      which = "column", show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  row_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use),
                                      which = "row", show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  ha1 = rowAnnotation(Strength = anno_barplot(rowSums(abs(mat)),
                                              border = FALSE, gp = gpar(fill = color.use, col = color.use)),
                      show_annotation_name = FALSE)
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(abs(mat)),
                                                  border = FALSE, gp = gpar(fill = color.use, col = color.use)),
                          show_annotation_name = FALSE)
  if (sum(abs(mat) > 0) == 1) {
    color.heatmap.use = c("white", color.heatmap.use)
  }
  # else {
  #      mat[mat == 0] <- NA
  # }
  #color.heatmap.use <- c("white",color.heatmap.use)
  ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white",
                name = legend.name, bottom_annotation = col_annotation,
                left_annotation = row_annotation, top_annotation = ha2,
                right_annotation = ha1, cluster_rows = cluster.rows,
                cluster_columns = cluster.cols, row_names_side = "left",
                row_names_rot = 0, row_names_gp = gpar(fontsize = font.size),
                column_names_gp = gpar(fontsize = font.size), column_title = title.name,
                column_title_gp = gpar(fontsize = font.size.title), column_names_rot = 45,
                row_title = "Sources (Sender)", row_title_gp = gpar(fontsize = font.size.title),
                row_title_rot = 90, heatmap_legend_param = list(title_gp = gpar(fontsize = 8,
                                                                                fontface = "plain"), title_position = "leftcenter-rot",
                                                                border = NA, legend_height = unit(20, "mm"), labels_gp = gpar(fontsize = 8),
                                                                grid_width = unit(2, "mm")))
  return(ht1)
}







####################################################








#####################################################




SpatialPlot <- function (object, group.by = NULL, features = NULL, images = NULL,
                         cols = NULL, image.alpha = 1, crop = TRUE, slot = "data",
                         min.cutoff = NA, max.cutoff = NA, cells.highlight = NULL,
                         cols.highlight = c("#DE2D26", "grey50"), facet.highlight = FALSE,
                         label = FALSE, label.size = 5, label.color = "white", label.box = TRUE,
                         repel = FALSE, ncol = NULL, combine = TRUE, pt.size.factor = 1.6,
                         alpha = c(1, 1), stroke = 0.25, interactive = FALSE, do.identify = FALSE,
                         identify.ident = NULL, do.hover = FALSE, information = NULL)
{
  if (isTRUE(x = do.hover) || isTRUE(x = do.identify)) {
    warning("'do.hover' and 'do.identify' are deprecated as we are removing plotly-based interactive graphics, use 'interactive' instead for Shiny-based interactivity",
            call. = FALSE, immediate. = TRUE)
    interactive <- TRUE
  }
  if (!is.null(x = group.by) & !is.null(x = features)) {
    stop("Please specific either group.by or features, not both.")
  }
  images <- images %||% Images(object = object, assay = DefaultAssay(object = object))
  if (length(x = images) == 0) {
    images <- Images(object = object)
  }
  if (length(x = images) < 1) {
    stop("Could not find any spatial image information")
  }
  if (is.null(x = features)) {
    if (interactive) {
      return(ISpatialDimPlot(object = object, image = images[1],
                             group.by = group.by, alpha = alpha))
    }
    group.by <- group.by %||% "ident"
    object[["ident"]] <- Idents(object = object)
    data <- object[[group.by]]
    for (group in group.by) {
      if (!is.factor(x = data[, group])) {
        data[, group] <- factor(x = data[, group])
      }
    }
  }
  else {
    if (interactive) {
      return(ISpatialFeaturePlot(object = object, feature = features[1],
                                 image = images[1], slot = slot, alpha = alpha))
    }
    data <- FetchData(object = object, vars = features, slot = slot)
    features <- colnames(x = data)
    min.cutoff <- mapply(FUN = function(cutoff, feature) {
      return(ifelse(test = is.na(x = cutoff), yes = min(data[,
                                                             feature]), no = cutoff))
    }, cutoff = min.cutoff, feature = features)
    max.cutoff <- mapply(FUN = function(cutoff, feature) {
      return(ifelse(test = is.na(x = cutoff), yes = max(data[,
                                                             feature]), no = cutoff))
    }, cutoff = max.cutoff, feature = features)
    check.lengths <- unique(x = vapply(X = list(features,
                                                min.cutoff, max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
    if (length(x = check.lengths) != 1) {
      stop("There must be the same number of minimum and maximum cuttoffs as there are features")
    }
    data <- sapply(X = 1:ncol(x = data), FUN = function(index) {
      data.feature <- as.vector(x = data[, index])
      min.use <- Seurat:::SetQuantile(cutoff = min.cutoff[index],
                                      data.feature)
      max.use <- Seurat:::SetQuantile(cutoff = max.cutoff[index],
                                      data.feature)
      data.feature[data.feature < min.use] <- min.use
      data.feature[data.feature > max.use] <- max.use
      return(data.feature)
    })
    colnames(x = data) <- features
    rownames(x = data) <- Cells(x = object)
  }
  features <- colnames(x = data)
  colnames(x = data) <- features
  rownames(x = data) <- colnames(x = object)
  facet.highlight <- facet.highlight && (!is.null(x = cells.highlight) &&
                                           is.list(x = cells.highlight))
  if (do.hover) {
    if (length(x = images) > 1) {
      images <- images[1]
      warning("'do.hover' requires only one image, using image ",
              images, call. = FALSE, immediate. = TRUE)
    }
    if (length(x = features) > 1) {
      features <- features[1]
      type <- ifelse(test = is.null(x = group.by), yes = "feature",
                     no = "grouping")
      warning("'do.hover' requires only one ", type, ", using ",
              features, call. = FALSE, immediate. = TRUE)
    }
    if (facet.highlight) {
      warning("'do.hover' requires no faceting highlighted cells",
              call. = FALSE, immediate. = TRUE)
      facet.highlight <- FALSE
    }
  }
  if (facet.highlight) {
    if (length(x = images) > 1) {
      images <- images[1]
      warning("Faceting the highlight only works with a single image, using image ",
              images, call. = FALSE, immediate. = TRUE)
    }
    ncols <- length(x = cells.highlight)
  }
  else {
    ncols <- length(x = images)
  }
  plots <- vector(mode = "list", length = length(x = features) *
                    ncols)
  for (i in 1:ncols) {
    plot.idx <- i
    image.idx <- ifelse(test = facet.highlight, yes = 1,
                        no = i)
    image.use <- object[[images[[image.idx]]]]
    coordinates <- GetTissueCoordinates(object = image.use)
    highlight.use <- if (facet.highlight) {
      cells.highlight[i]
    }
    else {
      cells.highlight
    }
    for (j in 1:length(x = features)) {
      cols.unset <- is.factor(x = data[, features[j]]) &&
        is.null(x = cols)
      if (cols.unset) {
        cols <- hue_pal()(n = length(x = levels(x = data[,
                                                         features[j]])))
        names(x = cols) <- levels(x = data[, features[j]])
      }
      plot <- SingleSpatialPlot(data = cbind(coordinates,
                                             data[rownames(x = coordinates), features[j],
                                                  drop = FALSE]), image = image.use, image.alpha = image.alpha,
                                col.by = features[j], cols = cols, alpha.by = if (is.null(x = group.by)) {
                                  features[j]
                                }
                                else {
                                  NULL
                                }, pt.alpha = if (!is.null(x = group.by)) {
                                  alpha[j]
                                }
                                else {
                                  NULL
                                }, geom = if (inherits(x = image.use, what = "STARmap")) {
                                  "poly"
                                }
                                else {
                                  "spatial"
                                }, cells.highlight = highlight.use, cols.highlight = cols.highlight,
                                pt.size.factor = pt.size.factor, stroke = stroke,
                                crop = crop)
      if (is.null(x = group.by)) {
        plot <- plot + scale_fill_gradientn(name = features[j],
                                            colours = SpatialColors(n = 100)) + theme(legend.position = "top") +
          scale_alpha(range = alpha) + guides(alpha = FALSE)
      }
      else if (label) {
        plot <- LabelClusters(plot = plot, id = ifelse(test = is.null(x = cells.highlight),
                                                       yes = features[j], no = "highlight"), geom = if (inherits(x = image.use,
                                                                                                                 what = "STARmap")) {
                                                         "GeomPolygon"
                                                       }
                              else {
                                "GeomSpatial"
                              }, repel = repel, size = label.size, color = label.color,
                              box = label.box, position = "nearest")
      }
      if (j == 1 && length(x = images) > 1 && !facet.highlight) {
        plot <- plot + ggtitle(label = images[[image.idx]]) +
          theme(plot.title = element_text(hjust = 0.5))
      }
      if (facet.highlight) {
        plot <- plot + ggtitle(label = names(x = cells.highlight)[i]) +
          theme(plot.title = element_text(hjust = 0.5)) +
          NoLegend()
      }
      plots[[plot.idx]] <- plot
      plot.idx <- plot.idx + ncols
      if (cols.unset) {
        cols <- NULL
      }
    }
  }
  if (length(x = images) > 1 && combine) {
    plots <- wrap_plots(plots = plots, ncol = length(x = images))
  }
  else if (length(x = images == 1) && combine) {
    plots <- wrap_plots(plots = plots, ncol = ncol)
  }
  return(plots)
}


######################################################################







######################################################################

netVisual_heatmap <- function (object, comparison = c(1, 2), measure = c("count",
                                                                         "weight"), signaling = NULL, slot.name = c("netP", "net"),
                               color.use = NULL, color.heatmap = c("#2166ac", "#b2182b"),
                               title.name = NULL, width = NULL, height = NULL, font.size = 8,
                               font.size.title = 10, cluster.rows = FALSE, cluster.cols = FALSE,
                               sources.use = NULL, targets.use = NULL, remove.isolate = FALSE,
                               row.show = NULL, col.show = NULL)
{
  if (!is.null(measure)) {
    measure <- match.arg(measure)
  }
  slot.name <- match.arg(slot.name)
  if (is.list(object@net[[1]])) {
    message("Do heatmap based on a merged object \n")
    obj1 <- object@net[[comparison[1]]][[measure]]
    obj2 <- object@net[[comparison[2]]][[measure]]
    net.diff <- obj2 - obj1
    if (measure == "count") {
      if (is.null(title.name)) {
        title.name = "Differential number of interactions"
      }
    }
    else if (measure == "weight") {
      if (is.null(title.name)) {
        title.name = "Differential interaction strength"
      }
    }
    legend.name = "Relative values"
  }
  else {
    message("Do heatmap based on a single object \n")
    if (!is.null(signaling)) {
      net.diff <- slot(object, slot.name)$prob[, , signaling]
      if (is.null(title.name)) {
        title.name = paste0(signaling, " signaling network")
      }
      legend.name <- "Communication Prob."
    }
    else if (!is.null(measure)) {
      net.diff <- object@net[[measure]]
      if (measure == "count") {
        if (is.null(title.name)) {
          title.name = "Number of interactions"
        }
      }
      else if (measure == "weight") {
        if (is.null(title.name)) {
          title.name = "Interaction strength"
        }
      }
      legend.name <- title.name
    }
  }
  net <- net.diff
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source", "target")
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) {
        sources.use <- rownames(net.diff)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) {
        targets.use <- rownames(net.diff)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- rownames(net.diff)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]],
                                          df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }
  mat <- net
  if (is.null(color.use)) {
    color.use <- scPalette(ncol(mat))
  }
  names(color.use) <- colnames(mat)
  if (!is.null(row.show)) {
    mat <- mat[row.show, ]
  }
  if (!is.null(col.show)) {
    mat <- mat[, col.show]
    color.use <- color.use[col.show]
  }
  # if (min(mat) < 0) {
  #      color.heatmap.use = colorRamp3(c(min(mat), 0, max(mat)),
  #                                     c(color.heatmap[1], "#f7f7f7", color.heatmap[2]))
  #      colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*",
  #                                                                        "\\1", min(mat, na.rm = T))) + 1), 0, round(max(mat,
  #                                                                                                                        na.rm = T), digits = nchar(sub(".*\\.(0*).*", "\\1",
  #                                                                                                                                                       max(mat, na.rm = T))) + 1))
  # }
  # else {
  #      if (length(color.heatmap) == 3) {
  #           color.heatmap.use = colorRamp3(c(0, min(mat), max(mat)),
  #                                          color.heatmap)
  #      }
  #      else if (length(color.heatmap) == 2) {
  #           color.heatmap.use = colorRamp3(c(min(mat), max(mat)),
  #                                          color.heatmap)
  #      }
  #      else if (length(color.heatmap) == 1) {
  #           color.heatmap.use = (grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9,
  #                                                                                      name = color.heatmap))))(100)
  #      }
  #      colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*",
  #                                                                        "\\1", min(mat, na.rm = T))) + 1), round(max(mat,
  #                                                                                                                     na.rm = T), digits = nchar(sub(".*\\.(0*).*", "\\1",
  #                                                                                                                                                    max(mat, na.rm = T))) + 1))
  # }
  color.heatmap.fun <- colorRampPalette(c("royalblue4","#b2182b"))#("dodgerblue4","dodgerblue2","deepskyblue1","white","#b2182b")
  color.heatmap.use <- color.heatmap.fun(100)
  df <- data.frame(group = colnames(mat))
  rownames(df) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use),
                                      which = "column", show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  row_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use),
                                      which = "row", show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  ha1 = rowAnnotation(Strength = anno_barplot(rowSums(abs(mat)),
                                              border = FALSE, gp = gpar(fill = color.use, col = color.use)),
                      show_annotation_name = FALSE)
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(abs(mat)),
                                                  border = FALSE, gp = gpar(fill = color.use, col = color.use)),
                          show_annotation_name = FALSE)
  if (sum(abs(mat) > 0) == 1) {
    color.heatmap.use = c("white", color.heatmap.use)
  }
  # else {
  #      mat[mat == 0] <- NA
  # }
  #color.heatmap.use <- c("white",color.heatmap.use)
  ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white",
                name = legend.name, bottom_annotation = col_annotation,
                left_annotation = row_annotation, top_annotation = ha2,
                right_annotation = ha1, cluster_rows = cluster.rows,
                cluster_columns = cluster.cols, row_names_side = "left",
                row_names_rot = 0, row_names_gp = gpar(fontsize = font.size),
                column_names_gp = gpar(fontsize = font.size), column_title = title.name,
                column_title_gp = gpar(fontsize = font.size.title), column_names_rot = 45,
                row_title = "Sources (Sender)", row_title_gp = gpar(fontsize = font.size.title),
                row_title_rot = 90, heatmap_legend_param = list(title_gp = gpar(fontsize = 8,
                                                                                fontface = "plain"), title_position = "leftcenter-rot",
                                                                border = TRUE, border_gp = 0.2, legend_height = unit(20, "mm"), labels_gp = gpar(fontsize = 8),
                                                                grid_width = unit(2, "mm")),
                column_split=1:dim(mat)[[2]], row_split=1:dim(mat)[[1]],
                row_gap=unit(0.5, "mm"), column_gap=unit(0.5, "mm"))
  return(ht1)
}



######################################################################





######################################################################

computeEnrichmentScore <- function(df, measure = c("ligand", "signaling","LR-pair"), species = c('mouse','human'), color.use = NULL, color.name = "Dark2", n.color = 8,
                                   scale=c(4,.8), min.freq = 0, max.words = 200, random.order = FALSE, rot.per = 0,return.data = FALSE,seed = 1,...) {
  measure <- match.arg(measure)
  species <- match.arg(species)
  LRpairs <- as.character(unique(df$interaction_name))
  ES <- vector(length = length(LRpairs))
  for (i in 1:length(LRpairs)) {
    df.i <- subset(df, interaction_name == LRpairs[i])
    if (length(which(rowSums(is.na(df.i)) > 0)) > 0) {
      df.i <- df.i[-which(rowSums(is.na(df.i)) > 0), ,drop = FALSE]
    }
    ES[i] = mean(abs(df.i$ligand.logFC) * abs(df.i$receptor.logFC) *abs(df.i$ligand.pct.2-df.i$ligand.pct.1)*abs(df.i$receptor.pct.2-df.i$receptor.pct.1))
  }
  if (species == "mouse") {
    CellChatDB <- CellChatDB.mouse
  } else if (species == 'human') {
    CellChatDB <- CellChatDB.human
  }
  df.es <- CellChatDB$interaction[LRpairs, c("ligand",'receptor','pathway_name')]
  df.es$score <- ES
  # summarize the enrichment score
  df.es.ensemble <- df.es %>% group_by(ligand) %>% summarize(total = sum(score))  # avg = mean(score),
  
  set.seed(seed)
  if (is.null(color.use)) {
    color.use <- RColorBrewer::brewer.pal(n.color, color.name)
  }
  
  wordcloud::wordcloud(words = df.es.ensemble$ligand, freq = df.es.ensemble$total, min.freq = min.freq, max.words = max.words,scale=scale,
                       random.order = random.order, rot.per = rot.per, colors = color.use,...)
  if (return.data) {
    return(df.es.ensemble)
  }
}


############################################################






############################################################

plotGeneExpression <- function (object, features = NULL, signaling = NULL, enriched.only = TRUE,
                                type = c("violin", "dot", "bar"), color.use = NULL, group.by = NULL, text.size = 15,
                                ...)
{
  type <- match.arg(type)
  meta <- object@meta
  if (is.list(object@idents)) {
    meta$group.cellchat <- object@idents$joint
  }
  else {
    meta$group.cellchat <- object@idents
  }
  if (!identical(rownames(meta), colnames(object@data.signaling))) {
    cat("The cell barcodes in 'meta' is ", head(rownames(meta)),
        "\n")
    warning("The cell barcodes in 'meta' is different from those in the used data matrix.\n              We now simply assign the colnames in the data matrix to the rownames of 'mata'!")
    rownames(meta) <- colnames(object@data.signaling)
  }
  w10x <- Seurat::CreateSeuratObject(counts = object@data.signaling,
                                     meta.data = meta)
  if (is.null(group.by)) {
    group.by <- "group.cellchat"
  }
  Seurat::Idents(w10x) <- group.by
  if (!is.null(features) & !is.null(signaling)) {
    warning("`features` will be used when inputing both `features` and `signaling`!")
  }
  if (!is.null(features)) {
    feature.use <- features
  }
  else if (!is.null(signaling)) {
    res <- CellChat:::extractEnrichedLR(object, signaling = signaling,
                                        geneLR.return = TRUE, enriched.only = enriched.only)
    feature.use <- res$geneLR
  }
  if (type == "violin") {
    gg <- myStackedVlnPlot(w10x, features = feature.use, color.use = color.use, text.size = text.size,
                           ...)
  }
  else if (type == "dot") {
    gg <- dotPlot(w10x, features = feature.use, color.use = color.use,
                  ...)
  }
  else if (type == "bar") {
    gg <- barPlot(w10x, features = feature.use, color.use = color.use,
                  ...)
  }
  return(gg)
}


############################################################







############################################################

myStackedVlnPlot <- function (object, features, idents = NULL, split.by = NULL, color.use = NULL,
                              colors.ggplot = FALSE, show.median = FALSE, median.size = 1,
                              angle.x = 45, vjust.x = NULL, hjust.x = NULL, show.text.y = TRUE, text.size = 15,
                              line.size = NULL, pt.size = 0, plot.margin = margin(0, 0,
                                                                                  0, 0, "cm"), ...)
{
  options(warn = -1)
  if (is.null(color.use)) {
    numCluster <- length(levels(Seurat::Idents(object)))
    if (colors.ggplot) {
      color.use <- NULL
    }
    else {
      color.use <- scPalette(numCluster)
    }
  }
  if (is.null(vjust.x) | is.null(hjust.x)) {
    angle = c(0, 45, 90)
    hjust = c(0, 1, 1)
    vjust = c(0, 1, 0.5)
    vjust.x = vjust[angle == angle.x]
    hjust.x = hjust[angle == angle.x]
  }
  plot_list <- purrr::map(features, function(x) mymodify_vlnplot(object = object,
                                                                 features = x, idents = idents, split.by = split.by, cols = color.use,
                                                                 show.median = show.median, median.size = median.size,
                                                                 pt.size = pt.size, show.text.y = show.text.y, line.size = line.size, text.size = text.size,
                                                                 ...))
  plot_list[[length(plot_list)]] <- plot_list[[length(plot_list)]] +
    theme(axis.text.x = element_text(), axis.ticks.x = element_line()) +
    theme(axis.text.x = element_text(angle = angle.x, hjust = hjust.x,
                                     vjust = vjust.x)) +
    theme(axis.text.x = element_text(size = text.size), axis.text.y = element_text(size = 10))
  ymaxs <- purrr::map_dbl(plot_list, CellChat:::extract_max)
  plot_list <- purrr::map2(plot_list, ymaxs, function(x, y) x +
                             scale_y_continuous(breaks = c(y)) + expand_limits(y = y))
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1) +
    patchwork::plot_layout(guides = "collect")
  return(p)
}




############################################################







############################################################


mymodify_vlnplot <- function (object, features, idents = NULL, split.by = NULL, cols = NULL,
                              show.median = FALSE, median.size = 1, show.text.y = TRUE,
                              line.size = NULL, pt.size = 0, plot.margin = margin(0, 0,
                                                                                  0, 0, "cm"), text.size = 15, ...)
{
  options(warn = -1)
  p <- Seurat::VlnPlot(object, features = features, cols = cols,
                       pt.size = pt.size, idents = idents, split.by = split.by,
                       ...) + xlab("") + ylab(features) + ggtitle("")
  if (show.median) {
    p <- p + stat_summary(fun.y = median, geom = "point",
                          shape = 3, size = median.size)
  }
  p <- p + theme(text = element_text(size = text.size)) + theme(axis.line = element_line(size = line.size)) +
    theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
          axis.line.x = element_line(colour = "black", size = line.size),
          axis.line.y = element_line(colour = "black", size = line.size))
  p <- p + theme(plot.title = element_blank(), axis.title.x = element_blank(),
                 axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                 axis.title.y = element_text(size = rel(1), angle = 0),
                 axis.text.y = element_text(size = rel(1)), plot.margin = plot.margin) +
    theme(axis.text.y = element_text(size = 8))
  p <- p + theme(element_line(size = line.size))
  if (!show.text.y) {
    p <- p + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
  }
  return(p)
}







####################################################################




# seuObj -- spatial seurat object, myexpression -- vector or list in format c("Itga4 > 0", "Ccr2 < 1", ....)
dps <- function(seuObj = seuObj, myexpression = list()) {
  myexpression <- strsplit(myexpression, split = " ")
  
  mycells <- list()
  n <- 0
  for (i in myexpression) {
    i <- unlist(i)
    n <- n + 1
    if (i[[2]] == ">") {
      mycells[[n]] <- colnames(seuObj@assays$SCT@data)[seuObj@assays$SCT@data[i[[1]], ] > i[[3]]]
    } else if (i[[2]] == "<") {
      mycells[[n]] <- colnames(seuObj@assays$SCT@data)[seuObj@assays$SCT@data[i[[1]], ] < i[[3]]]
    } else if (i[[2]] == "==") {
      mycells[[n]] <- colnames(seuObj@assays$SCT@data)[seuObj@assays$SCT@data[i[[1]], ] == i[[3]]]
    }
  }
  
  intercells <- mycells[[1]]
  n <- 1
  for (i in 2:(length(mycells))) {
    n <- n + 1
    intercells <- intersect(intercells, mycells[[n]])
  }
  print(paste("ncell =",length(intercells)),sep='')
  #return(SpatialDimPlot(seuObj, cells.highlight = intercells))
  return(intercells)
}




###################################################################


myDoHeatmap <- function (object, features = NULL, cells = NULL, group.by = "ident",
                         group.bar = TRUE, group.colors = NULL, disp.min = -2.5, disp.max = NULL,
                         slot = "scale.data", assay = NULL, label = TRUE, size = 5.5,
                         hjust = 0, angle = 45, raster = TRUE, draw.lines = TRUE,
                         lines.width = NULL, group.bar.height = 0.02, combine = TRUE, myPallete = c( "slategray1", 
                                                                                                     "grey88",
                                                                                                    "coral",
                                                                                                    "firebrick"), ...)
{
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  features <- rev(x = unique(x = features))
  disp.max <- disp.max %||% ifelse(test = slot == "scale.data",
                                   yes = 2.5, no = 6)
  possible.features <- rownames(x = GetAssayData(object = object,
                                                 slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if (length(x = features) == 0) {
      stop("No requested features found in the ", slot,
           " slot for the ", assay, " assay.")
    }
    warning("The following features were omitted as they were not found in the ",
            slot, " slot for the ", assay, " assay: ", paste(bad.features,
                                                             collapse = ", "))
  }
  data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(object = object,
                                                             slot = slot)[features, cells, drop = FALSE])))
  object <- suppressMessages(expr = StashIdent(object = object,
                                               save.name = "ident"))
  group.by <- group.by %||% "ident"
  groups.use <- object[[group.by]][cells, , drop = FALSE]
  plots <- vector(mode = "list", length = ncol(x = groups.use))
  for (i in 1:ncol(x = groups.use)) {
    data.group <- data
    group.use <- groups.use[, i, drop = TRUE]
    if (!is.factor(x = group.use)) {
      group.use <- factor(x = group.use)
    }
    names(x = group.use) <- cells
    if (draw.lines) {
      lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) *
                                                0.0025)
      placeholder.cells <- sapply(X = 1:(length(x = levels(x = group.use)) *
                                           lines.width), FUN = function(x) {
                                             return(RandomName(length = 20))
                                           })
      placeholder.groups <- rep(x = levels(x = group.use),
                                times = lines.width)
      group.levels <- levels(x = group.use)
      names(x = placeholder.groups) <- placeholder.cells
      group.use <- as.vector(x = group.use)
      names(x = group.use) <- cells
      group.use <- factor(x = c(group.use, placeholder.groups),
                          levels = group.levels)
      na.data.group <- matrix(data = NA, nrow = length(x = placeholder.cells),
                              ncol = ncol(x = data.group), dimnames = list(placeholder.cells,
                                                                           colnames(x = data.group)))
      data.group <- rbind(data.group, na.data.group)
    }
    lgroup <- length(levels(group.use))
    plot <- SingleRasterMap(data = data.group, raster = raster,
                            disp.min = disp.min, disp.max = disp.max, feature.order = features,
                            cell.order = names(x = sort(x = group.use)), group.by = group.use, colors=myPallete)
    if (group.bar) {
      default.colors <- c(hue_pal()(length(x = levels(x = group.use))))
      if (!is.null(x = names(x = group.colors))) {
        cols <- unname(obj = group.colors[levels(x = group.use)])
      }
      else {
        cols <- group.colors[1:length(x = levels(x = group.use))] %||%
          default.colors
      }
      if (any(is.na(x = cols))) {
        cols[is.na(x = cols)] <- default.colors[is.na(x = cols)]
        cols <- col2hex(cols)
        col.dups <- sort(x = unique(x = which(x = duplicated(x = substr(x = cols,
                                                                        start = 1, stop = 7)))))
        through <- length(x = default.colors)
        while (length(x = col.dups) > 0) {
          pal.max <- length(x = col.dups) + through
          cols.extra <- hue_pal()(pal.max)[(through +
                                              1):pal.max]
          cols[col.dups] <- cols.extra
          col.dups <- sort(x = unique(x = which(x = duplicated(x = substr(x = cols,
                                                                          start = 1, stop = 7)))))
        }
      }
      group.use2 <- sort(x = group.use)
      if (draw.lines) {
        na.group <- RandomName(length = 20)
        levels(x = group.use2) <- c(levels(x = group.use2),
                                    na.group)
        group.use2[placeholder.cells] <- na.group
        cols <- c(cols, "#FFFFFF")
      }
      pbuild <- ggplot_build(plot = plot)
      names(x = cols) <- levels(x = group.use2)
      y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
      y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) +
        y.range * 0.015
      y.max <- y.pos + group.bar.height * y.range
      x.min <- min(pbuild$layout$panel_params[[1]]$x.range) +
        0.1
      x.max <- max(pbuild$layout$panel_params[[1]]$x.range) -
        0.1
      plot <- plot + annotation_raster(raster = t(x = cols[group.use2]),
                                       xmin = x.min, xmax = x.max, ymin = y.pos, ymax = y.max) +
        coord_cartesian(ylim = c(0, y.max), clip = "off") +
        scale_color_manual(values = cols[-length(x = cols)],
                           name = "Identity", na.translate = FALSE)
      if (label) {
        x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
        x.divs <- pbuild$layout$panel_params[[1]]$x.major %||%
          attr(x = pbuild$layout$panel_params[[1]]$x$get_breaks(),
               which = "pos")
        x <- data.frame(group = sort(x = group.use),
                        x = x.divs)
        label.x.pos <- tapply(X = x$x, INDEX = x$group,
                              FUN = function(y) {
                                if (isTRUE(x = draw.lines)) {
                                  mean(x = y[-length(x = y)])
                                }
                                else {
                                  mean(x = y)
                                }
                              })
        label.x.pos <- data.frame(group = names(x = label.x.pos),
                                  label.x.pos)
        plot <- plot + geom_text(stat = "identity", data = label.x.pos,
                                 aes_string(label = "group", x = "label.x.pos"),
                                 y = y.max + y.max * 0.03 * 0.5, angle = angle,
                                 hjust = hjust, size = size)
        plot <- suppressMessages(plot + coord_cartesian(ylim = c(0,
                                                                 y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use))) *
                                                                   size), clip = "off"))
      }
    }
    plot <- plot + theme(line = element_blank())
    plots[[i]] <- plot
  }
  if (combine) {
    plots <- wrap_plots(plots)
  }
  return(plots)
}






##########################################################



netVisual_heatmap <- function (object, comparison = c(1, 2), measure = c("count", 
                                                                         "weight"), signaling = NULL, slot.name = c("netP", "net"), 
                               color.use = NULL, color.heatmap = c("#2166ac", "#b2182b"), 
                               title.name = NULL, width = NULL, height = NULL, font.size = 8, 
                               font.size.title = 10, cluster.rows = FALSE, cluster.cols = FALSE, 
                               sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, 
                               row.show = NULL, col.show = NULL) 
{
  if (!is.null(measure)) {
    measure <- match.arg(measure)
  }
  slot.name <- match.arg(slot.name)
  if (is.list(object@net[[1]])) {
    message("Do heatmap based on a merged object \n")
    obj1 <- object@net[[comparison[1]]][[measure]]
    obj2 <- object@net[[comparison[2]]][[measure]]
    net.diff <- obj2 - obj1
    if (measure == "count") {
      if (is.null(title.name)) {
        title.name = "Differential number of interactions"
      }
    }
    else if (measure == "weight") {
      if (is.null(title.name)) {
        title.name = "Differential interaction strength"
      }
    }
    legend.name = "Relative values"
  }
  else {
    message("Do heatmap based on a single object \n")
    if (!is.null(signaling)) {
      net.diff <- slot(object, slot.name)$prob[, , signaling]
      if (is.null(title.name)) {
        title.name = paste0(signaling, " signaling network")
      }
      legend.name <- "Communication Prob."
    }
    else if (!is.null(measure)) {
      net.diff <- object@net[[measure]]
      if (measure == "count") {
        if (is.null(title.name)) {
          title.name = "Number of interactions"
        }
      }
      else if (measure == "weight") {
        if (is.null(title.name)) {
          title.name = "Interaction strength"
        }
      }
      legend.name <- title.name
    }
  }
  net <- net.diff
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source", "target")
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) {
        sources.use <- rownames(net.diff)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) {
        targets.use <- rownames(net.diff)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- rownames(net.diff)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], 
                                          df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
    if (length(idx) > 0) {
      net <- net[-idx, ]
      net <- net[, -idx]
    }
  }
  mat <- net
  if (is.null(color.use)) {
    color.use <- scPalette(ncol(mat))
  }
  names(color.use) <- colnames(mat)
  if (!is.null(row.show)) {
    mat <- mat[row.show, ]
  }
  if (!is.null(col.show)) {
    mat <- mat[, col.show]
    color.use <- color.use[col.show]
  }
  if (min(mat) < 0) {
    color.heatmap.use = colorRamp3(c(min(mat), 0, max(mat)), 
                                   c(color.heatmap[1], "#f7f7f7", color.heatmap[2]))
    colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*", 
                                                                      "\\1", min(mat, na.rm = T))) + 1), 0, round(max(mat, 
                                                                                                                      na.rm = T), digits = nchar(sub(".*\\.(0*).*", "\\1", 
                                                                                                                                                     max(mat, na.rm = T))) + 1))
  }
  else {
    if (length(color.heatmap) == 3) {
      color.heatmap.use = colorRamp3(c(0, min(mat), max(mat)), 
                                     color.heatmap)
    }
    else if (length(color.heatmap) == 2) {
      color.heatmap.use = colorRamp3(c(min(mat), max(mat)), 
                                     color.heatmap)
    }
    else if (length(color.heatmap) == 1) {
      color.heatmap.use = (grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, 
                                                                                 name = color.heatmap))))(100)
    }
    colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*", 
                                                                      "\\1", min(mat, na.rm = T))) + 1), round(max(mat, 
                                                                                                                   na.rm = T), digits = nchar(sub(".*\\.(0*).*", "\\1", 
                                                                                                                                                  max(mat, na.rm = T))) + 1))
  }
  df <- data.frame(group = colnames(mat))
  rownames(df) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), 
                                      which = "column", show_legend = FALSE, show_annotation_name = FALSE, 
                                      simple_anno_size = grid::unit(0.2, "cm"))
  row_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), 
                                      which = "row", show_legend = FALSE, show_annotation_name = FALSE, 
                                      simple_anno_size = grid::unit(0.2, "cm"))
  ha1 = rowAnnotation(Strength = anno_barplot(rowSums(abs(mat)), 
                                              border = FALSE, gp = gpar(fill = color.use, col = color.use)), 
                      show_annotation_name = FALSE)
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(abs(mat)), 
                                                  border = FALSE, gp = gpar(fill = color.use, col = color.use)), 
                          show_annotation_name = FALSE)
  if (sum(abs(mat) > 0) == 1) {
    color.heatmap.use = c("white", color.heatmap.use)
  }
  else {
    mat[mat == 0] <- NA
  }
  library(circlize)
  ht1 = Heatmap(mat, col = colorRamp2(c(1,0,-1), c("#b2182b","white","#2166ac")), na_col = "white", 
                name = legend.name, bottom_annotation = col_annotation, 
                left_annotation = row_annotation, top_annotation = ha2, 
                right_annotation = ha1, cluster_rows = cluster.rows, 
                cluster_columns = cluster.rows, row_names_side = "left", 
                row_names_rot = 0, row_names_gp = gpar(fontsize = font.size), 
                column_names_gp = gpar(fontsize = font.size), column_title = title.name, 
                column_title_gp = gpar(fontsize = font.size.title), column_names_rot = 45, 
                row_title = "Sources (Sender)", row_title_gp = gpar(fontsize = font.size.title), 
                row_title_rot = 90, heatmap_legend_param = list(title_gp = gpar(fontsize = 8, 
                                                                                fontface = "plain"), title_position = "leftcenter-rot", 
                                                                border = NA, legend_height = unit(20, "mm"), labels_gp = gpar(fontsize = 8), 
                                                                grid_width = unit(2, "mm")))
  return(ht1)
}







#####################################################



mymarkers <- function(seuobject=NULL, clusterOfInterest=NULL, versus=NULL){
  library(rio)
  markers <- FindMarkers(seuobject, clusterOfInterest, versus, min.pct=0.2, logfc.threshold=0.5)
  markers <- markers[match(sort(markers$avg_log2FC,decreasing=T),markers$avg_log2FC),]
  markers <- markers[which(markers$p_val_adj < 0.05),]
  markers <- markers[c(1:25,(length(markers[[1]])-25):(length(markers[[1]]))),]
  write.table(markers,"markers.csv",sep=',')
  convert("markers.csv","markers.xlsx")
  return(markers)
}



#####################################################




FindAllMarkersBulk <- function (seurat=mk, clus_ident='lab', sample_ident='sampletype', expfilt_counts = 10, 
                                expfilt_freq = 0.5, n_top_genes = 50, pct.in = 20, out_dir = "FindMarkersBulk_outs", 
                                alpha = 0.1, assay = "SCT", heatmap=FALSE) {
  
  start <- Sys.time()
  coef <- variable <- value <- NULL
  dir.create(out_dir, showWarnings = FALSE)
  Idents(seurat) <- clus_ident
  clusters <- unique(Idents(seurat))
  clusters <- sort(clusters)
  pdf(paste(out_dir, "/cells_per_clus_HM.pdf", sep = ""))
  pheatmap(table(seurat@meta.data[, clus_ident], seurat@meta.data[, 
                                                                  sample_ident]), display_numbers = T, cluster_rows = F, 
           cluster_cols = F, fontsize_number = 4)
  dev.off()
  wilcox <- wilcoxauc(seurat, assay = "data", seurat_assay = assay, 
                      group_by = clus_ident)
  top_markers <- vector()
  mylist <- list()
  n=0
  for (cluster in clusters) {
    n=n+1
    groups <- seurat@meta.data[, c(clus_ident, sample_ident)]
    groups$iscluster <- as.vector(groups[[clus_ident]])
    groups$iscluster[which(groups$iscluster != cluster)] <- "other"
    pb <- Matrix.utils:::aggregate.Matrix(t(seurat@assays[[assay]]@counts), 
                           groupings = groups[, 2:3], fun = "sum")
    splitf <- sapply(stringr::str_split(rownames(pb), pattern = "_", 
                                        n = 2), `[`, 2)
    splitf[which(splitf != cluster)] <- "other"
    cluster_counts <- as.data.frame(t(as.matrix(pb)))
    cluster_metadata <- data.frame(sample_ident = sapply(stringr::str_split(colnames(cluster_counts), 
                                                                            pattern = "_", n = 2), `[`, 1), iscluster = sapply(stringr::str_split(colnames(cluster_counts), 
                                                                                                                                                  pattern = "_", n = 2), `[`, 2))
    cluster_metadata$iscluster <- factor(cluster_metadata$iscluster, 
                                         levels = c("other", as.vector(cluster)))
    dds <- DESeqDataSetFromMatrix(cluster_counts, colData = cluster_metadata, 
                                  design = ~iscluster)
    keep <- rowSums(counts(dds) >= expfilt_counts) >= expfilt_freq * 
      nrow(cluster_metadata)
    dds <- dds[keep, ]
    vst <- varianceStabilizingTransformation(dds)
    dds <- DESeq(dds, test = "LRT", reduced = ~1)
    res <- results(dds, alpha = 0.05)
    res_shrink <- lfcShrink(dds, coef = colnames(coef(dds))[2], 
                            res = res, type = "apeglm")
    res_shrink <- as.data.frame(res_shrink)
    res_shrink$sig <- rep("Not significant", nrow(res_shrink))
    res_shrink$sig[which(res_shrink$padj < alpha)] <- "Significant"
    d <- plotCounts(dds, gene = which.min(res$padj), intgroup = "iscluster", 
                    returnData = TRUE)
    exp_gene_logfc <- res[which.min(res$padj), ]$log2FoldChange
    exp_gene_padj <- res[which.min(res$padj), ]$padj
    gg_counts <- cluster_counts[, sort(colnames(cluster_counts))]
    gg_counts <- melt(log10(gg_counts))
    pdf(paste(out_dir, "/cluster_", cluster, "_diagnostic_plots.pdf", 
              sep = ""))
    print(ggplot(gg_counts, aes(x = variable, y = value, 
                                fill = variable)) + geom_boxplot() + theme_bw() + 
            theme(axis.text.x = element_text(angle = 90), legend.position = "none") + 
            ylab("Log10(Counts)"))
    print(DESeq2::plotPCA(vst, intgroup = "iscluster") + 
            theme_classic())
    print(DESeq2::plotPCA(vst, intgroup = "sample_ident") + 
            theme_classic() + geom_text_repel(aes(label = sample_ident), 
                                              show.legend = FALSE))
    plotDispEsts(dds)
    plotMA(res)
    
    print(ggplot(res_shrink, aes(x = log2FoldChange, y = -log10(pvalue), 
                                 color = sig)) + geom_point(alpha = 0.7) + 
            scale_color_manual(values = c("grey40", 
                                          "blue")) + theme_classic())
    print(ggplot(d, aes(x = iscluster, y = count)) + geom_point(position = position_jitter(w = 0.1, 
                                                                                           h = 0)) + ggtitle(paste("Gene:", row.names(res)[which.min(res$padj)], 
                                                                                                                   "\nLog2FC = ", exp_gene_logfc, "\npadj = ", exp_gene_padj, 
                                                                                                                   sep = "")))
    dev.off()
    wilcox_sub <- wilcox[which(wilcox$group == cluster), 
    ]
    res_shrink$feature <- row.names(res_shrink)
    merged_res <- merge(x = res_shrink, y = wilcox_sub, by = "feature", 
                        sort = F)
    res_shrink$pct_in <- merged_res$pct_in
    res_shrink$pct_out <- merged_res$pct_out
    res_shrink$feature <- NULL
    write.csv(file = paste(out_dir, "/cluster_", cluster, 
                           "_results.csv", sep = ""), res_shrink)
    res_sig <- res_shrink[which(res_shrink$padj < 0.05 & 
                                  res_shrink$pct_in > pct.in), ]
    top_markers <- c(top_markers, row.names(res_sig[order(res_sig$log2FoldChange, 
                                                          decreasing = T), ])[1:n_top_genes])
    mylist[[n]] <- res_shrink
  }
  write.csv(top_markers, file = paste(out_dir, "/Top_markers.csv", 
                                      sep = ""), row.names = F, quote = F)
  
  if (heatmap==TRUE){
    pdf(file = paste(out_dir, "/Top_markers_HM.pdf", sep = ""))
    print(DoHeatmap(subset(seurat, downsample = 1000), features = top_markers$x, 
                    assay = 'SCT', slot = "scale.data", raster = F) + 
            scale_fill_gradient2(low = rev(c("#d1e5f0",  "#67a9cf", "#2166ac")), mid = "white", high = rev(c("#b2182b", 
                                                                                                             "#ef8a62", "#fddbc7")), midpoint = 0, guide = "colourbar", 
                                 aesthetics = "fill", na.value = "white") + theme(text = element_text(size = 5)))
    dev.off()
  }
  print(start)
  print(Sys.time())
  return(mylist)
}

















