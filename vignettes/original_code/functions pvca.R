plot.pvca<-function (pvcaObj, cex.percentage = 1, fname = NULL, ht = 4, 
          wd = 5, title = fname) 
{
  require(stringr)
  require(ggplot2)
  title <- mgsub(c("/", ".", "PVCA"), rep("", 3), fixed = TRUE, 
                 fname)
  labels <- gsub(":", " x ", pvcaObj$label)
  db <- data.frame(perc = t(pvcaObj$dat), labels = factor(labels, 
                                                          levels = labels))
  p <- ggplot(db, aes(x = labels, y = perc)) + 
    geom_bar(stat = "identity", fill = "blue") + theme_bw() + 
    scale_y_continuous(limits = c(0,  1.1)) + 
    geom_text(aes(x = labels, y = perc + 0.05, label = paste0(as.character(round(100 *  perc, 1)), "%"))) +
    labs(x = "Effects", y = "Weighted average proportion variance", title = paste("PVCA ", title)) + 
    scale_x_discrete(labels = function(x) str_wrap(x, width = 8))
  if (!is.null(fname)) {
    ggsave(file = paste0(fname, ".pdf"), height = ht, width = wd, 
           plot = p)
  }
  return(p)
}

getEigen<-function(eset){
  theDataMatrix <- exprs(eset)
  dataRowN <- nrow(theDataMatrix)
  dataColN <- ncol(theDataMatrix)
  
  theDataMatrixCentered <- matrix(data = 0, nrow = dataRowN, 
                                  ncol = dataColN)
  theDataMatrixCentered_transposed <- apply(theDataMatrix, 
                                            1, scale, center = TRUE, scale = FALSE)
  theDataMatrixCentered <- t(theDataMatrixCentered_transposed)
  theDataCor <- cor(theDataMatrixCentered)
  eigenData <- eigen(theDataCor)
  
  return(eigenData)
}


pvcaBatchAssess.MSF2<-function(eset, eigenData=NULL, batch.factors, threshold, include.inter="all") {
  require(lme4)
  theDataMatrix <- exprs(eset)
  dataRowN <- nrow(theDataMatrix)
  dataColN <- ncol(theDataMatrix)
  
  if (is.null(eigenData)){
    theDataMatrixCentered <- matrix(data = 0, nrow = dataRowN, 
                                    ncol = dataColN)
    theDataMatrixCentered_transposed <- apply(theDataMatrix, 
                                              1, scale, center = TRUE, scale = FALSE)
    theDataMatrixCentered <- t(theDataMatrixCentered_transposed)
    theDataCor <- cor(theDataMatrixCentered)
    eigenData <- eigen(theDataCor)
  }
  
  
  eigenValues <- eigenData$values
  ev_n <- length(eigenValues)
  eigenVectorsMatrix <- eigenData$vectors
  eigenValuesSum <- sum(eigenValues)
  percents_PCs <- eigenValues/eigenValuesSum
  expInfo <- pData(eset)[, batch.factors]
  exp_design <- as.data.frame(expInfo)
  expDesignRowN <- nrow(exp_design)
  expDesignColN <- ncol(exp_design)
  my_counter_2 <- 0
  my_sum_2 <- 1
  for (i in ev_n:1) {
    my_sum_2 = my_sum_2 - percents_PCs[i]
    if ((my_sum_2) <= threshold) {
      my_counter_2 = my_counter_2 + 1
    }
  }
  if (my_counter_2 < 3) {
    pc_n = 3
  }
  else {
    pc_n = my_counter_2
  }
  pc_data_matrix <- matrix(data = 0, nrow = (expDesignRowN * 
                                               pc_n), ncol = 1)
  mycounter = 0
  for (i in 1:pc_n) {
    for (j in 1:expDesignRowN) {
      mycounter <- mycounter + 1
      pc_data_matrix[mycounter, 1] = eigenVectorsMatrix[j, 
                                                        i]
    }
  }
  AAA <- exp_design[rep(1:expDesignRowN, pc_n), ]
  Data <- cbind(AAA, pc_data_matrix)
  variables <- c(colnames(exp_design))
  for (i in 1:length(variables)) {
    Data$variables[i] <- as.factor(Data$variables[i])
  }
  op <- options(warn = (-1))
  model.func <- c()
  index <- 1
  for (i in 1:length(variables)) {
    mod = paste("(1|", variables[i], ")", sep = "")
    model.func[index] = mod
    index = index + 1
  }
  for (i in 1:(length(variables) - 1)) {
    for (j in (i + 1):length(variables)) {
      mod = paste("(1|", variables[i], ":", variables[j], 
                    ")", sep = "")
      model.func[index] = mod
      index = index + 1
    }
  }
  
  if (include.inter!="all"){
    i.delete.RE <- setdiff(grep(":", model.func), grep(include.inter, 
                                                       model.func))
    delete.RE <- model.func[i.delete.RE]
    model.func <- setdiff(model.func, delete.RE)
  }
  
  effects_n = length(model.func) + 1
  randomEffectsMatrix <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
  function.mods <- paste(model.func, collapse = " + ")
  for (i in 1:pc_n) {
    y = (((i - 1) * expDesignRowN) + 1)
    funct <- paste("pc_data_matrix", function.mods, sep = " ~ ")
    Rm1ML <- lmer(funct, Data[y:(((i - 1) * expDesignRowN) + 
                                   expDesignRowN), ], REML = TRUE, verbose = FALSE, 
                  na.action = na.omit)
    randomEffects <- Rm1ML
    randomEffectsMatrix[i, ] <- c(unlist(VarCorr(Rm1ML)), 
                                  resid = sigma(Rm1ML)^2)
  }
  effectsNames <- c(names(getME(Rm1ML, "cnms")), "resid")
  randomEffectsMatrixStdze <- matrix(data = 0, nrow = pc_n, 
                                     ncol = effects_n)
  for (i in 1:pc_n) {
    mySum = sum(randomEffectsMatrix[i, ])
    for (j in 1:effects_n) {
      randomEffectsMatrixStdze[i, j] = randomEffectsMatrix[i, 
                                                           j]/mySum
    }
  }
  randomEffectsMatrixWtProp <- matrix(data = 0, nrow = pc_n, 
                                      ncol = effects_n)
  for (i in 1:pc_n) {
    weight = eigenValues[i]/eigenValuesSum
    for (j in 1:effects_n) {
      randomEffectsMatrixWtProp[i, j] = randomEffectsMatrixStdze[i, 
                                                                 j] * weight
    }
  }
  randomEffectsSums <- matrix(data = 0, nrow = 1, ncol = effects_n)
  randomEffectsSums <- colSums(randomEffectsMatrixWtProp)
  totalSum <- sum(randomEffectsSums)
  randomEffectsMatrixWtAveProp <- matrix(data = 0, nrow = 1, 
                                         ncol = effects_n)
  for (j in 1:effects_n) {
    randomEffectsMatrixWtAveProp[j] = randomEffectsSums[j]/totalSum
  }
  return(list(dat = randomEffectsMatrixWtAveProp, label = effectsNames))
}

pvcaBatchAssess.MSF<-function (eset, batch.factors, threshold, include.inter = NULL) 
{
  require(lme4)
  theDataMatrix <- exprs(eset)
  dataRowN <- nrow(theDataMatrix)
  dataColN <- ncol(theDataMatrix)
  theDataMatrixCentered <- matrix(data = 0, nrow = dataRowN, 
                                  ncol = dataColN)
  theDataMatrixCentered_transposed <- apply(theDataMatrix, 
                                            1, scale, center = TRUE, scale = FALSE)
  theDataMatrixCentered <- t(theDataMatrixCentered_transposed)
  theDataCor <- cor(theDataMatrixCentered)
  eigenData <- eigen(theDataCor)
  eigenValues <- eigenData$values
  ev_n <- length(eigenValues)
  eigenVectorsMatrix <- eigenData$vectors
  eigenValuesSum <- sum(eigenValues)
  percents_PCs <- eigenValues/eigenValuesSum
  expInfo <- pData(eset)[, batch.factors]
  exp_design <- as.data.frame(expInfo)
  expDesignRowN <- nrow(exp_design)
  expDesignColN <- ncol(exp_design)
  my_counter_2 <- 0
  my_sum_2 <- 1
  for (i in ev_n:1) {
    my_sum_2 = my_sum_2 - percents_PCs[i]
    if ((my_sum_2) <= threshold) {
      my_counter_2 = my_counter_2 + 1
    }
  }
  if (my_counter_2 < 3) {
    pc_n = 3
  }
  else {
    pc_n = my_counter_2
  }
  pc_data_matrix <- matrix(data = 0, nrow = (expDesignRowN * 
                                               pc_n), ncol = 1)
  mycounter = 0
  for (i in 1:pc_n) {
    for (j in 1:expDesignRowN) {
      mycounter <- mycounter + 1
      pc_data_matrix[mycounter, 1] = eigenVectorsMatrix[j, 
                                                        i]
    }
  }
  AAA <- exp_design[rep(1:expDesignRowN, pc_n), ]
  Data <- cbind(AAA, pc_data_matrix)
  variables <- c(colnames(exp_design))
  for (i in 1:length(variables)) {
    Data$variables[i] <- as.factor(Data$variables[i])
  }
  op <- options(warn = (-1))
  model.func <- c()
  index <- 1
  for (i in 1:length(variables)) {
    mod = paste("(1|", variables[i], ")", sep = "")
    model.func[index] = mod
    index = index + 1
  }
  for (i in 1:(length(variables) - 1)) {
    for (j in (i + 1):length(variables)) {
      mod = paste("(1|", variables[i], ":", variables[j], 
                    ")", sep = "")
      model.func[index] = mod
      index = index + 1
    }
  }
  i.delete.RE <- setdiff(grep(":", model.func), grep(include.inter, 
                                                     model.func))
  delete.RE <- model.func[i.delete.RE]
  model.func <- setdiff(model.func, delete.RE)
  effects_n = length(model.func) + 1
  randomEffectsMatrix <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
  function.mods <- paste(model.func, collapse = " + ")
  for (i in 1:pc_n) {
    y = (((i - 1) * expDesignRowN) + 1)
    funct <- paste("pc_data_matrix", function.mods, sep = " ~ ")
    Rm1ML <- lmer(funct, Data[y:(((i - 1) * expDesignRowN) + 
                                   expDesignRowN), ], REML = TRUE, verbose = FALSE, 
                  na.action = na.omit)
    randomEffects <- Rm1ML
    randomEffectsMatrix[i, ] <- c(unlist(VarCorr(Rm1ML)), 
                                  resid = sigma(Rm1ML)^2)
  }
  effectsNames <- c(names(getME(Rm1ML, "cnms")), "resid")
  randomEffectsMatrixStdze <- matrix(data = 0, nrow = pc_n, 
                                     ncol = effects_n)
  for (i in 1:pc_n) {
    mySum = sum(randomEffectsMatrix[i, ])
    for (j in 1:effects_n) {
      randomEffectsMatrixStdze[i, j] = randomEffectsMatrix[i, 
                                                           j]/mySum
    }
  }
  randomEffectsMatrixWtProp <- matrix(data = 0, nrow = pc_n, 
                                      ncol = effects_n)
  for (i in 1:pc_n) {
    weight = eigenValues[i]/eigenValuesSum
    for (j in 1:effects_n) {
      randomEffectsMatrixWtProp[i, j] = randomEffectsMatrixStdze[i, 
                                                                 j] * weight
    }
  }
  randomEffectsSums <- matrix(data = 0, nrow = 1, ncol = effects_n)
  randomEffectsSums <- colSums(randomEffectsMatrixWtProp)
  totalSum <- sum(randomEffectsSums)
  randomEffectsMatrixWtAveProp <- matrix(data = 0, nrow = 1, 
                                         ncol = effects_n)
  for (j in 1:effects_n) {
    randomEffectsMatrixWtAveProp[j] = randomEffectsSums[j]/totalSum
  }
  return(list(dat = randomEffectsMatrixWtAveProp, label = effectsNames))
}