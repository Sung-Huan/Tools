singleChoiceWithSubParamValue <- function(paramFile, groupInd, subParamName, subParamTag) {
    suppressWarnings(fileLines <- readLines(paramFile))
    subParamInd <- 0
    itemInd <- 0
    detectValue <- FALSE
    detectItem <- FALSE
    for (j in 1:length(fileLines)) {
        if (detectItem) {
            if (grepl('<Item>', fileLines[j])) {
                if (itemInd == valueInd) {
                    return(strsplit(strsplit(fileLines[j], '>')[[1]][[2]], '<')[[1]][[1]])
                    detectValue <- FALSE
                    detectItem <- FALSE
                }
                itemInd <- itemInd + 1
            }
        }
        if (detectValue) {
            if (grepl('<Value>', fileLines[j])) {
                valueInd <- as.numeric(strsplit(strsplit(fileLines[j], '>')[[1]][[2]], '<')[[1]][[1]])
                if (subParamTag == 'Value') {
                    return(valueInd)
                }
                detectValue <- FALSE
                detectItem <- TRUE
                itemInd <- 0
            }
        }
        if (grepl(subParamName, fileLines[j])) {
            if (subParamInd == groupInd) {
                detectValue <- TRUE
            }
            subParamInd <- subParamInd + 1
        }

    }
}

computeLog <- function(dfLog, resP, resAdjP, i) {
    if (!is.na(resP)) {
        if (resP == 0)
            logP <- -1
        else
            logP <- -log10(resP)
    } else {
        logP <- NaN
    }
    if (!is.na(resAdjP)) {
        if (resAdjP == 0)
            logadjP <- -1
        else
            logadjP <- -log10(resAdjP)
    } else {
        logadjP <- NaN
    }
    dfLog$'logP'[i] <- logP
    dfLog$'logadjP'[i] <- logadjP
    return(dfLog)
}

checkSignificant <- function(parameters, result, lfcs, program, fcInd, pInd, adjpInd) {
    adjP <- intParamValue(parameters, "Max. Adjusted p-value (FDR)")
    setAdjP <- boolParamValue(parameters, "Adjusted p-value (FDR)")
    P <- intParamValue(parameters, "Max. p-value")
    setP <- boolParamValue(parameters, "P-value")
    lfcMethod <- singleChoiceParamValue(parameters, "Log2 Fold Change")
    NameAdjP <- paste('-log(', colnames(result)[adjpInd], ')', sep = '')
    if (program == 'SAM')
        NameP <- '-log(p.values)'
    else
        NameP <- paste('-log(', colnames(result)[pInd], ')', sep = '')
    dfLog <- data.frame('logP' = numeric(nrow(result)), 'logadjP' = numeric(nrow(result)))
    significant <- c()
    for (i in 1:nrow(result)) {
        resAdjP <- result[i, adjpInd]
        resP <- result[i, pInd]
        logFC <- result[i, fcInd]
#        if (is.na(resAdjP) || is.na(resP)) {
#            errorStr <- paste("The matrix contains some values which are relatively small.",
#                              "Their p values become NA. Please normalize or filter the data.", sep = " ")
#            stop(errorStr, call.=FALSE)
#        }
        dfLog <- computeLog(dfLog, resP, resAdjP, i)
        if ((!is.na(resAdjP)) && (setAdjP)) {
            if (resAdjP <= adjP) {
                if (logFC > 0) {
                    significant[i] <- "+"
                } else if (logFC < 0) {
                    significant[i] <- "+"
                } else {
                    significant[i] <- ""
                }
            } else {
                significant[i] <- ""
            }
        } else {
            significant[i] <- ""
        }
        if ((!is.na(logFC)) && (lfcMethod != 'None')) {
            if ((logFC > 0) && (logFC < lfcs[['up']])) {
                significant[i] <- ""
            } else if ((logFC < 0) && (logFC > lfcs[['down']])) {
                significant[i] <- ""
            }
            if (!setAdjP) {
                if (logFC >= lfcs[['up']]) {
                    significant[i] <- "+"
                } else if (logFC <= lfcs[['down']]) {
                    significant[i] <- "+"
                }
            }
        } else if (!setAdjP) {
            significant[i] <- ""
        }
        if ((!is.na(resP)) && (setP)) {
            if (resP > P) {
                significant[i] <- ""
            }
            if ((!setAdjP) && (lfcMethod == 'None')) {
                if (resP <= P) {
                    if (logFC > 0) {
                        significant[i] <- "+"
                    } else if (logFC < 0) {
                        significant[i] <- "+"
                    } else {
                        significant[i] <- ""
                    }
                }
            }
        } else if ((!setAdjP) && (lfcMethod == 'None')) {
            significant[i] <- ""
        }
    }
    replaceP <- max(dfLog$'logP')+1
    replaceadjP <- max(dfLog$'logadjP')+1
    for (i in 1:nrow(dfLog)) {
        if (!is.na(dfLog$'logP'[i])) {
            if (dfLog$'logP'[i] == -1)
                dfLog$'logP'[i] <- replaceP
        }
        if (!is.na(dfLog$'logadjP'[i])) {
            if (dfLog$'logadjP'[i] == -1)
                dfLog$'logadjP'[i] <- replaceadjP
        }
    }
    names(dfLog)[1] <- NameP
    names(dfLog)[2] <- NameAdjP
    return(list(significant, dfLog))
}

getCutoff <- function(parameters, result, program) {
    lfcs <- list('up' = 0, 'down' = 0)
    lfcMethod <- singleChoiceParamValue(parameters, "Log2 Fold Change")
    lfcInd <- singleChoiceParamInd(parameters, "Log2 Fold Change")
    up <- singleChoiceWithSubParamValue(paramFile, lfcInd, "Up-regluation", "Value")
    down <- singleChoiceWithSubParamValue(paramFile, lfcInd, "Down-regluation", "Value")
    if (lfcMethod == 'Number') {
        lfcs[['up']] <- up
        lfcs[['down']] <- down
    } else if (lfcMethod == 'Percentage') {
        if ((program == 'EdgeR') || (program == 'Limma and DEqMS')) {
            lfcList <- result[, 1]
        } else if ((program == 'DESeq2') || (program == 'ROTS')) {
            lfcList <- result[, 2]
        } else if (program == 'SAM') {
            lfcList <- result[, 4]
        }
        lfcs[['up']] <- max(abs(lfcList)) * (up / 100)
        lfcs[['down']] <- max(abs(lfcList)) * (down / 100)
    }
    return(lfcs)
}

generateOutput <- function(parameters, sort_result, program, pair1, pair2, Valid, testMethod) {
    lfcs <- getCutoff(parameters, sort_result, program)
    if (((program == 'EdgeR') && (testMethod != 'Exact Test')) || (program == 'Limma and DEqMS')) {
        statData <- checkSignificant(parameters, sort_result, lfcs, program, 1, 4, 5)
    } else if (program == 'DESeq2') {
        statData <- checkSignificant(parameters, sort_result, lfcs, program, 2, 5, 6)
    } else if ((program == 'EdgeR') && (testMethod == 'Exact Test')) {
        statData <- checkSignificant(parameters, sort_result, lfcs, program, 1, 3, 4)
    } else if (program == 'SAM') {
        statData <- checkSignificant(parameters, sort_result, lfcs, program, 1, 5, 5)
    } else if (program == 'ROTS') {
        statData <- checkSignificant(parameters, sort_result, lfcs, program, 2, 3, 4)
    }
    Significant <- statData[1]
    names(Significant)[1] <- "Significant"
    dfLog <- statData[2]
    df <- as.data.frame(Significant)
    if (program != 'EdgeR') {
        dfV <- as.data.frame(Valid)
        combine <- cbind(sort_result, dfV)
        combine <- cbind(combine, df)
    } else {
        combine <- cbind(sort_result, df)
    }
    combine <- cbind(combine, dfLog)
    for (i in 1:length(names(combine))) {
        colnames(combine)[i] <- paste(pair1, 'vs', pair2, colnames(combine)[i], sep = "_")
    }
    if (program == 'SAM') {
        combine<-combine[-c(ncol(combine)-1)]
    }
    return(combine)
}

checkInteger <- function(data) {
    for (i in 1:nrow(data)) {
        for (j in 1:ncol(data)) {
            if ((data[i, j] < 0) || (is.na(data[i, j])) || ((data[i, j] %% 1) > 0)) {
                stopStr <- paste("The table contains not integer values ",
                                 "(decimal value or NA) which are not acceptable ",
                                 "for DESeq2. Please use 'DE analysis->Adjust Values' for ",
                                 "modifying your table.", sep = '')
                stop(stopStr, call. = FALSE)
            }
        }
    }
}

runEdgeRWithExactTest <- function(parameters, groupInd, testMethod, mdata, program, pair1, pair2, normMethod) {
    log2scale <- boolParamValue(parameters, "Data in Log2 scale")
    library('limma')
    library('edgeR')
    data_raw <- main(mdata)
    if (log2scale) {
        data_raw <- '^'(2,data_raw)
    }
    mobDataGroups <- c()
    for (i in 1:length(annotRows(mdata)[, groupInd + 1])) {
        mobDataGroups[i] <- annotRows(mdata)[, groupInd + 1][[i]]
    }
    d <- DGEList(counts = data_raw, group = factor(mobDataGroups))
    d <- calcNormFactors(d, method = normMethod)
    d <- estimateDisp(d)
    de <- exactTest(d)
    result <- topTags(de, n = Inf)
    result_table <- as.data.frame(result)
    sort_result <- result_table[order(as.numeric(row.names(result_table))),]
    combine <- generateOutput(parameters, sort_result, program, pair1, pair2, c(), testMethod)
    return(combine)
}

runEdgeRWithGLM <- function(parameters, groupInd, testMethod, mdata, program, pair1, pair2, samples, normMethod) {
    log2scale <- boolParamValue(parameters, "Data in Log2 scale")
    library('limma')
    library('edgeR')
    data_raw <- main(mdata)
    if (log2scale) {
        data_raw <- '^'(2,data_raw)
    }
    mobDataGroups <- c()
    for (i in 1:length(samples)) {
        mobDataGroups[i] <- samples[[i]]
    }
    d <- DGEList(counts = data_raw, group = factor(mobDataGroups))
    d <- calcNormFactors(d, method = normMethod)
    design.mat <- model.matrix(~0 + d$samples$group)
    colnames(design.mat) <- levels(d$samples$group)
    contrastPair <- c()
    for (i in 1:length(colnames(design.mat))) {
        if (colnames(design.mat)[[i]] == pair1) {
            contrastPair[i] <- -1
        } else if (colnames(design.mat)[[i]] == pair2) {
            contrastPair[i] <- 1
        } else {
            contrastPair[i] <- 0
        }
    }
    d2 <- estimateDisp(d, design.mat)
    if (testMethod == "Likelihood ratio test") {
        fit <- glmFit(d2, design.mat)
        lrt <- glmLRT(fit, contrast = contrastPair)
    } else if (testMethod == "Quasi-likelihood F-test") {
        fit <- glmQLFit(d2, design.mat)
        lrt <- glmQLFTest(fit, contrast = contrastPair)
    }
    result <- topTags(lrt, n = Inf)
    result_table <- as.data.frame(result)
    sort_result <- result_table[order(as.numeric(row.names(result_table))),]
    combine <- generateOutput(parameters, sort_result, program, pair1, pair2, c(), testMethod)
    return(combine)
}

checkValid <- function(res) {
    Valid <- c()
    for (i in 1:nrow(res)) {
        Valid[i] <- "+"
        tmp <- res[i,]
        if (anyNA(tmp)) {
            Valid[i] <- "-"
        }
    }
    return(Valid)
}

runDESeq2 <- function(parameters, groupInd, mdata, program, pair1, pair2, samples) {
    fitTypeInd <- singleChoiceParamValue(parameters, "Fit type")
    data_scale <- singleChoiceParamValue(parameters, "Data scale")
    data_raw <- main(mdata)
    if (data_scale == "Log2 scale") {
        multi <- intParamValue(parameters, "Multiplication for Log2 values")
        data_raw <- data_raw*multi
    } else {
        transform <- singleChoiceParamValue(parameters, "Transformation")
        if (transform == "Perform Log2") {
            data_raw <- log(data_raw, 2)
            multi <- intParamValue(parameters, "Multiplication")
            data_raw <- data_raw*multi
        } else if (transform == "Divided by constant") {
            divd <- intParamValue(parameters, "Divided by")
            data_raw <- data_raw/divd
        } else {
            forceRun <- boolParamValue(parameters, "Force to run")
            if (!forceRun) {
                if (max(data_raw) > 100000000) {
                    stop("Some values are higher than 1 billion, please transform your data. " + 
                    "You can also switch on Force to run, but error may occurs.", call. = FALSE)
                }
            }
        }
    }
    rounding <- boolParamValue(parameters, "Rounding")
    if (rounding)
        data_raw <- round(data_raw)
    checkInteger(data_raw)
    suppressMessages(library(DESeq2))
    sample <- data.frame(groups = samples)
    ds <- DESeqDataSetFromMatrix(countData = data_raw, colData = sample, design = ~groups)
    colnames(ds) <- colnames(counts)
    ds <- DESeq(ds, fitType = fitTypeInd)
    residual <- mcols(ds)$dispGeneEst - mcols(ds)$dispFit
    res <- results(ds, c('groups', pair2, pair1))
    res <- as.data.frame(res)
    Valid <- checkValid(res)
    combine <- generateOutput(parameters, res, program, pair1, pair2, Valid, "")
    return(combine)
}

runLimma_DEqMS <- function(parameters, groupInd, mdata, program, pair1, pair2, samples) {
    voomNeed <- boolParamValue(parameters, "Voom")
    trendNeed <- boolParamValue(parameters, "Trend")
    robustNeed <- boolParamValue(parameters, "Robust")
    spanValue <- intParamValue(parameters, "Span")
    normMethod <- singleChoiceParamValue(parameters, "Normalization method")
    DEqMSNeed <- boolParamValue(parameters, "DEqMS")
    library('limma')
    library('edgeR')
    if (robustNeed)
        library('statmod')
    data_raw <- main(mdata)
    mobDataGroups <- c()
    for (i in 1:length(samples)) {
        mobDataGroups[i] <- samples[[i]]
    }
    Group <- factor(mobDataGroups)
    design <- model.matrix(~0 + Group)
    if ((normMethod == "TMM") || (normMethod == "RLE") || (normMethod == "TMMwsp") || (normMethod == "upperquartile"))
    {
        d0 <- DGEList(data_raw)
        d0 <- calcNormFactors(d0, method = normMethod)
        data <- cpm(d0)
    } else {
        data <- normalizeBetweenArrays(data_raw, method = normMethod)
    }
    if (voomNeed) {
        data <- voom(data, design, span = spanValue)
    }
    fit <- lmFit(data, design)
    compareGroups <- paste('Group', pair2, '-', 'Group', pair1, sep = "")
    commandStr <- paste('makeContrasts(', compareGroups, ', levels = design)', sep = "")
    contrast.matrix <- eval(parse(text = commandStr))
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2, trend=trendNeed, robust=robustNeed)
    if (!DEqMSNeed) {
        res <- topTable(fit2, coef = 1, n = Inf)
    } else {
        library(dplyr)
        library(matrixStats)
        library(DEqMS)
        fitMethod <- singleChoiceParamValue(parameters, "Fit method")
        PSMstring <- singleChoiceParamValue(parameters, "Peptide/PSM counts")
        PSMstring <- paste(gsub("[^A-z0-9_]", ".", PSMstring), ".", sep = "")
        PSMcol <- colnames(annotCols(mdata))[grep(PSMstring, colnames(annotCols(mdata)))]
        if (length(PSMcol) == 0){
            PSMstring<-substr(PSMstring, 1, nchar(PSMstring)-1)
            PSMcol <- colnames(annotCols(mdata))[grep(PSMstring, colnames(annotCols(mdata)))]
        }
        PSM <- select(annotCols(mdata), PSMcol)
        pep.count.table = data.frame(count = rowMins(as.matrix(PSM)), row.names = rownames(data_raw))
        pep.count.table$count = pep.count.table$count+1
        fit2$count = pep.count.table[rownames(fit2$coefficients),"count"]
        t<-min(fit2$count)
        fit3 = spectraCounteBayes(fit2, fit.method=fitMethod)
        DEqMS.results = outputResult(fit3,coef_col = 1)
        DEqMS.results <- DEqMS.results[order(as.numeric(row.names(DEqMS.results))),]
        res<- data.frame(logFC=DEqMS.results$logFC, AveExpr=DEqMS.results$AveExpr,
                         t=DEqMS.results$sca.t, P.Value=DEqMS.results$sca.P.Value,
                         adj.pval=DEqMS.results$sca.adj.pval)
    }
    result_table <- as.data.frame(res)
    sort_result <- result_table[order(as.numeric(row.names(result_table))),]
    Valid <- checkValid(sort_result)
    combine <- generateOutput(parameters, sort_result, program, pair1, pair2, Valid, "")
    return(combine)
}

runSAM <- function(parameters, groupInd, mdata, program, pair1, pair2, samples) {
    if (singleChoiceParamValue(parameters, "Method selection") == "SAM")
        omicsType = "array"
    else
        omicsType = "seq"
    if (singleChoiceParamValue(parameters, "Test statistic") == "Standard (t - statistic)")
        tStat = "standard"
    else
        tStat = "wilcoxon"
    if (singleChoiceParamValue(parameters, "Time summary type") == "Slope")
        timeType = "slope"
    else
        timeType = "signed.area"
    if (singleChoiceParamValue(parameters, "Regression method") == "Standard (linear least squares)")
        regress = "standard"
    else
        regress = "ranks"
    knnNei <- intParamValue(parameters, "KNN neighbors")
    nResamp <- intParamValue(parameters, "Number of resamples")
    nResampPerm <- intParamValue(parameters, "Number of resamples for permutations")
    nPerm <- intParamValue(parameters, "Number of permutations")
    respType <- singleChoiceParamValue(parameters, "Data type")
    seedNum <- intParamValue(parameters, "Seed")
    logScale <- boolParamValue(parameters, "Log2")
    roundUp <- boolParamValue(parameters, "Rounding up")
    normMethod <- singleChoiceParamValue(parameters, "Normalization")
    library('samr')
    library('limma')
    data_raw <- main(mdata)
    data_raw <- normalizeBetweenArrays(data_raw, method = normMethod)
    set.seed(seedNum)
    mobDataGroups <- c()
    if (respType == "Two class unpaired"){
        existGroup <- c()
        start<-TRUE
        for (i in 1:length(samples)) {
            if (samples[[i]] %in% existGroup) {
                ind<-which(existGroup == samples[[i]])
                mobDataGroups <- append(mobDataGroups, ind)
            } else {
                existGroup <- append(existGroup, samples[[i]])
                if (!start) {
                    mobDataGroups <- append(mobDataGroups, 2)
                } else{
                    start<-FALSE
                    mobDataGroups <- append(mobDataGroups, 1)
                }
            }
        }
    }
    pairs<-c()
    if (respType == "Two class paired"){
        existGroup <- c()
        start<-TRUE
        nnum <- 3
        pnum <- 0
        for (i in 1:length(samples)) {
            if (samples[[i]] %in% existGroup) {
                ind<-which(existGroup == samples[[i]])
                if (ind == 1){
                   ind <- ind - nnum
                   nnum<-nnum+1
                } else {
                   ind <- ind + pnum
                   pnum<-pnum+1
                }
                mobDataGroups <- append(mobDataGroups, ind)
            } else {
                existGroup <- append(existGroup, samples[[i]])
                if (!start) {
                    mobDataGroups <- append(mobDataGroups, 1)
                } else{
                    start<-FALSE
                    mobDataGroups <- append(mobDataGroups, -1)
                }
            }
        }
    }
    Group <- factor(mobDataGroups)
    if ((omicsType == "seq") && (logScale)) {
        data_raw <- '^'(2,data_raw)
    }
    if ((omicsType == "seq") && (roundUp)) {
        data_raw <- round(data_raw)
    }
    data_raw<-as.matrix(data_raw)
    colnames(data_raw)<-NULL
    rownames(data_raw)<-NULL
    if (omicsType == "array") {
        data=list(x=data_raw,y=mobDataGroups, geneid=as.character(1:nrow(data_raw)),
                  genenames=paste("g",as.character(1:nrow(data_raw)),sep=""),
                  logged2=logScale)
    } else {
        data=list(x=data_raw,y=mobDataGroups, geneid=as.character(1:nrow(data_raw)),
                  genenames=paste("g",as.character(1:nrow(data_raw)),sep=""))
    }
    delta<-intParamValue(parameters, "Delta")
    log <- capture.output({samr.obj<-samr(data, resp.type=respType,
            nperms=nPerm, testStatistic=tStat,
            assay.type=omicsType, time.summary.type=timeType,
            regression.method=regress, knn.neighbors=knnNei, random.seed=seedNum,
            nresamp=nResamp, nresamp.perm=nResampPerm)})
    log <- capture.output({delta.table <- samr.compute.delta.table(samr.obj)})
    siggenes.table<-samr.compute.siggenes.table(samr.obj, delta, data, delta.table, all.genes=TRUE)
    up <- as.data.frame(siggenes.table$genes.up)
    down <- as.data.frame(siggenes.table$genes.lo)
    result<-rbind(up, down)
    sort_result <- result[order(as.numeric(result$Row)),]
    sort_result$Row<-NULL
    sort_result$'Gene ID'<-NULL
    sort_result$'Gene Name'<-NULL
    sort_result[] <- lapply(sort_result, as.numeric)
    sort_result[ncol(sort_result)] = sort_result[ncol(sort_result)]/100
    sort_result[ncol(sort_result)-1] = log(sort_result[ncol(sort_result)-1],2)
    colnames(sort_result)[ncol(sort_result)]<-'q.values'
    rownames(sort_result) <- 1:nrow(sort_result)
    Valid <- checkValid(sort_result)
    combine <- generateOutput(parameters, sort_result, program, pair1, pair2, Valid, "")
    return(combine)
}

runROTS <- function(parameters, groupInd, mdata, program, pair1, pair2, samples) {
    b <- intParamValue(parameters, "B")
    k <- stringParamValue(parameters, "K")
    paired <- boolParamValue(parameters, "Paired")
    logged <- boolParamValue(parameters, "Log2")
    seedNum <- intParamValue(parameters, "Seed")
    normalMethod <- singleChoiceParamValue(parameters, "Norm. method")
    library(ROTS)
    library('limma')
    data_raw <- main(mdata)
    data_raw <- normalizeBetweenArrays(data_raw, method = normalMethod)
    Group <- c()
    existGroup <- c()
    start<-TRUE
    for (i in 1:length(samples)) {
        if (samples[[i]] %in% existGroup) {
            ind<-which(existGroup == samples[[i]])
            if (ind == 2){
                ind <- ind - 2
            }
            Group <- append(Group, ind)
        } else {
            existGroup <- append(existGroup, samples[[i]])
            if (!start) {
                Group <- append(Group, 0)
            } else{
                start<-FALSE
                Group <- append(Group, 1)
            }
        }
    }
    if (k == "Default"){
        log <- capture.output({rots.out <- ROTS(data = data_raw, groups = Group,
                               B = b, seed = seedNum)})
    } else {
        k <- as.numeric(k)
        if (is.na(k))
            stop("K need to be an integer.", call.=FALSE)
        else if (((k*10) %% 10) > 0)
            stop("K need to be an integer.", call.=FALSE)
        log <- capture.output({rots.out <- ROTS(data = data_raw, groups = Group,
                              B = b, K = k , seed = seedNum)})
    }
    result<-data.frame(d=rots.out$d, logfc=rots.out$logfc, pvalue=rots.out$pvalue, FDR=rots.out$FDR)
    Valid <- checkValid(result)
    combine <- generateOutput(parameters, result, program, pair1, pair2, Valid, "")
}

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 3) {
  stop("Should provide two arguments: paramFile inFile outFile", call.=FALSE)
}

paramFile <- args[1]
inFile <- args[2]
outFile <- args[3]
library(PerseusR)
parameters <- parseParameters(paramFile)
program <- singleChoiceParamValue(parameters, "Program")
is.installed <- function(mypkg) is.element(mypkg, installed.packages()[, 1])
if ((program == "EdgeR") || (program == "Limma and DEqMS")) {
    if (!is.installed('edgeR')) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager", repos='http://cran.us.r-project.org')
        BiocManager::install("edgeR")
    }
    if (!is.installed('limma')) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager", repos='http://cran.us.r-project.org')
        BiocManager::install("limma")
    }
    if (!is.installed('statmod')) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager", repos='http://cran.us.r-project.org')
        BiocManager::install("statmod")
    }
    if (!is.installed('DEqMS')) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager", repos='http://cran.us.r-project.org')
        BiocManager::install("DEqMS")
    }
    if (!is.installed('matrixStats')) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager", repos='http://cran.us.r-project.org')
        BiocManager::install("matrixStats")
    }
    if (!is.installed('dplyr')) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager", repos='http://cran.us.r-project.org')
        BiocManager::install("dplyr")
    }
} else if (program == "DESeq2") {
    if (!is.installed('DESeq2')) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager", repos='http://cran.us.r-project.org')
        BiocManager::install("DESeq2")
    }
} else if (program == "SAM") {
    if (!is.installed('shiny')) {
        install.packages("shiny", dependencies = TRUE, repos='http://cran.us.r-project.org', quiet = TRUE)
    }
    if (!is.installed('GSA')) {
        install.packages("GSA", dependencies = TRUE, repos='http://cran.us.r-project.org', quiet = TRUE)
    }
    if (!is.installed('openxlsx')) {
        install.packages("openxlsx", dependencies = TRUE, repos='http://cran.us.r-project.org', quiet = TRUE)
    }
    if (!is.installed('shinyFiles')) {
        install.packages("shinyFiles", dependencies = TRUE, repos='http://cran.us.r-project.org', quiet = TRUE)
    }
    if (!is.installed('matrixStats')) {
        install.packages("matrixStats", dependencies = TRUE, repos='http://cran.us.r-project.org', quiet = TRUE)
    }
    if (!is.installed('impute')) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager", repos='http://cran.us.r-project.org')
        BiocManager::install("impute")
    }
    if (!is.installed('limma')) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager", repos='http://cran.us.r-project.org')
        BiocManager::install("limma")
    }
    if (!is.installed('samr')) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager", repos='http://cran.us.r-project.org')
        BiocManager::install("samr")
    }
} else if (program == "ROTS") {
    if (!is.installed('ROTS')) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager", repos='http://cran.us.r-project.org')
        BiocManager::install("ROTS")
    }
    if (!is.installed('limma')) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager", repos='http://cran.us.r-project.org')
        BiocManager::install("limma")
    }
}
mdata <- read.perseus(inFile)
group <- singleChoiceParamValue(parameters, "Group")
groupInd <- singleChoiceParamInd(parameters, "Group")
pair1 <- singleChoiceWithSubParamValue(paramFile, groupInd, "Reference sample", "Item")
pair2 <- singleChoiceWithSubParamValue(paramFile, groupInd, "Target sample", "Item")
fdr <- boolParamValue(parameters, "P-value")
if (pair1 == pair2) {
    stop("Please select different groups for comparison.", call. = FALSE)
}
reps <- list('1' = 0, '2' = 0)
samples <- c()
numSample <- 1
for (sampleInd in annotRows(mdata)[, groupInd + 1]) {
    samples[numSample] = sampleInd
    if (sampleInd == pair1) {
        reps[['1']] <- reps[['1']] + 1
    } else if (sampleInd == pair2) {
        reps[['2']] <- reps[['2']] + 1
    }
    numSample <- numSample + 1
}
noRep <- FALSE
for (rep in reps) {
    if (rep <= 1) {
        noRep <- TRUE
    }
}
if (program == "EdgeR") {
    testMethod <- singleChoiceParamValue(parameters, "Test method")
    normMethod <- singleChoiceParamValue(parameters, "Normalization method (EdgeR)")
    if ((noRep) || (testMethod == 'Exact Test')) {
        testMethod <- 'Exact Test'
        combine <- runEdgeRWithExactTest(parameters, groupInd, testMethod, mdata, program, pair1, pair2, normMethod)
    } else {
        combine <- runEdgeRWithGLM(parameters, groupInd, testMethod, mdata, program, pair1, pair2, samples, normMethod)
    }
} else if (program == "DESeq2") {
    combine <- runDESeq2(parameters, groupInd, mdata, program, pair1, pair2, samples)
} else if (program == "Limma and DEqMS") {
    if (noRep) {
        stop("Replicates are needed", call. = FALSE)
    } else {
        combine <- runLimma_DEqMS(parameters, groupInd, mdata, program, pair1, pair2, samples)
    }
} else if (program == "SAM") {
    if (noRep) {
        stop("Replicates are needed", call. = FALSE)
    } else {
        combine <- runSAM(parameters, groupInd, mdata, program, pair1, pair2, samples)
    }
} else if (program == "ROTS") {
    if (noRep) {
        stop("Replicates are needed", call. = FALSE)
    } else {
        combine <- runROTS(parameters, groupInd, mdata, program, pair1, pair2, samples)
    }
}
f_annoCols <- cbind(annotCols(mdata), combine)
outdata <- matrixData(main = main(mdata), annotCols = f_annoCols,
                      annotRows = annotRows(mdata))
print(paste('writing to', outFile))
write.perseus(outdata, outFile)