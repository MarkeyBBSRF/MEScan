
################################################################################

#The log-likelihood function

funLogLikelihood <- function(parVec, mutationMat, q = NULL){
	gamma <- parVec[1]
	delta <- parVec[-1]
	M <- ncol(mutationMat)
	if(is.null(q))q <- as.vector(table(factor(rowSums(mutationMat), level = 0:M)))
	if(length(delta) == 1){
		J <- 1:M
		logL <- M * q[1] * log(1 - delta) + q[1] * log(1 - gamma) + sum(q[-1] * (log(J / M) + (J - 1) * log(delta) + (M - J) * log(1 - delta) + log(M / J * delta * (1 - gamma) + gamma)))
	}else if(length(delta) == M){
		tempA <- matrix(NA, 2^M, M)
		for(i in 1:M)tempA[, i] <- rep(rep(c(FALSE, TRUE), each = 2^(M - i)), 2^(i - 1))
		tempProbMat <- matrix(rep(delta, each = 2^M), 2^M, M)
		tempProbMat[!tempA] <- 1 - tempProbMat[!tempA]
		tempProd <- apply(tempProbMat, 1, prod)
		tempMat <- rep(delta / sum(delta), each = 2^M) * tempA * (tempProd / tempProbMat)
		tempProb <- (1 - gamma) * tempProd + gamma * rowSums(tempMat)
		prob <- tempProb[colSums(t(mutationMat) * 2^((M - 1):0)) + 1]
		logL <- sum(log(prob))
	}else stop("Invalid Value of delta")
	return(-logL)
}

################################################################################

#The gradient function of the log-likelihood function, called in the optim function for speeding up

funGradient <- function(parVec, mutationMat, q = NULL){
	gamma <- parVec[1]
	delta <- parVec[-1]
	M <- ncol(mutationMat)
	if(is.null(q))q <- as.vector(table(factor(rowSums(mutationMat), level = 0:M)))
	if(length(delta) == 1)delta <- rep(delta, M)
	if(length(delta) != M)stop("Invalid Value of delta")
	sumDelta <- sum(delta)
	tempA <- matrix(NA, 2^M, M)
	for(i in 1:M)tempA[, i] <- rep(rep(c(FALSE, TRUE), each = 2^(M - i)), 2^(i - 1))
	tempProbMat1 <- tempProbMat2 <- matrix(delta, 2^M, M, byrow = TRUE)
	tempProbMat2[!tempA] <- 1 - tempProbMat1[!tempA]
	temp1 <- apply(tempProbMat2, 1, prod)
	temp2 <- temp1 / tempProbMat2
	tempProbMat3 <- tempProbMat1 / sumDelta * tempA * temp2
	temp3 <- rowSums(tempProbMat3)
	tempProb <- (1 - gamma) * temp1 + gamma * temp3
	tempGradGamma <- (temp3 - temp1) / tempProb
	tempGradMat <- matrix(NA, 2^M, M)
	temp4 <- (1 - gamma) * (tempA * 2 - 1) * temp2
	temp5 <- (sumDelta - delta) / sumDelta^2
	temp6 <- rep(temp5, each = 2^M) * tempA * temp2
	for(i in 1:M){
		temp7 <- -tempProbMat1[, -i, drop = FALSE] / sumDelta^2 * tempA[, -i, drop = FALSE] * temp2[, -i, drop = FALSE] + tempProbMat1[, -i, drop = FALSE] / sumDelta * tempA[, -i, drop = FALSE] * (temp2[, -i, drop = FALSE] / tempProbMat2[, i] * (tempA[, i] * 2 - 1))
		tempGradMat[, i] <- temp4[, i] + gamma * (temp6[, i] + rowSums(temp7))
	}
	tempGradMat <- tempGradMat / tempProb
	tempIdx <- colSums(t(mutationMat) * 2^((M - 1):0)) + 1
	grad <- c(sum(tempGradGamma[tempIdx]), colSums(tempGradMat[tempIdx, ]))
	return(-grad)
}

################################################################################

#Function for the maximum likelihood estimation of the parameters

funEstimate <- function(mutationMat, tol = 1e-7){
	N <- nrow(mutationMat)
	M <- ncol(mutationMat)
	tempRowSums <- rowSums(mutationMat)
	q <- as.vector(table(factor(tempRowSums, level = 0:M)))
	piHat <- colMeans(mutationMat)
	probMat <- matrix(rep(piHat, each = N), N, M)
	tempMat <- log(probMat^mutationMat * (1 - probMat)^(!mutationMat))
	logL0Each <- rowSums(tempMat)
	logL0 <- sum(logL0Each)
	for(eps in 10^(-((-log10(tol)):3))){
		tempOptim <- try(optim(c(mean(tempRowSums > 0), piHat), funLogLikelihood, gr = funGradient, method = "L-BFGS-B", lower = c(0, rep(eps, M)), upper = rep(1 - eps, M + 1), mutationMat = mutationMat, q = q), silent = TRUE)
		if(class(tempOptim) == "list" && tempOptim$convergenc == 0 && all(tempOptim$par > 0) && all(tempOptim$par < 1))break
	}
	if(class(tempOptim) == "try-error")stop("ERROR in OPTIM FUNCTION!")
	logL1 <- -tempOptim$value
	gammaHat <- tempOptim$par[1]
	deltaHat <- tempOptim$par[-1]
	if(logL1 < logL0){gammaHat <- 0; deltaHat <- piHat; logL1 <- logL0}
	S <- -2 * (logL0 - logL1)
	return(list(pi = piHat, gamma = gammaHat, delta = deltaHat, logL0 = logL0, logL1 = logL1, S = S))
}

################################################################################

#Function for adding one gene to a group of genes based on maximizing the LRT statistic

funAdd1 <- function(geneStart, mutationMat, detail = TRUE){
	geneAll <- colnames(mutationMat)
	M0 <- length(geneStart)
	if(M0 < 1)stop("Must start from at least 1 gene!")
	if(any(!geneStart %in% geneAll))stop("Invalid gene names!")
	geneCandidate <- geneAll[!geneAll %in% geneStart]
	M1 <- length(geneCandidate)
	if(M1 < 1)stop("Candidate gene set must have at least one gene that is not in the start set!")
	S <- rep(NA, M1)
	for(i in 1:M1){
		S[i] <- funEstimate(mutationMat[, c(geneStart, geneCandidate[i])])$S
		if(detail)cat(i, "\t", geneStart, "+", geneCandidate[i], S[i], "\n")
	}
	tempIdx <- which.max(S)
	geneEnd <- c(geneStart, geneCandidate[tempIdx])
	return(list(gene = geneEnd, S = S[tempIdx]))
}

################################################################################

#Funtion for maximizing the LRT statistic based on forward steps from the most significant gene pairs

funMaxS <- function(mutationMat, nPairStart = 10, maxSize = 6, detail = TRUE){
	geneAll <- colnames(mutationMat)
	M <- ncol(mutationMat)
	tempCombn <- combn(M, 2)
	nPair <- ncol(tempCombn)
	pairS <- rep(NA, nPair); names(pairS) <- 1:nPair
	for(i in 1:nPair){
		pairS[i] <- funEstimate(mutationMat[, tempCombn[, i]])$S
		if(detail)cat(i, "\t", geneAll[tempCombn[, i]], pairS[i], "\n")
	}
	startPairS <- sort(pairS, decreasing = T)[1:nPairStart]
	startCombn <- tempCombn[, as.integer(names(startPairS))]
	maxSMat <- matrix(NA, nPairStart, maxSize - 1, dimnames = list(1:nPairStart, 2:maxSize))
	geneMat <- matrix(NA, nPairStart, maxSize, dimnames = list(1:nPairStart, 1:maxSize))
	maxSMat[, "2"] <- startPairS
	geneMat[, "1"] <- geneAll[startCombn[1, ]]
	geneMat[, "2"] <- geneAll[startCombn[2, ]]
	for(i in 1:nPairStart){
		tempGeneSet <- geneAll[startCombn[, i]]
		for(j in 2:(maxSize - 1)){
			tempAdd <- funAdd1(tempGeneSet, mutationMat, detail = detail)
			tempGeneSet <- tempAdd$gene
			maxSMat[i, j] <- tempAdd$S
			geneMat[i, j + 1] <- tempGeneSet[j + 1]
		}
	}
	return(list(maxSMat = maxSMat, geneMat = geneMat, startCombn = startCombn, startPairS = startPairS))
}

################################################################################

#Function for calculating the maximum LRT statistic of the permutated mutation matrix

funMaxSSimu <- function(mutationMat, nSimu = 1, nPairStart = 10, maxSize = 3, detail = TRUE){
	N <- nrow(mutationMat)
	M <- ncol(mutationMat)
	maxSArray <- array(NA, c(nPairStart, maxSize - 1, nSimu))
	for(s in 1:nSimu){
		permutationMat <- mutationMat
		for(j in 1:M)permutationMat[, j] <- mutationMat[sample(N), j]
		maxSArray[, , s] <- funMaxS(permutationMat, nPairStart = nPairStart, maxSize = maxSize, detail = detail)$maxSMat
		if(detail)cat("Simulation", s, "done", "\n")
	}
	maxSSimu <- apply(maxSArray, 3:2, max)
}

################################################################################

#Function for the global test whether MEGS exists in the mutation matrix

funGlobalTest <- function(mutationMat, maxSSimu = NULL, nSimu = 1000, nPairStart = 10, maxSize = 6, detail = TRUE){
      N <- nrow(mutationMat)
	M <- ncol(mutationMat)
	if(!is.null(maxSSimu)){
		nSimu <- nrow(maxSSimu)
		maxSize <- ncol(maxSSimu) + 1
	}else{
		maxSSimu <- funMaxSSimu(mutationMat, nSimu = nSimu, nPairStart = nPairStart, maxSize = maxSize, detail = detail)
	}
	resultReal <- funMaxS(mutationMat, nPairStart = nPairStart, maxSize = maxSize, detail = detail)
	maxSMatReal <- resultReal$maxSMat
	maxSReal <- apply(maxSMatReal, 2, max)
	rankMat <- apply(-rbind(maxSReal, maxSSimu), 2, rank, ties.method = "max")
	minRankVecReal <- min(rankMat[1, ] - 1)
	minRankVecSimu <- apply(apply(-maxSSimu, 2, rank, ties.method = "max"), 1, min)
	p <- mean(minRankVecSimu <= minRankVecReal)
	return(list(p = p, q = minRankVecSimu / nSimu, maxSReal = maxSReal, maxSSimu = maxSSimu, startPairIdx = resultReal$startCombn, startPairS = resultReal$startPairS))
}

################################################################################

#function for adding one gene to a group of genes based on q values

funAdd1RealData <- function(geneStartList, mutationMat, maxSSimu, eps = 0.005, level = 0.05, nPerm = 1e2, detail = TRUE){
	maxSize <- ncol(maxSSimu) + 1
	geneAll <- colnames(mutationMat)
	N <- nrow(mutationMat)
	geneEndList <- list(); k <- 1
	for(i in 1:length(geneStartList)){
		geneStart <- geneStartList[[i]]
		n0 <- length(geneStart)
		if(is.null(attr(geneStart, "S")) || is.null(attr(geneStart, "q"))){
			tempS <- funEstimate(mutationMat[, geneStart])$S
			tempQ <- mean(maxSSimu[, n0 - 1] >= tempS)
			attr(geneStart, "S") <- tempS
			attr(geneStart, "q") <- tempQ
		}
		if(n0 == maxSize)attr(geneStart, "change") <- FALSE
		if(!is.null(attr(geneStart, "change")) && !attr(geneStart, "change")){
			geneEnd <- geneStart
			geneEndList[[k]] <- geneEnd; k <- k + 1
			next
		}
		S0 <- attr(geneStart, "S")
		q0 <- attr(geneStart, "q")
		geneCandidate <- geneAll[!geneAll %in% geneStart]
		n1 <- length(geneCandidate)
		S1Vec <- q1Vec <- rep(NA, n1)
		for(j in 1:n1){
			S1Vec[j] <- funEstimate(mutationMat[, c(geneStart, geneCandidate[j])])$S
			q1Vec[j] <- mean(maxSSimu[, n0] >= S1Vec[j])
			if(detail)cat(i, j, geneStart, S0, q0, "+", geneCandidate[j], S1Vec[j], q1Vec[j], "\n")
		}
		tempIdx1 <- which(abs(q1Vec - q0) < eps)
		tempIdx2 <- which(q1Vec <= q0 - eps)
		if(length(tempIdx1) > 0){
			SPermMat <- matrix(NA, nPerm, n1)
			maxSPermVec <- rep(NA, nPerm)
			for(ii in 1:nPerm){
				tempPermIdx <- sample(N)
				for(jj in 1:n1){
					APerm <- mutationMat[, c(geneStart, geneCandidate[jj])]; APerm[, n0 + 1]  <- APerm[tempPermIdx, n0 + 1]
					SPermMat[ii, jj] <- funEstimate(APerm)$S
				}
				maxSPermVec[ii] <- max(SPermMat[ii, ])
				pMin <- mean(maxSPermVec >= max(S1Vec[tempIdx1]), na.rm = T)
#				if(detail)cat("Permutation", ii, ":\tmax(S1) =", max(S1Vec[tempIdx1]), "\tmax(SPerm) =", maxSPermVec[ii], "\tpMin =", pMin, "\n")
				if(ii >= nPerm / 5 && pMin > 2 * level)break
			}
			pVec <- rep(NA, length(tempIdx1))
			for(ii in 1:length(tempIdx1))pVec[ii] <- mean(maxSPermVec >= S1Vec[tempIdx1[ii]], na.rm = TRUE)
#			if(detail)cat(pVec, "\n")
			tempIdx2 <- sort(c(tempIdx2, tempIdx1[which(pVec <= level)]))
		}
		if(length(tempIdx2) == 0){
			geneEnd <- geneStart
			attr(geneEnd, "change") <- FALSE
			geneEndList[[k]] <- geneEnd; k <- k + 1
		}else{
			for(j in 1:length(tempIdx2)){
				geneEnd <- c(geneStart, geneCandidate[tempIdx2[j]])
				attr(geneEnd, "S") <- S1Vec[tempIdx2[j]]
				attr(geneEnd, "q") <- q1Vec[tempIdx2[j]]
				attr(geneEnd, "change") <- TRUE
				geneEndList[[k]] <- geneEnd; k <- k + 1
			}
		}
	}
	return(geneEndList)
}

################################################################################

#function for droping one gene from a group of genes based on q values

funDrop1RealData <- function(geneStartList, mutationMat, maxSSimu, eps = 0.005, level = 0.05, nPerm = 1e2, permResult = list(), detail = TRUE){
	minSize <- 2
	geneAll <- colnames(mutationMat)
	N <- nrow(mutationMat)
	geneEndList <- list(); k <- 1
	for(i in 1:length(geneStartList)){
		geneStart <- geneStartList[[i]]
		n0 <- length(geneStart)
		if(is.null(attr(geneStart, "S")) || is.null(attr(geneStart, "q"))){
			tempS <- funEstimate(mutationMat[, geneStart])$S
			tempQ <- mean(maxSSimu[, n0 - 1] >= tempS)
			attr(geneStart, "S") <- tempS
			attr(geneStart, "q") <- tempQ
		}
		if(n0 == minSize)attr(geneStart, "change") <- FALSE
		if(!is.null(attr(geneStart, "change")) && !attr(geneStart, "change")){
			geneEnd <- geneStart
			geneEndList[[k]] <- geneEnd; k <- k + 1
			next
		}
		S0 <- attr(geneStart, "S")
		q0 <- attr(geneStart, "q")
		S1Vec <- q1Vec <- rep(NA, n0)
		for(j in 1:n0){
			S1Vec[j] <- funEstimate(mutationMat[, geneStart[-j]])$S
			q1Vec[j] <- mean(maxSSimu[, n0 - 2] >= S1Vec[j])
			if(detail)cat(i, j, geneStart, S0, q0, "-", geneStart[j], S1Vec[j], q1Vec[j], "\n")
		}
		tempIdx1 <- which(abs(q1Vec - q0) < eps)
		tempIdx2 <- which(q1Vec <= q0 - eps)
		if(length(tempIdx1) > 0){
			pVec <- rep(NA, length(tempIdx1))
			for(j in 1:length(tempIdx1)){
				tempName <- paste(sort(geneStart[-tempIdx1[j]]), collapse = "_")
				if(!is.null(permResult[[tempName]]) && length(permResult[[tempName]]) >= nPerm){
					pVec[j] <- mean(permResult[[tempName]] >= S0, na.rm = TRUE)
				}else{
					geneCandidate <- geneAll[!geneAll %in% geneStart[-tempIdx1[j]]]
					n1 <- length(geneCandidate)
					SPermMat <- matrix(NA, nPerm, n1)
					maxSPermVec <- rep(NA, nPerm)
					for(ii in 1:nPerm){
						tempPermIdx <- sample(N)
						for(jj in 1:n1){
							APerm <- mutationMat[, c(geneStart[-tempIdx1[j]], geneCandidate[jj])]; APerm[, n0]  <- APerm[tempPermIdx, n0]
							SPermMat[ii, jj] <- funEstimate(APerm)$S
						}
						maxSPermVec[ii] <- max(SPermMat[ii, ])
						p <- mean(maxSPermVec >= S0, na.rm = T)
#						if(detail)cat(i, j, "Permutation", ii, ":\tS0 =", S0, "\tmax(SPerm) =", maxSPermVec[ii], "\tp =", p, "\n")
						if(ii >= nPerm / 5 && p > 2 * level)break
					}
					pVec[j] <- p
					permResult[[tempName]] <- c(permResult[[tempName]], maxSPermVec[!is.na(maxSPermVec)])
				}
			}
#			if(detail)cat(pVec, "\n")
			tempIdx2 <- sort(c(tempIdx2, tempIdx1[which(pVec > level)]))
		}
		if(length(tempIdx2) == 0){
			geneEnd <- geneStart
			attr(geneEnd, "change") <- FALSE
			geneEndList[[k]] <- geneEnd; k <- k + 1
		}else{
			for(j in 1:length(tempIdx2)){
				geneEnd <- geneStart[-tempIdx2[j]]
				attr(geneEnd, "S") <- S1Vec[tempIdx2[j]]
				attr(geneEnd, "q") <- q1Vec[tempIdx2[j]]
				attr(geneEnd, "change") <- ifelse(n0 > 3, TRUE, FALSE)
				geneEndList[[k]] <- geneEnd; k <- k + 1
			}
		}
	}
	attr(geneEndList, "permResult") <- permResult
	return(geneEndList)
}

################################################################################

#function for select the significant MEGS

funSelect <- function(mutationMat, maxSSimu = NULL, nSimu = 1000, nPairStart = 10, maxSize = 3, level = 0.05, detail = TRUE){
	geneAll <- colnames(mutationMat)
	M <- ncol(mutationMat)
	test <- funGlobalTest(mutationMat, maxSSimu = maxSSimu, nSimu = nSimu, nPairStart = nPairStart, maxSize = maxSize, detail = detail)
	if(test$p > level){
		MEGSList <- list()
	}else{
		MEGSList <- list(); length(MEGSList) <- nPairStart
		for(i in 1:nPairStart){
			tempGeneSet <- geneAll[test$startPairIdx[, i]]
			attr(tempGeneSet, "S") <- test$startPairS[i]
			attr(tempGeneSet, "q") <- mean(maxSSimu[, 1] >= test$startPairS[i])
			attr(tempGeneSet, "change") <- TRUE
			MEGSList[[i]] <- tempGeneSet
		}
		for(i in 2:(maxSize - 1)){
			MEGSList <- funAdd1RealData(MEGSList, mutationMat, maxSSimu, detail = detail)
			MEGSList <- MEGSList[!duplicated(unlist(lapply(MEGSList, function(x)paste(sort(x), collapse = "_"))))]
		}
		permResult <- list()
		for(i in maxSize:(2 + 1)){
			MEGSList <- funDrop1RealData(MEGSList, mutationMat, maxSSimu, detail = detail, permResult = permResult)
			permResult <- attr(MEGSList, "permResult")
			MEGSList <- MEGSList[!duplicated(unlist(lapply(MEGSList, function(x)paste(sort(x), collapse = "_"))))]
#			if(detail)cat(i, length(MEGSList), "\n")
		}
		tempFun <- function(x){
			attr(x, "pNominal") <- 0.5 * pchisq(attr(x, "S"), 1, lower.tail = FALSE)
			attr(x, "pCorrected") <- mean(test$q <= attr(x, "q"))
			attr(x, "coverage") <- mean(rowSums(mutationMat[, x]) > 0)
			attr(x, "q") <- NULL
			attr(x, "change") <- NULL
			return(x)
		}
		MEGSList <- lapply(MEGSList, tempFun)
		MEGSList <- MEGSList[order(unlist(lapply(MEGSList, attr, which = "pCorrected")))]
		MEGSList <- MEGSList[unlist(lapply(MEGSList, attr, which = "pCorrected")) <= level]
	}
	return(list(p = test$p, MEGSList = MEGSList))
}

################################################################################

#Function for converting the MEGS list to a data frame and write it to a specify file

funPrintMEGS <- function(MEGSList, outputFile = NULL){
	gene <- unlist(lapply(MEGSList, paste, collapse = " "))
	LRT <- unlist(lapply(MEGSList, attr, which = "S"))
	pNominal <- unlist(lapply(MEGSList, attr, which = "pNominal"))
	pCorrected <- unlist(lapply(MEGSList, attr, which = "pCorrected"))
	coverage <- unlist(lapply(MEGSList, attr, which = "coverage"))
	MEGSDF <- data.frame(gene = gene, coverage = coverage, LRT = LRT, pNominal = pNominal, pCorrected = pCorrected)
	if(!is.null(outputFile))write.table(MEGSDF, file = outputFile, sep = "\t", quote = FALSE, row.names = FALSE)
	return(MEGSDF)
}

################################################################################

#Function for Visualization

funPlotMutationMat <- function(mutationMat, main = NULL, reorder = TRUE){
	N <- nrow(mutationMat)
	M <- ncol(mutationMat)
	if(reorder){
		mutationMat <- mutationMat[, order(colSums(mutationMat), decreasing = TRUE)]
		for(i in ncol(mutationMat):1)mutationMat <- mutationMat[order(mutationMat[, i]),]
	}
	image(1:M, 0:N, t(!mutationMat), axes = F, xlim = c(0.5, M + 0.5), ylim = c(0, N), main = main, xlab = "", ylab = "")
	abline(v = 0:M + 0.5, h = c(0, N))
	axis(1, 1:M, colnames(mutationMat))
	axis(2)
}

################################################################################

#Function for generating the figures of the significant MEGS

funPlotMEGS <- function(MEGSList, mutationMat, outputDir, type = "pdf"){
	if(!file.exists(outputDir))dir.create(outputDir)
	if(type == "pdf"){
		for(i in 1:length(MEGSList)){
			outputFile <- paste0(outputDir, "/", "gene_set_", i, ".pdf")
			pdf(outputFile, width = 5, height = 10)
			par(cex.main = 1.5, cex.axis = 1.2)
			funPlotMutationMat(mutationMat[, MEGSList[[i]]], main = paste0("Coverage = ", signif(attr(MEGSList[[i]], "coverage"), 3), "\n", "Corrected p = ", signif(attr(MEGSList[[i]], "pCorrected"), 3)))
			dev.off()
		}
	}else if(type == "png"){
		for(i in 1:length(MEGSList)){
			outputFile <- paste0(outputDir, "/", "gene_set_", i, ".png")
			png(outputFile, width = 600, height = 1200)
			par(cex.main = 1.5, cex.axis = 1.2)
			funPlotMutationMat(mutationMat[, MEGSList[[i]]], main = paste0("Coverage = ", signif(attr(MEGSList[[i]], "coverage"), 3), "\n", "Corrected p = ", signif(attr(MEGSList[[i]], "pCorrected"), 3)))
			dev.off()
		}
	}else if(type == "jpg"){
		for(i in 1:length(MEGSList)){
			outputFile <- paste0(outputDir, "/", "gene_set_", i, ".jpg")
			jpeg(outputFile, width = 600, height = 1200)
			par(cex.main = 1.5, cex.axis = 1.2)
			funPlotMutationMat(mutationMat[, MEGSList[[i]]], main = paste0("Coverage = ", signif(attr(MEGSList[[i]], "coverage"), 3), "\n", "Corrected p = ", signif(attr(MEGSList[[i]], "pCorrected"), 3)))
			dev.off()
		}
	}else stop("Invalid value of the file type!")
}

################################################################################

#The main function of MEGSA

funMEGSA <- function(mutationMatFile, maxSSimuFile = NULL, resultTableFile = NULL, figureDir = NULL, nSimu = 1000, nPairStart = 10, maxSize = 6, level = 0.05, detail = TRUE, type = "pdf"){
	mutationMat <- as.matrix(read.table(mutationMatFile, header = TRUE, row.names = 1) != 0)
	if(is.null(maxSSimuFile)){
		maxSSimu <- NULL
	}else maxSSimu <- as.matrix(read.table(maxSSimuFile))
	resultMEGSA <- funSelect(mutationMat, maxSSimu, nSimu = nSimu, nPairStart = nPairStart, maxSize = maxSize, level = level, detail = detail)
	MEGSList <- resultMEGSA$MEGSList
	if(!is.null(resultTableFile))funPrintMEGS(MEGSList, outputFile = resultTableFile)
	if(!is.null(figureDir))funPlotMEGS(MEGSList, mutationMat, outputDir = figureDir, type = type)
	return(resultMEGSA)
}

################################################################################


