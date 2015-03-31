#03-31-2015
#C. Willis

########################################################################################
########################################################################################
# Parallel version of 'comdist' function from Picante
# clN = number of nodes to use for the cluster

comdist.pl = function (comm, dis, clN = 2, abundance.weighted = FALSE) 
{
	require(doParallel)
    dat <- match.comm.dist(comm, dis)
    x <- as.matrix(dat$comm)
    dis <- as.matrix(dat$dist)
    if (!abundance.weighted) {
        x <- decostand(x, method = "pa")
    }
    N <- dim(x)[1]
    S <- dim(x)[2]
    x <- decostand(x, method = "total", MARGIN = 1)
    comdist <- matrix(nrow = N, ncol = N)
    for (l in 1:N) {
    cl <- makeCluster(clN) # create parellel clusters
	registerDoParallel(cl)
    col = foreach(k = 1:N,.inorder=T,.combine='c') %dopar% {
            r = sum(dis * outer(as.vector(t(x[k,])), as.vector(t(x[l, ]))))
        	return(r)
        }
    comdist[,l] = col
    stopCluster(cl)
    }
    row.names(comdist) <- row.names(x)
    colnames(comdist) <- row.names(x)
    return(as.dist(comdist))
}

########################################################################################
#Parallelize version of 'comdistnt' function from Picante

comdistnt.pl =function (comm, dis,clN=2, abundance.weighted = FALSE, exclude.conspecifics = FALSE) 
{
	require(doParallel)
    dat <- match.comm.dist(comm, dis)
    comm <- dat$comm
    dis <- dat$dist
    N <- dim(comm)[1]
    comm <- decostand(comm, method = "total", MARGIN = 1)
    out <- matrix(nrow = N, ncol = N)
    for (i in 1:N) {
    cl <- makeCluster(clN) # create parellel clusters
	registerDoParallel(cl)

    col = foreach(j = 1:N,.inorder=T,.combine='c') %dopar% {
            sppInSample1 <- colnames(comm[i, comm[i, ] > 0, drop = FALSE])
            sppInSample2 <- colnames(comm[j, comm[j, ] > 0, drop = FALSE])
            if ((length(sppInSample1) >= 1) && (length(sppInSample2) >= 
                1)) {
                sample.dis <- dis[sppInSample1, sppInSample2, 
                  drop = FALSE]
                if (exclude.conspecifics) {
                  sample.dis[sample.dis == 0] <- NA
                }
                sample1NT <- apply(sample.dis, 1, min, na.rm = TRUE)
                sample1NT[sample1NT == Inf] <- NA
                sample2NT <- apply(sample.dis, 2, min, na.rm = TRUE)
                sample2NT[sample2NT == Inf] <- NA
                if (abundance.weighted) {
                  sample1.weights <- as.numeric(comm[i, sppInSample1])
                  sample2.weights <- as.numeric(comm[j, sppInSample2])
                  if (any(is.na(sample1NT))) {
                    miss <- which(is.na(sample1NT))
                    sample1NT <- sample1NT[-miss]
                    sample1.weights <- sample1.weights[-miss]
                    sample1.weights <- sample1.weights/sum(sample1.weights)
                  }
                  if (any(is.na(sample2NT))) {
                    miss <- which(is.na(sample2NT))
                    sample2NT <- sample2NT[-miss]
                    sample2.weights <- sample2.weights[-miss]
                    sample2.weights <- sample2.weights/sum(sample2.weights)
                  }
                  sampleNT <- c(sample1NT, sample2NT)
                  sample.weights <- c(sample1.weights, sample2.weights)
                  #ORIGINAL: comdisnt[i, j] <- weighted.mean(sampleNT, sample.weights, na.rm = TRUE)
                  r <- weighted.mean(sampleNT, sample.weights, na.rm = TRUE)
                  return(r)
                }
                else {
                  #ORIGINAL: comdisnt[i, j] <- mean(c(sample1NT, sample2NT), na.rm = TRUE)
                  r = comdisnt[i, j] <- mean(c(sample1NT, sample2NT), na.rm = TRUE)
                  return(r)
                }
            }
            else {
                r <- NA
                return(r)
            }
        } # end foreach loop   
    	out[,i] = col
    	stopCluster(cl)
    }
    rownames(out) = colnames(out) <- rownames(comm)
    return(as.dist(t(out)))
}

########################################################################################
########################################################################################
# parallelized 'pd' function from picante

pd.pl = function (samp, tree, clN=2, include.root = TRUE) {

require(doParallel)

    if (is.null(tree$edge.length)) {
        stop("Tree has no branch lengths, cannot compute pd")
    }
    species <- colnames(samp)
    SR <- rowSums(ifelse(samp > 0, 1, 0))
    nlocations = dim(samp)[1]
    nspecies = dim(samp)[2]
    #PDs = NULL
    
    cl <- makeCluster(clN) # create parellel clusters
	registerDoParallel(cl)
	
    #for (i in 1:nlocations) {
    PDs = foreach(i = 1:nlocations,.packages='picante',.inorder=T,.combine='rbind') %dopar% {
        present <- species[samp[i, ] > 0]
        treeabsent <- tree$tip.label[which(!(tree$tip.label %in% 
            present))]
        if (length(present) == 0) {
            PD <- 0
            return(PD)
        }
        else if (length(present) == 1) {
            if (!is.rooted(tree) || !include.root) {
                warning("Rooted tree and include.root=TRUE argument required to calculate PD of single-species sampunities. Single species sampunity assigned PD value of NA.")
                PD <- NA
                return(PD)
            }
            else {
                PD <- node.age(tree)$ages[which(tree$edge[, 2] == 
                  which(tree$tip.label == present))]
                return(PD)
            }
        }
        else if (length(treeabsent) == 0) {
            PD <- sum(tree$edge.length)
            return(PD)
        }
        else {
            sub.tree <- drop.tip(tree, treeabsent)
            if (include.root) {
                if (!is.rooted(tree)) {
                  stop("Rooted tree required to calculate PD with include.root=TRUE argument")
                }
                sub.tree.depth <- max(node.age(sub.tree)$ages)
                orig.tree.depth <- max(node.age(tree)$ages[which(tree$edge[, 
                  2] %in% which(tree$tip.label %in% present))])
                PD <- sum(sub.tree$edge.length) + (orig.tree.depth - 
                  sub.tree.depth)
            return(PD)	
            }
            else {
                PD <- sum(sub.tree$edge.length)
                return(PD)
            }
        }
        #PDs <- c(PDs, PD)
        }
	stopCluster(cl)
	
    PDout <- data.frame(PD = PDs, SR = SR)
    rownames(PDout) <- rownames(samp)
    return(PDout)
}


########################################################################################
########################################################################################
# ses.pd function parallelized
# requires pd.pl to run, see above

 ses.pd.pl = function (samp, tree, cl.size=2,null.model = c("taxa.labels", "richness", 
    "frequency", "sample.pool", "phylogeny.pool", "independentswap", 
    "trialswap"), runs = 999, iterations = 1000, ...) 
{
    pd.obs <- as.vector(pd.pl(samp, tree,clN=cl.size, ...)$PD)
    null.model <- match.arg(null.model)
    pd.rand <- switch(null.model, taxa.labels = t(replicate(runs, 
        as.vector(pd.pl(samp, tipShuffle(tree),clN=cl.size, ...)$PD))), richness = t(replicate(runs, 
        as.vector(pd.pl(randomizeMatrix(samp, null.model = "richness"), 
            tree,clN=cl.size, ...)$PD))), frequency = t(replicate(runs, as.vector(pd.pl(randomizeMatrix(samp, 
        null.model = "frequency"), tree,clN=cl.size, ...)$PD))), sample.pool = t(replicate(runs, 
        as.vector(pd.pl(randomizeMatrix(samp, null.model = "richness"), 
            tree,clN=cl.size, ...)$PD))), phylogeny.pool = t(replicate(runs, 
        as.vector(pd.pl(randomizeMatrix(samp, null.model = "richness"), 
            tipShuffle(tree),clN=cl.size, ...)$PD))), independentswap = t(replicate(runs, 
        as.vector(pd.pl(randomizeMatrix(samp, null.model = "independentswap", 
            iterations), tree,clN=cl.size, ...)$PD))), trialswap = t(replicate(runs, 
        as.vector(pd.pl(randomizeMatrix(samp, null.model = "trialswap", 
            iterations), tree,clN=cl.size, ...)$PD))))
    pd.rand.mean <- apply(X = pd.rand, MARGIN = 2, FUN = mean, 
        na.rm = TRUE)
    pd.rand.sd <- apply(X = pd.rand, MARGIN = 2, FUN = sd, na.rm = TRUE)
    pd.obs.z <- (pd.obs - pd.rand.mean)/pd.rand.sd
    pd.obs.rank <- apply(X = rbind(pd.obs, pd.rand), MARGIN = 2, 
        FUN = rank)[1, ]
    pd.obs.rank <- ifelse(is.na(pd.rand.mean), NA, pd.obs.rank)
    data.frame(ntaxa = specnumber(samp), pd.obs, pd.rand.mean, 
        pd.rand.sd, pd.obs.rank, pd.obs.z, pd.obs.p = pd.obs.rank/(runs + 
            1), runs = runs, row.names = row.names(samp))
}
