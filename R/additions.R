#' Compute maximum correlation for two images, taking into account a sequence of
#' rotation angles and all integer-valued translations (updated)
#'
#' This is an updated version of calculateCCFmax(). It normalizes the images
#' after rotation (mean 0, variance 1), and uses a slightly different grid
#' search to reduce runtime.
#'
#' @param image1 first image matrix
#' @param image2 second image matrix. The two images do not need to be the same
#'   size.
#'
#' @return A list with four items: 1) maximum correlation taking into account
#'   many possible rotations and translations of the second image, 2) the
#'   corresponding dx, the horizontal translation, 3) dy, the vertical
#'   translation, and 4) theta, the corresponding rotation angle.
#'
#' @export
calculateCCFmaxSearch <- function(image1, image2) {
    image1_small <- EBImage::resize(image1, w = floor(dim(image1)[1]/4),
        h = floor(dim(image1)[2]/4))
    image2_small <- EBImage::resize(image2, w = floor(dim(image2)[1]/4),
        h = floor(dim(image2)[2]/4))
    thetas <- c(seq(from = -175, to = 180, by = 5), -7.5, -2.5, 2.5, 7.5)
    allResults <- data.frame(thetas = thetas, dx = NA, dy = NA,
        corr = NA)
    for (i in 1:length(thetas)) {
        rotated <- cartridges3D:::bilinearInterpolation(image2_small,
            thetas[i])
        mean1 <- mean(image1_small[image1_small != 0])
        image1_small[image1_small != 0] <- image1_small[image1_small !=
            0] - mean1
        mean2 <- mean(rotated[rotated != 0])
        rotated[rotated != 0] <- rotated[rotated != 0] - mean2
        fact <- sqrt(sum(image1_small^2))
        image1_small <- image1_small/fact
        fact <- sqrt(sum(rotated^2))
        rotated <- rotated/fact
        out <- cartridges3D:::comparison(image1_small, rotated)
        allResults[i, ] <- c(thetas[i], out$dx, out$dy, out$corr)
    }
    tmp <- thetas[which.max(allResults$corr)]
    if (tmp %in% seq(from = -10, to = 10, by = 2.5)) {
        fineThetas <- seq(from = tmp - 2, to = tmp + 2, by = 0.5)
    } else {
        fineThetas <- seq(from = tmp - 4, to = tmp + 4, by = 1)
    }
    fineResults <- data.frame(fineThetas = fineThetas, dx = NA,
        dy = NA, corr = NA)
    for (i in 1:length(fineThetas)) {
        rotated <- cartridges3D:::bilinearInterpolation(image2_small,
            fineThetas[i])
        mean1 <- mean(image1_small[image1_small != 0])
        image1_small[image1_small != 0] <- image1_small[image1_small !=
            0] - mean1
        mean2 <- mean(rotated[rotated != 0])
        rotated[rotated != 0] <- rotated[rotated != 0] - mean2
        out <- cartridges3D:::comparison(image1_small, rotated)
        fineResults[i, ] <- c(fineThetas[i], out$dx, out$dy,
            out$corr)
    }
    index <- which.max(fineResults$corr)
    ret <- list(corr = fineResults$corr[index], dx = fineResults$dx[index],
        dy = fineResults$dy[index], theta = fineResults$fineThetas[index])
    return(ret)
}


#' Given a data set with all pairwise comparisons including A-B and B-A
#' comparisons, take the larger of the similarity score
#'
#' @param allResults needs to have columns `compare`, `newImage`, and a column
#'   with similarity scores (see next parameter)
#' @param similarityCol name of column with similarity scores, e.g. `"corr"`
#'
#' @return a data frame with one row per comparison, adding a column for the B-A
#'   match and a column `corrMax`
#'
#' @export

removeDups <- function(allResults, similarityCol) {
    out <- dplyr::inner_join(allResults[, c("compare", "newImage", similarityCol, "match")], allResults[, c("compare", "newImage", similarityCol)], by = c("compare" = "newImage", "newImage" = "compare")) # corr.x and corr.y

    out$corrMax <- NA
    for (i in 1:nrow(out)) {
        tmp <- out[i, c("compare", "newImage")]
        tmp <- sort(tmp)
        out[i, c("compare", "newImage")] <- tmp
        out$corrMax[i] <- max(out[i, c(paste0(similarityCol, ".x"), paste0(similarityCol, ".y"))])
    }
    out <- out[duplicated(out[, c("compare", "newImage", "corrMax")]) == FALSE, ]
    return(out)
}

#' Hierarchical clustering to generate final clusters
#'
#' Given similarities/predictions for each pair, generate transitive closures
#' using hierarchical clustering, returning pairwise match prediction
#'
#' @param pairs pairs with `predCol` column, which are similarity
#'   scores/predictions (higher = more similar)
#' @param predCol name of column with similarities/predictions, input as
#'   character, e.g. "myPreds"
#' @param myCutoff similarity cutoff (>= myCutoff corresponds to a "match")
#' @param myMethod type of linkage, can be "single", "complete", "average",
#'   "minimax"
#' @param hash1 name of column with first item in comparison, e.g. "image1"
#' @param hash2 name of column with second item in comparison, e.g. "image2"
#' @return vector (logical) indicating if pair is matched or not, length is the
#'   same as `pairs`
#' @export

linksAnalysis <- function(pairs, predCol, myCutoff, myMethod, hash1, hash2) {
    tmp <- which(pairs[, predCol] >= myCutoff)
    if (length(tmp) == 0) {
        return(rep(0, nrow(pairs)))
    } else {
        forHierarchical <- cbind(pairs[tmp, c(hash1, hash2)], preds = pairs[tmp, predCol])
        # cat("Original number of links: ", nrow(forHierarchical), "\n")
        hashes <- unique(c(forHierarchical[, hash1], forHierarchical[, hash2]))
        hashes <- sort(hashes)

        distMat <- matrix(NA, nrow = length(hashes), ncol = length(hashes))
        distMat[lower.tri(distMat)] <- 1

        for (ii in 1:nrow(forHierarchical)) {
            tmpi <- which(hashes == forHierarchical[ii, hash1])
            tmpj <- which(hashes == forHierarchical[ii, hash2])
            i <- max(tmpi, tmpj)
            j <- min(tmpi, tmpj)
            distMat[i, j] <- 1 - forHierarchical$preds[ii]
        }

        distObj <- as.dist(distMat, diag = FALSE, upper = FALSE)

        if (myMethod != "minimax") {
            hcluster <- hclust(distObj, method = myMethod)
            clustersAll <- cutree(hcluster, h = 1 - myCutoff)
        } else if (myMethod == "minimax") {
            hcluster <- protoclust::protoclust(distObj)
            clustersAll <- protoclust::protocut(hcluster, h = 1 - myCutoff)$cl
        }

        out <- data.frame(clustersAll) %>% dplyr::group_by(clustersAll) %>% dplyr::summarize(clusterSize = length(clustersAll))
        # cat("Cluster Sizes: ")
        # print(table(out$clusterSize))

        adjacencyMatrix <- GraphAT::clust2Mat(clustersAll)
        tmp <- sum(adjacencyMatrix)/2
        # cat("Final number of links: ", tmp, "\n")
        # cat("Links added: ", tmp - nrow(forHierarchical), "\n")

        ################ need to get it to return the link col
        adjacencyMatrix[lower.tri(adjacencyMatrix)] <- 0
        getHashes <- data.frame(which(adjacencyMatrix == 1, arr.ind = TRUE))

        getHashes[, hash1] <- hashes[getHashes[, "row"]]
        getHashes[, hash2] <- hashes[getHashes[, "col"]]

        tmp <- dplyr::left_join(pairs[, c(hash1, hash2)], getHashes[, c("row", hash1, hash2)], by = c(hash1, hash2))
        tmp$link <- 0
        tmp$link[!is.na(tmp$row)] <- 1
        return(tmp$link)
    }
}


