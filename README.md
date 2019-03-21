
<!-- README.md is generated from README.Rmd. Please edit that file -->
cartridges3D
============

`cartridges` adapted for 3D topographies.

Installation
------------

You can install `cartridges3D` from github with:

``` r
# install.packages("devtools")
devtools::install_github("xhtai/cartridges3D")
```

Usage
-----

This is an example of how to use the package. First we pre-process the data, selecting breechface marks using RANSAC, leveling the image, removing circular symmetry and filtering the image. Next in the alignment step we need to find rotation and translation parameters -- this can be done using a grid search or using the Lucas-Kanade algorithm. The former is more accurate but slower.

The following is some example code for running this on multiple images.

Pre-processing:

``` r
library(cartridges3D)
fileList <- system("ls /media/xtai/4AF561E807105059/3D/data/Fadul/cc/*.x3p", intern = TRUE)

basis701 <- getBasisFunctions(701) # just run this once
for (i in 1:length(fileList)) {
    cat(i, ", ")
    processed <- allPreprocess3D(fileList[i])
    outName <- sub(".x3p", ".Rdata", fileList[i])
    outName <- sub("/cc/", "/processed_9-22/", outName)
    save(processed, file = outName)
}
```

Comparisons:

``` r
fileList <- system("ls /media/xtai/4AF561E807105059/3D/data/Fadul/processed_9-22/*.Rdata", intern = TRUE)

allPairwise <- function(imageName) {
    load(imageName)
    newImage <- processed
    out <- data.frame(compare = fileList, corr = NA, dx = NA, dy = NA, theta = NA, stringsAsFactors = FALSE)
    index <- which(out$compare == imageName)
    out <- out[-index, ]
    for (i in 1:nrow(out)) {
        #cat(i, ", ")
        load(out$compare[i])
        out[i, 2:5] <- calculateCCFmaxSearch(processed, newImage) # new image gets rotated
    }
    return(out)
}

Sys.time()
for (j in 1:length(fileList)) {
    cat(j, ", ")
    compare <- allPairwise(fileList[j])
    outName <- sub("/processed_9-22/", "/results_9-24/", fileList[j])
    save(compare, file = outName)
}
Sys.time()
```

Consolidating results and removing duplicates:

``` r
load("./allMetadata_3-17-2019.Rdata")

allResults <- c()
for (i in 1:(length(fileList))) {
    # if (i %% 10000 == 0) cat(i, ", ")
    load(fileList[i])
    compare$newImage <- fileList[i]
    allResults <- rbind(allResults, compare)
}
allResults$compare <- unlist(lapply(strsplit(allResults$compare, split = "/"), FUN = function(x) x[3]))
allResults <- dplyr::left_join(allResults, metadata[, c("rdataName", "GunNumber")], by = c("compare" = "rdataName"))
names(allResults)[length(names(allResults))] <- "compareGun"
allResults <- dplyr::left_join(allResults, metadata[, c("rdataName", "GunNumber")], by = c("newImage" = "rdataName"))
names(allResults)[length(names(allResults))] <- "newGun"

allResults$match <- as.numeric(allResults$compareGun == allResults$newGun)

removedDups <- removeDups(allResults, "corr")
```

Calculate precision and recall:

``` r
fg <- removedDups$corrMax[removedDups$match == 1]
bg <- removedDups$corrMax[removedDups$match == 0]
# PR Curve
pr <- PRROC::pr.curve(scores.class0 = fg, scores.class1 = bg, curve = TRUE)
pr <- pr$x
names(pr) <- c("recall", "precision", "cutoff")
```

Hierarchical clustering:

``` r
cutoffs <- seq(from = 0, to = 1, by = .01) 
for (i in 1:length(cutoffs)) { 
    removedDups[paste0("single", cutoffs[i])] <- linksAnalysis(removedDups, "corrMax", cutoffs[i], "single", "compare", "newImage")
    removedDups[paste0("complete", cutoffs[i])] <- linksAnalysis(removedDups, "corrMax", cutoffs[i], "complete", "compare", "newImage")
    removedDups[paste0("average", cutoffs[i])] <- linksAnalysis(removedDups, "corrMax", cutoffs[i], "average", "compare", "newImage")
    removedDups[paste0("minimax", cutoffs[i])] <- linksAnalysis(removedDups, "corrMax", cutoffs[i], "minimax", "compare", "newImage")
}

########### precision-recall
precRecall <- data.frame(link = rep(c("single", "complete", "average", "minimax"), each = length(cutoffs)), cutoff = rep(cutoffs, 4), precision = NA, recall = NA, stringsAsFactors = FALSE)

for (i in 1:nrow(precRecall)) {
    numerator <- sum(removedDups[paste0(precRecall$link[i], precRecall$cutoff[i])] >= .5 & removedDups$match == 1) # preds are 0 or 1 so doesn't matter that i used .5
    denom <- sum(removedDups[paste0(paste0(precRecall$link[i], precRecall$cutoff[i]))] >= .5)
    if (denom == 0) {
        precRecall$precision[i] <- 1
    } else {
        precRecall$precision[i] <- numerator/denom
    }
    precRecall$recall[i] <- numerator/sum(removedDups$match == 1)
}
```

Results
-------

These are example results using the Fadul data set.

![](allResults_distributionCCFmax_9-24.png)

License
-------

The `cartridges3D` package is licensed under GPLv3 (<http://www.gnu.org/licenses/gpl.html>).
