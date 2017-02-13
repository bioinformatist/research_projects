PC <-
  # Function for PCA by eigen or svd method. Still has some problems.
  function(X,
           method = 'svd',
           scaled = FALSE,
           center = FALSE,
           graph = FALSE,
           rm.na = TRUE,
           print.results = TRUE) {
    if (any(is.na(X))) {
      tmp <- X
      if (rm.na == TRUE) {
        X <- na.omit(data.frame(X))
        X <- as.matrix(X)
      }
      else{
        X[is.na(X)] = matrix(
          apply(X, 2, mean, na.rm = TRUE),
          ncol = ncol(X),
          nrow = nrow(X),
          byrow = TRUE
        )[is.na(X)]
      }
    }
    else{
      tmp <- X
    }
    if (method == 'eigen') {
      if (scaled == 1) {
        X1 = cor(X)
        X2 = scale(X)
      }
      else{
        X1 = cov(X)
        X2 = scale(X, scale = F)
      }
      total.var <- sum(diag(cov(X2)))
      values <-
        eigen(X1)$values
      vectors <- eigen(X1)$vectors
      sdev = sqrt(values)
    }
    if (method == 'svd') {
      if (sum(scaled, center) > 1) {
        X2 <- scale(X)
      } else{
        if (scaled == 1) {
          X2 = scale(X, center = F)
        } else{
          if (center == 1) {
            X2 = scale(X, scale = F)
          } else{
            X2 = X
          }
        }
      }
      total.var <- sum(diag(cov(X2)))
      var <- nrow(X2) - 1
      vectors <- svd(X2)$v
      sdev = svd(X2)$d / sqrt(var)
      values <- sdev * sdev
      prop.var <- rep(NA, ncol(X))
      cum.var <- rep(NA, ncol(X))
      scores <- X2 %*% vectors
      namex <- as.character(1:ncol(X))
      scorenames <- rep(NA, ncol(X))
      for (i in 1:(ncol(X))) {
        scorenames[i] <- do.call(paste, c("PC", as.list(namex[i]), sep = ""))
      }
      colnames(scores) <- scorenames
      rownames(vectors) <- colnames(X)
      colnames(vectors) <- scorenames
      for (i in 1:ncol(X)) {
        prop.var[i] <- var(scores[, i]) / total.var
      }
      for (i in 1:ncol(X)) {
        cum.var[i] <- sum(prop.var[1:i])
      }
      importance <- t(matrix(c(sdev, prop.var, cum.var), ncol = 3))
      importance <- as.table(importance)
      colnames(importance) <- scorenames
      rownames(importance) <-
        c("Standard Deviation",
          "Proportion of Variance",
          "Cumulative Proportion")
      z <-
        list(
          values = values,
          vectors = vectors,
          scores = scores,
          importance = importance
          ,
          sdev = sdev
        )
      if (graph == 1) {
        biplot(
          scores[, 1:2],
          vectors[, 1:2],
          main = "Biplot of Data",
          xlab = do.call
          (paste, c("PC1 (", as.list(
            round(z$importance[2, 1] * 100, 2)
          ), "%)", sep = ""))
          ,
          ylab = do.call(paste, c("PC2
                                  (", as.list(
                                    round(z$importance[2, 2] * 100, 2)
                                  ), "%)", sep = "")),
          cex = 0.7
        )
        windows()
        screeplot(z, type = 'l', main = 'Screeplot of Components')
        abline(1, 0, col = 'red', lty = 2)
      }
      if (print.results == TRUE) {
        if (method == "eigen") {
          print("PCA Analysis Using Spectral Decomposition")
        }
        if (method == "svd") {
          print("PCA Analyis Using Singular Value Decomposition")
        }
        if (any(is.na(tmp))) {
          if (rm.na == TRUE) {
            print("Warning:One or more rows of data were omitted from
                  analysis")
          }
          if (rm.na == FALSE) {
            print("Warning: Mean of the variable was used for Missing
                  values")
          }
        }
        print(importance)
      }
      z <-
        list(
          values = values,
          vectors = vectors,
          scores = scores,
          importance = importance,
          sdev = sdev
        )
    }
  }