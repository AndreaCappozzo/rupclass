raedda_l_model <- function(X_train,
                           class_train,
                           alpha_train,
                           # The proportion of obs in the Xtrain to be trimmed
                           model_name,
                           ctrl_init,
                           ...) {

  N_train <- nrow(X_train)
  ltrain <- mclust::unmap(class_train)
  d <- ncol(X_train)
  G <- ncol(ltrain)
  class_train <- as.factor(class_train)
  classLabel <- levels(class_train)

  # Core algorithm ----------------------------------------------------------


  if (alpha_train != 0) {

    nsamp <- ctrl_init$n_samp
    max_iter_init <- ctrl_init$max_iter

    N_train_trim <- N_train - ceiling(N_train * alpha_train)
    robust_result <- patterned_MCD(nsamp = nsamp, # I perform the estimation starting from nsamp J_g subsets of sample size (d+1), inspired to what done in \cite{Hubert2018}
                                   X_train,
                                   class_train,
                                   ltrain,
                                   alpha_train,
                                   model_name, max_iter_init)
    ll <- robust_result$ll
    fitm <- robust_result$fitm

  } else {
    N_train_trim <- NULL
    fitm <- mclust::mstep(modelName = model_name,
                          data = X_train,
                          z = ltrain)
    mTau.train <-
      matrix(log(fitm$parameters$pro),
             nrow(X_train),
             fitm$G,
             byrow = TRUE)
    lDensity.train <- do.call(mclust::cdens, c(list(
      data = X_train,
      logarithm = TRUE
    ), fitm))

    sum.train <- mTau.train + lDensity.train
    mat.train <- ltrain * sum.train
    ll <- sum(mat.train)
  }


  # Checking if errors in the procedure -------------------------------------

  fitetrain <-
    tryCatch(
      do.call(mclust::estep, c(list(data = X_train), fitm)),
      error = function(e) {
        list(z = NA)
      }
    )
  emptyz <- ifelse(all(!is.na(fitetrain$z)), yes = FALSE, no = TRUE)

  # Results Collection ------------------------------------------------------

  if (!emptyz) {
    res <- list()
    res$N_train <- N_train
    res$N_train_after_trimming <- N_train_trim
    res$alpha_train <- alpha_train
    res$d <- d
    res$G <- G
    res$model_name <- model_name
    res$parameters <- fitm$parameters
    #Name the different groups
    names(res$parameters$pro) <- classLabel
    if (d == 1) {
      names(res$parameters$mean) <- classLabel
    } else {
      colnames(res$parameters$mean) <- classLabel
    }

    ztrain <- fitetrain$z
    cltrain <-
      factor(sapply(map(ztrain), function(i)
        classLabel[i]), levels = classLabel) # I classify a posteriori also the trimmed units
    pos_trimmed_train <- NULL
    cltrain_after_trimming <- NULL
    if (alpha_train != 0) {
      D_Xtrain_cond <- do.call(mclust::cdens, c(list(
        data = X_train, # computing the component density, this is done because I am interested in trimming also obs that might have
        logarithm = T
      ), fitm)) # been WRONGLY assinged to a class
      ind_D_Xtrain_cdens <-
        cbind(1:N_train, mclust::map(ltrain)) # matrix with obs and respective group from ltrain
      D_Xtrain <-
        D_Xtrain_cond[ind_D_Xtrain_cdens] # I compute D_g conditioning ON the fact that I know the supposed true class
      pos_trimmed_train <-
        # I trim the conditional density \phi(x_n; \mu_g, \Sigma_g) when x_n comes from group g
        which(D_Xtrain <= (sort(D_Xtrain, decreasing = F)
                           [[ceiling(N_train * alpha_train)]]))
      cltrain_after_trimming <- cltrain
      cltrain_after_trimming <-
        factor(cltrain, levels = c(classLabel, "0"))
      cltrain_after_trimming[pos_trimmed_train] <- "0"
    }

    res$train <- list()
    res$train$z <- ztrain
    res$train$cl <- cltrain
    res$train$cl_after_trimming <- cltrain_after_trimming
    res$train$obs_trimmed <- pos_trimmed_train
    res$train$alpha_train <- alpha_train
    bic_all <-
      2 * ll - mclust::nMclustParams(
        modelName = model_name,
        d = d,
        G = G,
        noise = FALSE,
        equalPro = FALSE
      ) * log(fitm$n) # if alpha_train !=0 this is trimmed BIC
    res$ll <- ll
    res$bic <- bic_all
  }
  else {
    res <- list()
    res$error <-
      "Either training groups too small or trimming level too large for this model type"
    res$ll <- NA
    res$bic <- NA
  }
  res
}
