
rupclass_model <- function (X_train,
                            class_train,
                            X_test,
                            model_name,
                            alpha_train,
                            alpha_test,
                            # The proportion of obs in the Xtest to be trimmed
                            # The proportion of obs in the Xtraining to be trimmed
                            restr_factor,
                            ctrl_EM,
                            ctrl_restr,
                            ctrl_init,
                            verbose,
                            ...)
{
  logLik_raedda <-
    function(X_train, ltrain, X_test, fitm) {
      #for computing the ll within the EM algorithm
      mTau.train <- matrix(log(fitm$parameters$pro), nrow(X_train),
                           fitm$G, byrow = TRUE)
      lDensity.train <-
        tryCatch(
          do.call(mclust::cdens, c(
            list(data = X_train,
                 logarithm = TRUE), fitm
          )),
          error = function(e)
            - Inf
        )
      sum.train <-
        tryCatch(
          mTau.train + lDensity.train,
          error = function(e)
            - Inf
        )
      mat.train <- ltrain * sum.train
      ll.train <- sum(mat.train)
      ll.test <-
        tryCatch(
          sum(do.call(mclust::dens, c(
            list(data = X_test, logarithm = TRUE),
            fitm
          ))),
          error = function(e)
            - Inf,
          warning = function(w)
            - Inf
        )
      ll <- ll.train + ll.test
      ll <-
        ifelse(is.nan(ll) |
                 is.na(ll), -Inf, ll) #If llis NA or Nan it proceeds till the next estep and then the loop breaks
      ll
    }

  if(ncol(X_train)!=ncol(X_test)){
    stop("Training and Test set must have the same number of variables")
  }
  N_train <- nrow(X_train)
  N_test <- nrow(X_test)
  ltrain <- mclust::unmap(class_train)
  d <- ncol(X_train)
  G <-
    dim(ltrain)[2] # the number of groups in the training set
  if (G == 1 & !(model_name %in% c("EII", "EEI", "EEE"))) {
    warning('With G=1, only EII, EEI and EEE models are to be considered, others are redundant \n')
  }
  class_train <- as.factor(class_train)
  classLabel <- levels(class_train)

# Initialization ----------------------------------------------------------


  if (alpha_train != 0) {

    nsamp <- ctrl_init$n_samp
    max_iter_init <- ctrl_init$max_iter

    N_train_trim <- N_train - ceiling(N_train * alpha_train)
    fitm <-
      fitm_TRAIN <-
      patterned_MCD(nsamp = nsamp,
                    # the starting values for the EM algorithm are obtained from the estimation starting with nsamp J_g subsets of sample size (d+1), inspired to what done in \cite{Hubert2018}
                    X_train,
                    class_train,
                    ltrain,
                    alpha_train,
                    model_name,
                    max_iter_init)$fitm
  } else {
    N_train_trim <- NULL
    fitm <- fitm_TRAIN <-
      mclust::mstep(modelName = model_name,
                    data = X_train,
                    z = ltrain) #the starting values for the EM algorithm are the parameters obtained using only the training set
    fitm$X <-
      X_train # this is used in the initialization of raedda_transductive_model: for VVE and EVE models for which the constraint is not satisfied the data from which the M-step is computed are needed for the MM algorithm

  }
  #the robust initialization gives an idea of the magnitude of the eigenvalue-ratio in the known groups, this information can be exploited for avoiding setting in advance the value of restr.factor
  eigenvalues_known_groups <- tryCatch(apply(fitm$parameters$variance$sigma, 3, function(s2) eigen(s2, only.values = T)$val), error=function(e) 1)
  restr_factor_train <- abs(max(eigenvalues_known_groups)/min(eigenvalues_known_groups)) #in order to avoid spurious solution, I can set restr.factor to be no larger than the eigenvalue-ratio in the known groups. Abs is added for avoiding numerical problems

  if(alpha_test!=0){
    N_test_trim <- N_test-ceiling(N_test * alpha_test)
  } else {
    N_test_trim <- NULL
  }
  Xtrain_fit <- X_train
  Xtest_fit <- X_test
  ltrain_fit <- ltrain

if(is.null(restr_factor)){
  restr_factor <- restr_factor_train
}
#   else if (restr.factor<restr.factor_train){
#   warning(paste("The selected value for restr.factor =", restr.factor,"is smaller than the differences among group scatters in the known groups. You may want to increase it."))
# }
  if (!any(is.na(fitm$parameters))) {
    # if there are NA or NULL something went wrong and I will not compute the constrained maximization
    suppressWarnings(fitm <-
      constr_Sigma(
        fitm = fitm,
        restr_factor = restr_factor,
        extra_groups = 1:fitm$G, # it identifies the transductive approach
        ctrl_restr = ctrl_restr
      )) # this performs constrained estimation of Sigma according to the selected model, that is, the initial values in the EM algorithm satisfy the eigenvalues-ratio
  }

# EM algorithm ------------------------------------------------------------

  EM_tol <- ctrl_EM$tol
  EM_max_iter <- ctrl_EM$max_iter
  aitken <- ctrl_EM$aitken

  llold <- -Inf
  ll <- -Inf
  llstore <- rep(-Inf, 3) # for aitken
  criterion <- TRUE
  iter <- 0

  while (criterion) {
    iter <- iter + 1
    fite <-
      tryCatch(do.call(mclust::estep, c(list(data = X_test), fitm)), error = function(e) {
        list(z = NA)
      }) #expectation step, it remains the same
    emptyz <- TRUE
    if (all(!is.na(fite$z)))
    {
      emptyz <- FALSE
      z <- fite$z
      z_fit <- fite$z

      # Concentration Step Xtest
      if (alpha_test != 0){
        D <-
          do.call(mclust::dens, c(list(
            data = X_test, #compute the Density for Parameterized MVN Mixtures for each obs
            logarithm = T
          ), fitm)) #I temporarily discard those alpha_test% of obs whose density is lowest
        pos_trimmed_test <-
          which(D <= (sort(D, decreasing = F)[[ceiling(N_test * alpha_test)]]))
        z_fit <- z[-pos_trimmed_test, ,drop=F]
        Xtest_fit <- X_test[-pos_trimmed_test, ,drop=F]
      }

      # Concentration Step Xtrain
      if(alpha_train!=0){
        D_Xtrain_cond <-
          do.call(mclust::cdens, c(list(
            data = X_train, # computing the component density, this is done because I am interested in trimming also obs that might have
            logarithm = T
          ), fitm)) #been WRONGLY assinged to a class
        ind_D_Xtrain_cdens <- cbind(1:N_train,mclust::map(ltrain)) # matrix with obs and respective group from ltrain
        D_Xtrain <-
          D_Xtrain_cond[ind_D_Xtrain_cdens] # I compute D_g conditioning ON the fact that I know the supposed true class
        pos_trimmed_train <-
          which(D_Xtrain <= (sort(D_Xtrain, decreasing = F)
                             [[ceiling(N_train * alpha_train)]]))
        ltrain_fit <- ltrain[-pos_trimmed_train, ,drop=F]
        Xtrain_fit <- X_train[-pos_trimmed_train, , drop=F]
      }

      Xall <- rbind(Xtrain_fit, Xtest_fit)
      zall <- rbind(ltrain_fit, z_fit)
      fitm <-
        mclust::mstep(modelName = model_name,
              data = Xall,
              z = zall)
      if (!any(is.na(fitm$parameters$variance$sigma))) {
        # if there are NA something went wrong and I will not compute the constrained maximization
        if(fitm$modelName=="VVE"|fitm$modelName=="EVE"){
          fitm$X <- Xall # I add the data on which the M-step is computed since I need them for the MM
        }
        suppressWarnings(fitm <-
          constr_Sigma(
            fitm = fitm,
            restr_factor = restr_factor,
            extra_groups = 1:fitm$G, # it identifies the transductive approach
            ctrl_restr = ctrl_restr
          ) )# this performs constrained estimation of Sigma according to the selected model
      }

      ll <-
        logLik_raedda(
          X_train = Xtrain_fit,
          ltrain = ltrain_fit,
          X_test = Xtest_fit,
          fitm = fitm
        )

      llstore <- c(llstore[-1], ll)
      if (aitken) {
        criterion <- (Aitken(llstore)$linf - ll) > EM_tol
      } else {
        criterion <- (ll - llold) > EM_tol
      }
      criterion <- (criterion) & (iter < EM_max_iter)
      llold <- ll
    }
    else {
      criterion <- FALSE
    }
  }

# Results Collection ------------------------------------------------------


  res <- list()
  if (!emptyz  & ll != -Inf) {
    res$N_train <- N_train
    res$N_train_after_trimming <- N_train_trim
    res$N_test <- N_test
    res$N_test_after_trimming <- N_test_trim
    res$d <- d
    res$G <- G
    res$iter <- iter
    res$converged <- (iter < EM_max_iter)
    res$model_name <- model_name
    res$restr_factor <- restr_factor
    res$parameters <- fitm$parameters

    # eigen_final <- apply(res$parameters$variance$sigma, 3, function(x) eigen(x, only.values = T)$val)
    #
    # if(max(eigen_final)/min(eigen_final)<=restr.factor){
    #   warning(paste("The result is artificially constrained due to restr.factor= ", restr.factor), call. = FALSE)
    # }

    #Name the different groups
    names(res$parameters$pro) <- classLabel
    if(d==1){
      names(res$parameters$mean) <- classLabel
    } else {
    colnames(res$parameters$mean) <- classLabel
    }

    fitetrain <- do.call(mclust::estep, c(list(data = X_train),
                                    fitm))
    ztrain <- fitetrain$z
    cltrain <-
      factor(sapply(mclust::map(ztrain), function(i)
        classLabel[i]), levels = classLabel) # I classify a posteriori also the trimmed units
    # if alpha_train>0
    pos_trimmed_train = NULL
    cltrain_after_trimming = NULL
    if (alpha_train != 0) {
      D_Xtrain_cond <-
        do.call(mclust::cdens, c(list(
          data = X_train, # computing the component density, this is done because I am interested in trimming also obs that might have
          logarithm = T
        ), fitm)) #been WRONGLY assinged to a class
      ind_D_Xtrain_cdens <- cbind(1:N_train, mclust::map(ltrain)) # matrix with obs and respective group from ltrain
      D_Xtrain <- D_Xtrain_cond[ind_D_Xtrain_cdens] # I compute D_g conditioning ON the fact that I know the supposed true class
      pos_trimmed_train <-
        which(D_Xtrain <= (sort(D_Xtrain, decreasing = F)
                           [[ceiling(N_train * alpha_train)]]))
      cltrain_after_trimming <- cltrain
      cltrain_after_trimming <-
        factor(cltrain, levels = c(classLabel, "0"))
      cltrain_after_trimming[pos_trimmed_train] <- "0"
    }
    ##
    res$train <- list()
    res$train$z <- ztrain
    res$train$cl <- cltrain
    res$train$cl_after_trimming <- cltrain_after_trimming
    res$train$obs_trimmed <- pos_trimmed_train
    res$train$alpha_train <- alpha_train
    cl <-
      factor(sapply(map(z), function(i)
        classLabel[i]), levels = classLabel)
    pos_trimmed_test = NULL
    cltest_after_trimming = NULL
    if (alpha_test != 0) {
      D <-
        do.call(mclust::dens, c(list(
          data = X_test, #compute the Density for Parameterized MVN Mixtures for each obs
          logarithm = T
        ), fitm))
      pos_trimmed_test <-
        which(D <= (sort(D, decreasing = F)[[ceiling(N_test * alpha_test)]]))
      cltest_after_trimming <-
        factor(cl, levels = c(classLabel, "0"))
      cltest_after_trimming[pos_trimmed_test] <- "0"
    }
    res$test <- list()
    res$test$z <- z
    res$test$cl <- cl
    res$test$cl_after_trimming <- cltest_after_trimming
    res$test$obs_trimmed <- pos_trimmed_test
    res$test$alpha_test <- alpha_test
    bic.all <-
      robust_bic_raedda_t(
        modelName = model_name,
        loglik = ll,
        n = fitm$n,
        d = fitm$d,
        G = fitm$G,
        restr_factor = restr_factor
      )
    res$ll <- ll
    res$bic <- bic.all
  } else {
    res$error <- "Groups too small for this model type"
    res$ll <- NA
    res$bic <- NA
  }
  res
}
