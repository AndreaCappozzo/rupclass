redda <-
  function(X_train,
           class_train,
           alpha_train = 0,
           model_names = NULL,
           ctrl_init = control_init_rupclass(),
           verbose = interactive(),
           ...) {
    X_train <- data.matrix(X_train)
    if (is.null(model_names)) {
      if(ncol(X_train)==1){
        model_names <- c("E", "V")
      } else {
        model_names <- mclust::mclust.options("emModelNames")
      }
    }
    if (verbose) {
      cat("fitting ...\n")
      flush.console()
      pbarL <- txtProgressBar(min = 0,
                              max = length(model_names),
                              style = 3)
      on.exit(close(pbarL))
      ipbarL <- 0
    }
    RES <- list()
    bestBIC <- -Inf
    RES[["Best"]] <- list()
    for (model_name in model_names) {
      RES[[model_name]] <- list()
      RES[[model_name]] <- raedda_l_model(X_train,
                                           class_train,
                                           alpha_train,
                                           model_name,
                                           ctrl_init = ctrl_init,
                                           ...)
      if (!is.na(RES[[model_name]]$bic)) {
        if (RES[[model_name]]$bic > bestBIC) {
          RES[["Best"]] <- RES[[model_name]]
          bestBIC <- RES[[model_name]]$bic
        }
      }
      if (verbose) {
        ipbarL <- ipbarL + 1
        setTxtProgressBar(pbarL, ipbarL)
      }
    }
    RES$Best$train$X_train <-
      X_train #I also return input training data
    if (verbose) {
      cat("\nA",
          RES$Best$model_name,
          "patterned model was selected")
    }
    class(RES) <- "redda"
    RES
  }
