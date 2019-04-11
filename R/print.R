# Printing Methods for the different estimation procedures ------------------------------

print.rupclass <- function(x, ...)
{
  cat("\n")
  txt <- paste(" ", "Robust Updating Classification Rule \n")
  sep <- paste0(rep("-", max(nchar(txt)) + 1),
                collapse = "")
  cat(sep, "\n")
  cat(txt)
  cat(sep, "\n")
  cat("\n")
  cat( paste0("  ", "Model = ", x$Best$model_name, "\n") )
  cat( paste0("  ", "Training trimming level = ", x$Best$train$alpha_train, "\n") )
  cat( paste0("  ", "Test trimming level = ", x$Best$test$alpha_test, "\n") )
  cat( paste0("  ", "Log-likelihood= ", round(x$Best$ll, 3), "\n") )
  cat( paste0("  Robust BIC= ", round(x$Best$bic, 3), "\n") )
}

print.redda <- function(x, ...)
{
  cat("\n")
  txt <- paste(" ", "Robust Eigenvalue Decomposition Discriminant Analysis \n")
  sep <- paste0(rep("-", max(nchar(txt)) + 1),
                collapse = "")
  cat(sep, "\n")
  cat(txt)
  cat(sep, "\n")
  cat("\n")
  cat( paste0("  ", "Model = ", x$Best$model_name, "\n") )
  cat( paste0("  ", "Training trimming level = ", x$Best$train$alpha_train, "\n") )
  cat( paste0("  ", "Log-likelihood= ", round(x$Best$ll, 3), "\n") )
  cat( paste0("  Robust BIC= ", round(x$Best$bic, 3), "\n") )
}
