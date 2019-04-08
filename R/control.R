
# Utility control functions -----------------------------------------------

control_init_rupclass <- function(n_samp = 50,
                         max_iter = 1e02) {
  # Initialization control paramaters
  list(
    n_samp = n_samp,
    max_iter = max_iter
  )
}

control_EM_rupclass <- function(tol = 1e-05,
                    max_iter = 1e02,
                    aitken = TRUE)
  # EM control parameters
{
  list(tol = tol,
       max_iter = max_iter,
       aitken = aitken)
}

control_restr_rupclass <- function(tol = 1e-10, max_iter = 1e02, MM_tol = 1e-10, MM_max_iter = 1e04)
  # Constrained maximization control parameters
  #MM is related to Majorization - Minimization algorithm (used for VVE, EVE models)
{
  list(
    tol = tol,
    max_iter = max_iter,
    MM_tol = MM_tol,
    MM_max_iter = MM_max_iter
  )
}
