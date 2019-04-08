
# ROBUST BIC related functions --------------------------------------------



n_raedda_t_params <- function(modelName, d, G, restr_factor) {
  # Transductive approach: gives the number of estimated parameters for parameterizations
  # of the Gaussian mixture model that are used in MCLUST considering
  # a penalty term that takes into account the higher model complexity
  # that a higher restr_factor entails.

  # For models that do not require a restriction, we keep the original MCLUST function
  if (modelName %in% c("E", "EII", "EEI", "EEE", "EEV")) {
    return(mclust::nMclustParams(
      modelName = modelName,
      d = d,
      G = G
    ))
  }
  #alpha= #param related to mixing proportion and mean vector
  #gamma= #param related to orthogonal rotation
  #delta= #param related to eigenvalues
  alpha <- G * d + G - 1
  if (modelName == "V") {
    delta <-
      (G - 1) * (1 - 1 / restr_factor) + 1 # Cerioli et al 2018 univariate case (d=1)
    return(alpha + delta)
  }
  NAME <- unlist(strsplit(modelName, split = ""))
  #Orientation
  if (NAME[3] == "I") {
    gamma <- 0
  } else if (NAME[3] == "E") {
    gamma <- d * (d - 1) / 2
  } else{
    gamma <- G * d * (d - 1) / 2
  }
  # Eigenvalues
  if (NAME[1] == "E" & NAME[2] == "V") {
    delta <- G * d - (G - 1)
  } else if (NAME[1] == "V" & NAME[2] == "E") {
    delta <- G + d - 1
  } else if (NAME[1] == "V" & NAME[2] == "V") {
    delta <- G * d
  } else if (NAME[1] == "V" & NAME[2] == "I") {
    delta <- G
  }
  return(alpha + gamma + (delta - 1) * (1 - 1 / restr_factor) + 1) # Cerioli et al 2018
}


# RAEDDA transductive
robust_bic_raedda_t <-
  function (modelName,
            loglik,
            n,
            d,
            G,
            restr_factor,
            #This function computes the Robust BIC as defined in Cerioli et al 2018, according to the specified modelName
            ...) {
    nparams <-
      n_raedda_t_params(
        modelName = modelName,
        d = d,
        G = G,
        restr_factor = restr_factor
      )
    2 * loglik - nparams * log(n)
  }
