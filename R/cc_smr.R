cc_smr <- function(smr_params, v_e, alpha, seed=NULL){

  # CHANGE THIS SO IT RETURNS ONE CC.

  # first get the critical constants for the smr procedure
  idx = which(!is.na(smr_params$p)) # which of the contrasts are product contrasts,
  # inferred from the smr_params list.
  ccs[idx] <- mapply(smr_crit, alpha=alpha[idx],
                     p=smr_params$p[idx],
                     q=smr_params$q[idx],
                     n=rep(v_e, times=length(idx)),#smr_params$nu1[idx],
                     MoreArgs=list(n_sim=100000, seed=seed))
  # now take the square root of the critical constants, multiplied by nu1
  # see Bird & Hadzi-Pavlovic (2005), page 360,
  # for the conversion of SMR into a critical constant
  ccs[idx] <- sqrt(ccs[idx] * smr_params$nu1[idx])
  ccs
}
