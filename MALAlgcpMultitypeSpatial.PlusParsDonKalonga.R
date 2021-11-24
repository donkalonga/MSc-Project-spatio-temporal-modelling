MALAlgcpMultitypeSpatial.PlusParsDonKalonga <- function (mcmcloop, inits, adaptivescheme, M, N, Mext, Next, 
          mcens, ncens, formulaList, zml, Zmat, model.priorsList, model.initsList, 
          fftgrid, spatial.covmodelList, nis, cellarea, spatialvals, 
          cellInside, MCMCdiag, gradtrunc, gridfun, gridav, marks, 
          ntypes, d) 
{
  SpatialOnlyMode <- TRUE
  SpatialPlusParameters <- TRUE
  SpatioTemporalPlusParameters <- FALSE
  MultiTypeMode <- TRUE
  cellInsideLogical <- as.logical(cellInside)
  cidx <- ntypes + 1
  nlevs <- length(formulaList)
  M <- M
  N <- N
  temporal.fitted <- Inf
  GFinitialise(gridfun)
  GAinitialise(gridav)
  nsamp <- 0
  icount <- 0
  MCMCacc <- 0
  y.mean <- array(0, dim = c(M, N, ntypes + 1))
  y.var <- array(0, dim = c(M, N, ntypes + 1))
  EY.mean <- array(0, dim = c(M, N, ntypes + 1))
  EY.var <- array(0, dim = c(M, N, ntypes + 1))
  etainvtrans <- list()
  etaval <- list()
  etainv <- list()
  cp <- list()
  for (i in 1:(ntypes + 1)) {
    etainvtrans[[i]] <- model.priorsList[[i]]$etaprior$inverse_transform
    if (is.null(model.initsList)) {
      etaval[[i]] <- model.priorsList[[i]]$etaprior$mean
    }
    else {
      etaval[[i]] <- model.initsList[[i]]$etainit
    }
    etainv[[i]] <- etainvtrans[[i]](etaval[[i]])
    cp[[i]] <- CovParameters(list(sigma = etainv[[i]][1], 
                                  phi = etainv[[i]][2]))
  }
  mlev <- sapply(formulaList, function(x) {
    variablesinformula(x)[1]
  })
  SigmaBeta <- list()
  betavals <- list()
  dfr <- attr(Zmat, "data.frame")
  modls <- list()
  off <- c()
  rm(off)
  for (i in 1:ntypes) {
    dfr[[mlev[i]]] <- nis[[i]][cellInsideLogical]
    dfr$off <- log(cellarea * spatialvals[[i]][cellInsideLogical])
    dfr$off[is.infinite(dfr$off)] <- NA
    mod <- glm(formulaList[[i]], data = dfr, family = quasipoisson, 
               offset = off)
    modls[[i]] <- mod
    betavals[[i]] <- BetaParameters(coefficients(mod))
    SigmaBeta[[i]] <- vcov(mod)
  }
  Zlist <- list()
  Ztlist <- list()
  tm <- matrix(FALSE, Mext, Next)
  tm[1:M, 1:N] <- TRUE
  for (i in 1:ntypes) {
    Zlist[[i]] <- matrix(0, Next * Mext, ncol = ncol(zml[[i]]))
    Zlist[[i]][as.vector(tm), ] <- zml[[i]]
    Ztlist[[i]] <- t(Zlist[[i]])
  }
  neta <- list()
  nbeta <- list()
  for (i in 1:(ntypes + 1)) {
    neta[[i]] <- length(etaval[[i]])
  }
  for (i in 1:ntypes) {
    nbeta[[i]] <- length(betavals[[i]])
  }
  gammaVar <- list()
  CB <- list()
  etaCovMat <- list()
  GPls <- list()
  gmm <- 0
  for (k in 1:ntypes) {
    y <- log(nis[[k]]/(cellarea * spatialvals[[k]] * exp(matrix(as.vector(Zlist[[k]] %*% 
                                                                            betavals[[k]]), Mext, Next))))
    y[is.na(y) | is.infinite(y)] <- -cp[[k]]$sigma^2/2
    GPdummy <- GPrealisation(gamma = matrix(0, Mext, Next), 
                             fftgrid = fftgrid, covFunction = spatial.covmodelList[[k]], 
                             covParameters = cp[[k]], d = d)
    gammainit <- GammafromY(Y = y, rootQeigs = GPdummy$rootQeigs, 
                            mu = cp[[k]]$mu)
    GP <- GPrealisation(gamma = gammainit, fftgrid = fftgrid, 
                        covFunction = spatial.covmodelList[[k]], covParameters = cp[[k]], 
                        d = d)
    gmm <- gmm + gammainit
    GPls[[k]] <- GP
    eMat <- cellarea * spatialvals[[k]] * exp(matrix(as.vector(Zlist[[k]] %*% 
                                                                 betavals[[k]]), nrow(GP$gamma), ncol(GP$gamma))) * 
      GP$expY
    fi <- (1/length(GP$Y)) * Re(fft(GP$invrootQeigs * fft((1/length(GP$Y)) * 
                                                            Re(fft(GP$invrootQeigs * fft(eMat, inverse = TRUE))), 
                                                          inverse = TRUE)))
    gammaVar[[k]] <- 2 * (1.65^2/((Mext * Next)^(1/3))) * 
      (-1)/(-1 - fi)
    EZ <- as.vector(cellarea * spatialvals[[k]] * exp(matrix(as.vector(Zlist[[k]] %*% 
                                                                         betavals[[k]]), nrow(GP$gamma), ncol(GP$gamma)) + 
                                                        GP$Y))
    CB[[k]] <- matrix(NA, nbeta[[k]], nbeta[[k]])
    for (i in 1:nbeta[[k]]) {
      for (j in 1:nbeta[[k]]) {
        CB[[k]][i, j] <- sum(Zlist[[k]][, i] * Zlist[[k]][, 
                                                          j] * EZ) + model.priorsList[[k]]$betaprior$precision[i, 
                                                                                                               j]
      }
    }
  }
  GPls[[ntypes + 1]] <- GPrealisation(gamma = matrix(0, Mext, 
                                                     Next), fftgrid = fftgrid, covFunction = spatial.covmodelList[[ntypes + 
                                                                                                                     1]], covParameters = cp[[ntypes + 1]], d = d)
  if (TRUE) {
    tmpls <- list()
    additional_scaleconst <- (0.234/0.574)
    for (k in 1:(ntypes + 1)) {
      lensq <- 10
      sqsigma <- seq(model.priorsList[[k]]$etaprior$mean[1] - 
                       2 * sqrt(model.priorsList[[k]]$etaprior$variance[1, 
                                                                        1]), model.priorsList[[k]]$etaprior$mean[1] + 
                       2 * sqrt(model.priorsList[[k]]$etaprior$variance[1, 
                                                                        1]), length.out = lensq)
      sqphi <- seq(model.priorsList[[k]]$etaprior$mean[2] - 
                     2 * sqrt(model.priorsList[[k]]$etaprior$variance[2, 
                                                                      2]), model.priorsList[[k]]$etaprior$mean[2] + 
                     2 * sqrt(model.priorsList[[k]]$etaprior$variance[2, 
                                                                      2]), length.out = lensq)
      ltarmat <- matrix(NA, lensq, lensq)
      for (i in 1:lensq) {
        for (j in 1:lensq) {
          cpA <- CovParameters(list(sigma = exp(sqsigma[i]), 
                                    phi = exp(sqphi[j])))
          gpA <- GPls
          etatmp <- etaval
          etatmp[[k]] <- c(sqsigma[i], sqphi[j])
          gpA[[k]] <- GPrealisation(gamma = GPls[[k]]$gamma, 
                                    fftgrid = fftgrid, covFunction = spatial.covmodelList[[k]], 
                                    covParameters = cpA, d = d)
          matent <- try(target.and.grad.MultitypespatialPlusPars(GPlist = gpA, 
                                                                 priorlist = model.priorsList, Zlist = Zlist, 
                                                                 Ztlist = Ztlist, eta = etatmp, beta = betavals, 
                                                                 nis = nis, cellarea = fftgrid$cellarea, spatial = spatialvals, 
                                                                 gradtrunc = gradtrunc)$logtarget)
          if (class(matent) != "try-error") {
            ltarmat[i, j] <- matent
          }
        }
      }
      ltarmat[is.infinite(ltarmat)] <- NA
      dffit <- data.frame(ltar = as.vector(ltarmat))
      exgr <- expand.grid(sqsigma, sqphi)
      dffit$sigma <- exgr[, 1]
      dffit$sigma2 <- exgr[, 1]^2
      dffit$phi <- exgr[, 2]
      dffit$phi2 <- exgr[, 2]^2
      dffit$sigmaphi <- exgr[, 1] * exgr[, 2]
      try(tarmod <- lm(ltar ~ sigma2 + sigma + phi2 + phi + 
                         sigmaphi, data = dffit))
      try(coef <- coefficients(tarmod))
      try(etaCovMat[[k]] <- additional_scaleconst * solve((-1) * 
                                                            matrix(c(2 * coef[2], coef[6], coef[6], 2 * coef[4]), 
                                                                   2, 2)))
    }
  }
  gammaVar[[cidx]] <- 0
  for (i in 1:ntypes) {
    gammaVar[[cidx]] <- gammaVar[[cidx]] + (1/ntypes) * gammaVar[[i]]
  }
  gammaVar[[cidx]] <- gammaVar[[cidx]]/ntypes
  rootgammaVar <- lapply(gammaVar, sqrt)
  sigma_eta <- list()
  SIGMA_ETA <- list()
  Q_eta <- list()
  sigma_beta <- list()
  Q_beta <- list()
  sigma_eta_chol <- list()
  sigma_beta_chol <- list()
  for (i in 1:(ntypes + 1)) {
    sigma_eta[[i]] <- (2.38^2/length(etaval[[i]])) * etaCovMat[[i]]
    SIGMA_ETA[[i]] <- sigma_eta[[i]]
    Q_eta[[i]] <- solve(sigma_eta[[i]])
    sigma_eta_chol[[i]] <- t(RandomFieldsUtils::chol(sigma_eta[[i]]))
    print(sigma_eta[[i]])
  }
  for (i in 1:ntypes) {
    sigma_beta[[i]] <- (1.65^2/(length(betavals[[i]])^(1/3))) * 
      solve(CB[[i]])
    Q_beta[[i]] <- as.matrix(solve(sigma_beta[[i]]))
    sigma_beta_chol[[i]] <- as.matrix(t(RandomFieldsUtils::chol(sigma_beta[[i]])))
    print(sigma_beta[[i]])
  }
  betarec <- c()
  etarec <- c()
  adapt_h <- adaptivescheme
  h <- initialiseAMCMC(adapt_h)
  gammacount <- etacount <- betacount <- 1
  GPlist <- list()
  for (i in 1:(ntypes + 1)) {
    GPlist[[i]] <- GPrealisation(gamma = matrix(0, Mext, 
                                                Next), fftgrid = fftgrid, covFunction = spatial.covmodelList[[i]], 
                                 covParameters = cp[[i]], d = d)
  }
  oldtags <- target.and.grad.MultitypespatialPlusPars(GPlist = GPlist, 
                                                      priorlist = model.priorsList, Zlist = Zlist, Ztlist = Ztlist, 
                                                      eta = etaval, beta = betavals, nis = nis, cellarea = fftgrid$cellarea, 
                                                      spatial = spatialvals, gradtrunc = gradtrunc)
  trigger1 <- TRUE
  tarrec <- oldtags$logtarget
  hallrec <- h
  reject_its <- c()
  while (nextStep(mcmcloop)) {
    propmeans_gamma <- list()
    propmeans_eta <- list()
    propeta <- list()
    propetainv <- list()
    propcp <- list()
    propGPlist <- list()
    propmeans_beta <- list()
    propbeta <- list()
    for (i in (ntypes + 1):1) {
      propmeans_gamma[[i]] <- GPlist[[i]]$gamma + (h/2) * 
        gammaVar[[i]] * oldtags$gradgamma[[i]]
      propmeans_eta[[i]] <- etaval[[i]]
      propeta[[i]] <- as.vector(propmeans_eta[[i]] + sqrt(h) * 
                                  sigma_eta_chol[[i]] %*% rnorm(neta[[i]]))
      propetainv[[i]] <- etainvtrans[[i]](propeta[[i]])
      propcp[[i]] <- CovParameters(list(sigma = propetainv[[i]][1], 
                                        phi = propetainv[[i]][2]))
      propGPlist[[i]] <- GPrealisation(gamma = propmeans_gamma[[i]] + 
                                         sqrt(h) * rootgammaVar[[i]] * rnorm(Mext * Next), 
                                       fftgrid = fftgrid, covFunction = spatial.covmodelList[[i]], 
                                       covParameters = propcp[[i]], d = d)
    }
    for (i in 1:ntypes) {
      propmeans_beta[[i]] <- betavals[[i]] + (h/2) * sigma_beta[[i]] %*% 
        oldtags$gradbeta[[i]]
      propbeta[[i]] <- BetaParameters(as.vector(propmeans_beta[[i]] + 
                                                  sqrt(h) * sigma_beta_chol[[i]] %*% rnorm(nbeta[[i]])))
    }
    proptags <- target.and.grad.MultitypespatialPlusPars(GPlist = propGPlist, 
                                                         priorlist = model.priorsList, Zlist = Zlist, Ztlist = Ztlist, 
                                                         eta = propeta, beta = propbeta, nis = nis, cellarea = fftgrid$cellarea, 
                                                         spatial = spatialvals, gradtrunc = gradtrunc)
    revpropmeans_gamma <- list()
    revpropmeans_eta <- list()
    revpropmeans_beta <- list()
    gamma_acceptance_contrib <- 0
    eta_acceptance_contrib <- 0
    beta_acceptance_contrib <- 0
    for (i in (ntypes + 1):1) {
      revpropmeans_gamma[[i]] <- propGPlist[[i]]$gamma + 
        (h/2) * gammaVar[[i]] * proptags$gradgamma[[i]]
      revpropmeans_eta[[i]] <- propeta[[i]]
      gamma_acceptance_contrib <- gamma_acceptance_contrib - 
        sum((GPlist[[i]]$gamma - revpropmeans_gamma[[i]])^2/(2 * 
                                                               h * gammaVar[[i]])) + sum((propGPlist[[i]]$gamma - 
                                                                                            propmeans_gamma[[i]])^2/(2 * h * gammaVar[[i]]))
    }
    for (i in 1:ntypes) {
      revpropmeans_beta[[i]] <- propbeta[[i]] + (h/2) * 
        sigma_beta[[i]] %*% proptags$gradbeta[[i]]
      beta_acceptance_contrib <- beta_acceptance_contrib + 
        (-(0.5/h) * t(betavals[[i]] - revpropmeans_beta[[i]]) %*% 
           Q_beta[[i]] %*% (betavals[[i]] - revpropmeans_beta[[i]])) - 
        (-(0.5/h) * t(propbeta[[i]] - propmeans_beta[[i]]) %*% 
           Q_beta[[i]] %*% (propbeta[[i]] - propmeans_beta[[i]]))
    }
    ac <- exp(proptags$logtarget - oldtags$logtarget + gamma_acceptance_contrib + 
                eta_acceptance_contrib + beta_acceptance_contrib)
    ac <- min(ac, 1)
    icount <- icount + 1
    MCMCacc <- ((icount - 1)/icount) * MCMCacc + ac/icount
    if (proptags$logtarget == -Inf | is.na(ac) | is.nan(ac)) {
      warning("One possible cause of this warning is that there may be evidence in the data for quite large values of the spatial correlation parameter. If this is the case, then this warning means that the MCMC chain has wandered into a region of the phi-space that causes the variance matrix of Y (computed by the DFT) to become non positive-definite. One possible way of rectifying this issue is to restart the chain using a larger value of 'ext' in the call to lgcpPredictMultitypeSpatialPlusPars. You should do this if the warning message is repeated many times. If this warning message appears at all then you should be warned that inference from this run may not be reliable: the proposed move at this iteration was rejected.", 
              immediate. = TRUE)
      ac <- 0
      reject_its <- c(reject_its, iteration(mcmcloop))
    }
    if (ac > runif(1)) {
      GPlist <- propGPlist
      oldtags <- proptags
      etaval <- propeta
      betavals <- propbeta
    }
    if (iteration(mcmcloop)%%100 == 0) {
      print(h)
    }
    h <- updateAMCMC(adapt_h)
    if (is.retain(mcmcloop)) {
      hallrec <- c(hallrec, h)
      tarrec <- c(tarrec, oldtags$logtarget)
      betarec <- rbind(betarec, unlist(betavals))
      etarec <- rbind(etarec, unlist(etaval))
      Y <- GPlist2array(GPlist = GPlist, element = "Y")
      expY <- GPlist2array(GPlist = GPlist, element = "expY")
      nsamp <- nsamp + 1
      y.mean <- ((nsamp - 1)/nsamp) * y.mean + Y[1:M, 1:N, 
      ]/nsamp
      EY.mean <- ((nsamp - 1)/nsamp) * EY.mean + expY[1:M, 
                                                      1:N, ]/nsamp
      if (nsamp > 1) {
        y.var <- ((nsamp - 2)/(nsamp - 1)) * y.var + 
          (nsamp/(nsamp - 1)^2) * (y.mean - Y[1:M, 1:N, 
          ])^2
        EY.var <- ((nsamp - 2)/(nsamp - 1)) * EY.var + 
          (nsamp/(nsamp - 1)^2) * (EY.mean - expY[1:M, 
                                                  1:N, ])^2
      }
      GFupdate(gridfun)
      GAupdate(gridav)
    }
  }
  retlist <- list(lasth = rev(hallrec)[1], lastGAM = oldtags$Gamma)
  retlist$mcmcacc <- MCMCacc
  retlist$y.mean <- lgcpgrid(y.mean, xvals = mcens[1:M], yvals = ncens[1:N])
  retlist$y.var <- lgcpgrid(y.var, xvals = mcens[1:M], yvals = ncens[1:N])
  retlist$EY.mean <- lgcpgrid(EY.mean, xvals = mcens[1:M], 
                              yvals = ncens[1:N])
  retlist$EY.var <- lgcpgrid(EY.var, xvals = mcens[1:M], yvals = ncens[1:N])
  retlist$gridfunction <- GFreturnvalue(gridfun)
  retlist$gridaverage <- GAreturnvalue(gridav)
  retlist$mcmcinfo <- mcmcloop
  retlist$gradtrunc <- gradtrunc
  retlist$etarec <- etarec
  retlist$betarec <- betarec
  retlist$glmfit <- modls
  retlist$Zlist <- Zlist
  retlist$hallrec <- hallrec
  retlist$tarrec <- tarrec
  retlist$cinslogical <- cellInsideLogical
  retlist$reject_its <- reject_its
  return(retlist)
}