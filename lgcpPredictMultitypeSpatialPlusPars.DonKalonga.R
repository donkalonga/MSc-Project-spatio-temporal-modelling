lgcpPredictMultitypeSpatialPlusPars.DonKalonga <- function (formulaList, sd, typemark = NULL, Zmat = NULL, model.priorsList, 
          model.initsList = NULL, spatial.covmodelList, cellwidth = NULL, 
          poisson.offset = NULL, mcmc.control, output.control = setoutput(), 
          gradtrunc = Inf, ext = 2, inclusion = "touching") 
{
  spatial.offsetList <- poisson.offset
  regionalcovariates <- NULL
  pixelcovariates <- NULL
  gridsize <- NULL
  formula <- aggregateformulaList(formulaList)
  starttime <- Sys.time()
  if (is.null(typemark)) {
    if (!is.null(dim(sd$marks))) {
      stop("Since sd$marks is not a vector, you must specify which of the marks defines the point type by setting the argument typemark.")
    }
  }
  else {
    if (!any(names(sd$marks) == typemark)) {
      stop(paste("None of the marks is named", typemark))
    }
  }
  if (is.null(dim(sd$marks))) {
    if (!is.factor(sd$marks)) {
      warning(paste("Converting", typemark, "to a factor"), 
              .immediate = TRUE)
    }
    marks <- as.factor(sd$marks)
  }
  else {
    if (!is.factor(sd$marks$typemark)) {
      warning(paste("Converting", typemark, "to a factor"), 
              .immediate = TRUE)
    }
    marks <- as.factor(sd$marks$typemark)
  }
  ntypes <- length(formulaList)
  if (sd$window$type == "rectangle") {
    sd$window <- as.polygonal(sd$window)
  }
  if (is.null(cellwidth) & is.null(gridsize)) {
    stop("Either cell width OR grid size must be specified")
  }
  if (!is.null(cellwidth) & !is.null(gridsize)) {
    stop("Either cell width OR grid size must be specified")
  }
  if (!all(sapply(gridsize, is.pow2))) {
    stop("All elements of gridsize must be a power of 2")
  }
  if (!is.null(gridsize)) {
    approxcw <- diff(sd$window$xrange)/gridsize[1]
    cwseq <- seq(approxcw/2, 2 * approxcw, length.out = 500)
    cwfun <- function(cw) {
      ow <- selectObsWindow(sd, cw)
      return(c(ow$M, ow$N))
    }
    gsmat <- t(sapply(cwseq, cwfun))
    tf <- apply(gsmat, 1, function(x) {
      return(all(x == gridsize))
    })
    if (sum(tf) == 0) {
      stop("No sensible observation window found: either change gridsize, or specify cellwidth instead")
    }
    else {
      cellwidth <- cwseq[min(which(tf))]
    }
  }
  if (!is.null(gradtrunc)) {
    if (gradtrunc < 0) {
      stop("gradtrunc must be non-negative")
    }
  }
  if (mcmc.control$burnin > mcmc.control$mala.length) {
    stop("Number of burnin iterations must be less than the total number of iterations")
  }
  ow <- selectObsWindow(sd, cellwidth)
  sd <- ow$xyt
  M <- ow$M
  N <- ow$N
  if (M * N >= (256^2)) {
    Sys.sleep(1)
    cat("\n")
    warning("USING LARGE FFT GRID: COMPUTATION MAY BE SLOW ON SOME MACHINES ...", 
            .immediate = TRUE)
    cat("\n")
  }
  cat(paste("FFT Grid size: [", ext * M, " , ", 
            ext * N, "]\n", sep = ""))
  Sys.sleep(1)
  rm(ow)
  spatial <- list()
  if (is.null(spatial.offsetList)) {
    for (i in 1:ntypes) {
      spatial[[i]] <- spatialAtRisk(list(X = seq(sd$window$xrange[1], 
                                                 sd$window$xrange[2], length.out = M), Y = seq(sd$window$yrange[1], 
                                                                                               sd$window$yrange[2], length.out = N), Zm = matrix(1, 
                                                                                                                                                 M, N)))
    }
  }
  else if (class(spatial.offsetList) == "SpatialAtRisk") {
    for (i in 1:ntypes) {
      spatial[[i]] <- spatial.offsetList
    }
  }
  else if (class(spatial.offsetList) == "list") {
    for (i in 1:ntypes) {
      if (!any(class(spatial.offsetList[[i]]) == "spatialAtRisk")) {
        spatial[[i]] <- spatialAtRisk(spatial.offsetList[[i]])
      }
      else {
        spatial[[i]] <- spatial.offsetList[[i]]
      }
    }
  }
  else {
    stop("Invalid spatial.offsetList")
  }
  funk <- function(spatial) {
    if (any(class(spatial) == "fromXYZ")) {
      spatial$Zm <- spatial$Zm * attr(spatial, "NC")
    }
    if (any(class(spatial) == "fromSPDF")) {
      spatial$atrisk <- spatial$atrisk * attr(spatial, 
                                              "NC")
      spatial$spdf$atrisk <- spatial$atrisk
    }
    return(spatial)
  }
  spatial <- lapply(spatial, funk)
  study.region <- sd$window
  if (!is.null(attr(Zmat, "gridobj"))) {
    gridobj <- attr(Zmat, "gridobj")
  }
  else {
    gridobj <- genFFTgrid(study.region = study.region, M = M, 
                          N = N, ext = ext, inclusion = inclusion)
  }
  del1 <- gridobj$del1
  del2 <- gridobj$del2
  Mext <- gridobj$Mext
  Next <- gridobj$Next
  mcens <- gridobj$mcens
  ncens <- gridobj$ncens
  cellarea <- gridobj$cellarea
  cellInside <- gridobj$cellInside
  x <- gridobj$mcens
  y <- gridobj$ncens
  xidx <- rep(1:Mext, Next)
  yidx <- rep(1:Next, each = Mext)
  dxidx <- pmin(abs(xidx - xidx[1]), Mext - abs(xidx - xidx[1]))
  dyidx <- pmin(abs(yidx - yidx[1]), Next - abs(yidx - yidx[1]))
  d <- sqrt(((x[2] - x[1]) * dxidx)^2 + ((y[2] - y[1]) * dyidx)^2)
  if (is.null(Zmat)) {
    Zmat <- cov.interp.fft(formula = formula, W = study.region, 
                           regionalcovariates = regionalcovariates, pixelcovariates = pixelcovariates, 
                           mcens = mcens[1:M], ncens = ncens[1:N], cellInside = cellInside[1:M, 
                                                                                           1:N])
  }
  else {
    if (!isTRUE(all.equal(mcens[1:M], attr(Zmat, "mcens"))) | 
        !isTRUE(all.equal(ncens[1:N], attr(Zmat, "ncens")))) {
      stop("FFT grid and Zmat are on different grids. Please recompute Zmat using 'getZmat'.")
    }
  }
  zml <- getZmats(Zmat = Zmat, formulaList = formulaList)
  spatialvals <- list()
  for (i in 1:ntypes) {
    spatialvals[[i]] <- fftinterpolate(spatial[[i]], mcens, 
                                       ncens, ext = ext)
    spatialvals[[i]] <- spatialvals[[i]] * cellInside
  }
  mLoop = mcmcLoop(N = mcmc.control$mala.length, burnin = mcmc.control$burnin, 
                   thin = mcmc.control$retain, progressor = mcmcProgressTextBar)
  nsamp <- floor((mLoop$N - mLoop$burnin)/mLoop$thin)
  if (!is.null(output.control$gridfunction) & class(output.control$gridfunction)[1] == 
      "dump2dir") {
    cat("WARNING: disk space required for saving is approximately ", 
        round(nsamp * object.size(array(runif((M) * (N) * 
                                                (ntypes + 1)), dim = c((M), (N), ntypes + 1)))/1024^2, 
              2), " Mb, ", sep = "")
    if (!output.control$gridfunction$forceSave) {
      m <- menu(c("yes", "no"), title = "continue?")
      if (m == 1) {
        cat("Note: to bypass this menu, set forceSave=TRUE in dump2dir\n")
        Sys.sleep(2)
      }
      else {
        stop("Stopped")
      }
    }
  }
  mlev <- sapply(formulaList, function(x) {
    variablesinformula(x)[1]
  })
  nis <- list()
  ct1 <- 0
  ct2 <- 0
  for (i in 1:ntypes) {
    nis[[i]] <- getCounts(xyt = sd[marks == mlev[i], ], M = M, 
                          N = N, ext = ext)
    ct1 <- ct1 + sum(nis[[i]])
    nis[[i]] <- nis[[i]] * cellInside
    ct2 <- ct2 + sum(nis[[i]])
  }
  if (ct2 < ct1) {
    warning(paste(ct1 - ct2, " data points lost due to discretisation.", 
                  sep = ""), immediate. = TRUE)
  }
  gridfun <- output.control$gridfunction
  if (is.null(gridfun)) {
    gridfun <- nullFunction()
  }
  gridav <- output.control$gridmeans
  if (is.null(gridav)) {
    gridav <- nullAverage()
  }
  lg <- MALAlgcpMultitypeSpatial.PlusParsDonKalonga(mcmcloop = mLoop, 
                                          inits = mcmc.control$inits, adaptivescheme = mcmc.control$adaptivescheme, 
                                          M = M, N = N, Mext = Mext, Next = Next, mcens = mcens, 
                                          ncens = ncens, formulaList = formulaList, zml = zml, 
                                          Zmat = Zmat, model.priorsList = model.priorsList, model.initsList = model.initsList, 
                                          fftgrid = gridobj, spatial.covmodelList = spatial.covmodelList, 
                                          nis = nis, cellarea = cellarea, spatialvals = spatialvals, 
                                          cellInside = cellInside, MCMCdiag = mcmc.control$MCMCdiag, 
                                          gradtrunc = gradtrunc, gridfun = gridfun, gridav = gridav, 
                                          marks = marks, ntypes = ntypes, d = d)
  endtime <- Sys.time()
  timetaken <- endtime - starttime
  lg$xyt <- sd
  lg$M <- M
  lg$N <- N
  lg$aggtimes <- NA
  lg$tdiffs <- NA
  lg$vars <- NA
  lg$spatial <- spatial
  lg$temporal <- NA
  lg$grid <- gridobj
  lg$nis <- nis
  lg$mcens <- mcens[1:M]
  lg$ncens <- ncens[1:N]
  lg$mcmcpars <- mcmc.control
  lg$timetaken <- timetaken
  lg$spatialonly <- TRUE
  lg$spatialonlyplusparameters <- TRUE
  lg$ext <- ext
  lg$cellInside <- cellInside[1:M, 1:N]
  lg$inclusion <- inclusion
  lg$poisson.offset <- spatialvals
  lg$priors <- model.priorsList
  lg$covFct <- spatial.covmodelList
  class(lg) <- c("lgcpPredictMultitypeSpatialPlusParameters", 
                 "lgcpPredict", "lgcpobject")
  return(lg)
}