## Simulations for "Inference for generalized additive mixed models via penalized marginal likelihood"
## Alex Stringer
## January 2025
## File 04: simulations with a more complex function

runinparallel <- .Platform$OS.type == "unix" # Run in parallel using "future"?

seed <- 553169 # Arbitrary
args <- commandArgs(TRUE) # Returns character(0) if interactive
if (length(args) > 0) {
  numsims <- args[1]
  simname <- args[2]
  seed <- seed * as.numeric(args[3]) # Use a different random seed for each version
} else {
  numsims <- 100
  simname <- "sims-complex-20250228-v1"
}

## Packages ##
pkgs <- c(
  "lme4",
  "mgcv",
  "parallel",
  "dplyr",
  "ggplot2",
  "tidyr",
  "remotes",
  "future.apply"
)
for (pkg in pkgs) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste0("Could not find package ", pkg, ", installing from CRAN.\n"))
    install.packages(pkg)
    require(pkg, character.only = TRUE, quietly = TRUE)
  }
}

# aghqmm package for PML implementation
remotes::install_github("awstringer1/aghqmm", force = TRUE)
# If you want to set up remotes, this tutorial is helpful:
# https://carpentries.github.io/sandpaper-docs/github-pat.html
library(aghqmm)

## Set paths ##
basepath <- getwd()
resultspath <- file.path(basepath, "results")
if (!dir.exists(resultspath)) dir.create(resultspath)
figurespath <- file.path(basepath, "figures")
if (!dir.exists(figurespath)) dir.create(figurespath)
simresultsname <- paste0(simname, ".RData")
simsprocessedname <- paste0(simname, ".csv")

cores <- future::availableCores() - 1
cat("System time: ", format(Sys.time(), "%H:%m:%S"), ", executing ", numsims, " simulations, file ", simresultsname, ", using ", cores, " cores.\n", sep = "")

## Generate data
ff <- function(x) {
  n <- length(x)
  out <- numeric(n)
  for (i in 1:n) {
    out[i] <- (1 / 10) * (6 * dbeta(x[i], 30, 17) + 4 * dbeta(x[i], 3, 11)) - 1
  }
  out
}
# curve(ff)

simulate_data <- function(m, n, sigma) {

  if (length(n) == 1) n <- rep(n, m)
  N <- sum(n)

  # generate the random effects
  u <- stats::rnorm(m, 0, sigma)

  # covariate
  id <- Reduce(c, mapply(function(x, y) rep(x, each = y), 1:m, n, SIMPLIFY = FALSE))
  x <- runif(N)

  # Linear predictor
  fmla <- y ~ (1 | id)
  df <- data.frame(id = id, x = x, y = 0)
  reterms <- lme4::glFormula(fmla, data = df, family = stats::binomial)
  Z <- t(reterms$reTrms$Zt)
  eta <- as.numeric(Z %*% u + ff(x))

  # Response
  pp <- 1 / (1 + exp(-eta))
  df$y <- stats::rbinom(N, 1, pp)

  df
}

## Function to fit the 3 models ##
fitnewmodel <- function(smoothformula, glmmformula, data, k = 5) {
  ## GAM ##
  thegam <- mgcv::bam(
    smoothformula,
    data = data,
    family = binomial(),
    method = "REML",
    discrete = FALSE, nthreads = 1
  )
  lambda <- thegam$sp # Scaled smoothing parameter
  betastart <- coef(thegam) # Includes intercept
  S <- thegam$smooth[[1]]$S[[1]] # Scaled penalty matrix
  penmat <- rbind(
    0, cbind(0, lambda * S)
  ) # Add row/column of zeroes for the intercept

  # prepare data
  response <- all.vars(glmmformula)[1]
  reterms <- lme4::glFormula(glmmformula,data=data,family=stats::binomial) # TODO: update this for other families
  idvar <- names(reterms$reTrms$cnms)[1]
  # check if response is 0/1 coded
  if(!all.equal(sort(unique(data[[response]])),c(0,1))) 
    stop(paste0("Response should be coded as numeric 0/1, yours is coded as",sort(unique(data[[response]])),". I don't know how to automatically convert it, sorry!"))
  
  # Spline design matrix
  X <- model.matrix(thegam)

  modeldata <- list(
    X = X,
    y = data[[response]],
    group = data[[idvar]],
    id = as.numeric(table(data[[idvar]]))
  )
  
  # Prepare the data for the likelihood functions
  yy <- with(modeldata,split(y,group))
  XX <- with(modeldata,lapply(split(X,group),matrix,ncol=ncol(X)))
  d <- length(Reduce(c,reterms$reTrms$cnms)) # NOTE: not really tested
  
  # Quadrature
  gg <- mvQuad::createNIGrid(d,'GHe',k)
  nn <- mvQuad::getNodes(gg)
  ww <- mvQuad::getWeights(gg)
  
  control <- aghqmm::aghqmm_control(onlynllgrad = TRUE)
  
  thetastart <- c(betastart, 0)
  optfun <- function(theta) aghqmm:::optimizegammscalar(theta, yy, XX, penmat, nn, ww, control)$nll
  optgrad <- function(theta) aghqmm:::optimizegammscalar(theta, yy, XX, penmat, nn, ww, control)$grad
  opt <- optim(thetastart, optfun, optgrad, method = "BFGS", hessian = TRUE)
  
  ## GAMM ##
  thegamm <- gamm4::gamm4(smoothformula, random = update(glmmformula, NULL ~.), data = data, family = binomial(), REML = TRUE)

  # Confidence intervals
  zquant <- stats::qnorm(.975)
  H <- opt$hessian
  Hinv <- solve(H)
  paramsd <- sqrt(diag(Hinv))

  s <- length(opt$par)
  sigmaest <- exp(-opt$par[s] / 2)
  
  betaidx <- 2:(s - 1) # No intercept
  betaest <- opt$par[betaidx]
  xx <- seq(0, 1, length.out = 1e03)
  Xpred <- mgcv::PredictMat(thegam$smooth[[1]], data.frame(x = xx))
  ffpred <- Xpred %*% betaest
  ffpredsd <- sqrt(diag(Xpred %*% Hinv[betaidx, betaidx] %*% t(Xpred)))
  ffpredlower <- ffpred - zquant * ffpredsd
  ffpredupper <- ffpred + zquant * ffpredsd
  
  # Return the results from the GAM too
  ffpredgamobj <- predict(thegam, newdata = data.frame(x = xx), se.fit = TRUE, type = "terms")
  ffpredgam <- ffpredgamobj$fit[ ,1]
  ffpredsdgam <- ffpredgamobj$se.fit[ ,1]
  ffpredlowergam <- ffpredgam - zquant * ffpredsdgam
  ffpreduppergam <- ffpredgam + zquant * ffpredsdgam

  # Return the results from the GAMM too
  ffpredgamm4obj <- predict(thegamm$gam, newdata = data.frame(x = xx), se.fit = TRUE, type = "terms")
  ffpredgamm4 <- ffpredgamm4obj$fit[ ,1]
  ffpredsdgamm4 <- ffpredgamm4obj$se.fit[ ,1]
  ffpredlowergamm4 <- ffpredgamm4 - zquant * ffpredsdgamm4
  ffpreduppergamm4 <- ffpredgamm4 + zquant * ffpredsdgamm4

  list(
    x = xx,
    sigmainterval = c("estimate" = sigmaest),
    sigmagamm = sqrt(summary(thegamm$mer)$varcor[1]$id[1]),
    ffprednew = data.frame("lower" = ffpredlower, "estimate" = ffpred, "upper" = ffpredupper),
    ffpredgam = data.frame("lower" = ffpredlowergam, "estimate" = ffpredgam, "upper" = ffpreduppergam),
    ffpredgamm = data.frame("lower" = ffpredlowergamm4, "estimate" = ffpredgamm4, "upper" = ffpreduppergamm4),
    lambdagam = thegam$sp,
    pmlopt = opt # Return the optimization object for diagnostics
  )
}

## Simulation study ##
# numsims <- 100 # Set via command line
m <- c(1000, 2000, 5000, 10000)
# n <- c(3, 9)
n <- c(3, 9) # Average group size
sigma <- c(1)
k <- c(5)
knots <- 20

# Create simulation objects
simstodoframe <- expand.grid(
  sim = 1:numsims,
  n = n,
  m = m,
  sigma = sigma,
  k = k,
  knots = knots
)
simlist <- list()
idx <- 1
for (i in 1:nrow(simstodoframe)) {
  lst <- list(
    idx = idx,
    n = simstodoframe[i, "n"],
    m = simstodoframe[i, "m"],
    sigma = simstodoframe[i, "sigma"],
    k = simstodoframe[i, "k"],
    knots = simstodoframe[i, "knots"]
  )

  simlist <- c(simlist, list(lst))
  idx <- idx + 1
}

dosim <- function(lst) {
  cat(simname, ": simulation ", lst$idx, " of ", length(simlist), ", m = ", lst$m, ", n = ", lst$n, ", k = ", lst$k, ", sigma = ", lst$sigma, "\n", sep = "")
  # Average sample size = n
  nvec <- sample(2:(2 * (lst$n - 1)), size = lst$m, replace = TRUE)
  df <- with(lst, simulate_data(m, nvec, sigma))
  res <- tryCatch(fitnewmodel(
    smoothformula = y ~ s(x, bs = "bs", k = lst$knots),
    glmmformula = y ~ (1 | id),
    data = df,
    k = lst$k
  ), error = function(e) e, warning = function(w) w)
  if (inherits(res, "condition")) {
    cat("ERROR: Simulation ", lst$idx, " of ", length(simlist), ", m = ", lst$m, ", n = ", lst$n, ", k = ", lst$k, "\n", sep = "")
    cat("Error message: ", res$message, "\n")
    return(NULL)
  }

  truef <- ff(res$x)
  biasnew <- mean(res$ffprednew[["estimate"]] - truef)
  biasgam <- mean(res$ffpredgam[["estimate"]] - truef)
  biasgamm <- mean(res$ffpredgamm[["estimate"]] - truef)
  
  covrnew <- with(res, mean(ffprednew[["lower"]] <= truef & truef <= ffprednew[["upper"]]))
  covrgam <- with(res, mean(ffpredgam[["lower"]] <= truef & truef <= ffpredgam[["upper"]]))
  covrgamm <- with(res, mean(ffpredgamm[["lower"]] <= truef & truef <= ffpredgamm[["upper"]]))

  biassigmagamm <- res$sigmagamm - lst$sigma
  biassigmanew <- res$sigmainterval["estimate"] - lst$sigma
  
  ressum <- c(
    "biasnew" = biasnew, "biasgam" = biasgam, "biasgamm" = biasgamm, "biassigmanew" = biassigmanew, "biassigmagamm" = biassigmagamm,
    "covrnew" = covrnew, "covrgam" = covrgam, "covrgamm" = covrgamm
  )

  c(unlist(lst), ressum)
}

## Simulation ##

if (runinparallel) {
  cores <- future::availableCores() - 1
  options(mc.cores = cores)
  RNGkind("L'Ecuyer-CMRG") # For reproducibility with parallel
  mc.reset.stream() # Reproducbility in parallel
}
# Do the simulations
tm <- Sys.time()
if (runinparallel) {
  plan(multisession) # gamm4 crashes when using forking, use multisession instead.
  sims <- future_lapply(simlist, dosim, future.seed = TRUE)
} else {
  sims <- lapply(simlist, dosim)
}
simtime <- as.numeric(difftime(Sys.time(), tm, units='secs'))
cat("Finished simulations, they took", simtime, "seconds.\n")
cat("Saving simulations...\n")
save(sims, file = file.path(resultspath, simresultsname))
cat("Saved simulations to file:",file.path(resultspath,simresultsname),"\n")
