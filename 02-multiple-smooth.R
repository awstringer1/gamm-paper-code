## Simulations for "Inference for generalized additive mixed models via penalized marginal likelihood"
## Multiple smooth terms
## Alex Stringer
## January 2025
## File 02: multiple smooth functions

runinparallel <- .Platform$OS.type == "unix" # Run in parallel?

seed <- 447113 # Arbitrary

args <- commandArgs(TRUE) # Returns character(0) if interactive
if (length(args) > 0) {
  numsims <- args[1]
  simname <- args[2]
  seed <- seed * as.numeric(args[3]) # Use a different random seed for each version
} else {
  numsims <- 2
  simname <- "sims-multiple-20250228-v1"
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
ff1 <- function(x) {
  n <- length(x)
  out <- numeric(n)
  for (i in 1:n) {
    out[i] <- sin(2 * pi * x[i])
  }
  out
}
ff2 <- function(x) {
  n <- length(x)
  out <- numeric(n)
  for (i in 1:n) {
    out[i] <- cos(2 * pi * x[i])
  }
  out
}

simulate_data <- function(m, n, sigma) {

  if (length(n) == 1) n <- rep(n, m)
  N <- sum(n)

  # generate the random effects
  u <- stats::rnorm(m, 0, sigma)

  # covariate
  id <- Reduce(c, mapply(function(x, y) rep(x, each = y), 1:m, n, SIMPLIFY = FALSE))
  x1 <- runif(N)
  x2 <- runif(N)
  
  # Linear predictor
  fmla <- y ~ (1 | id)
  df <- data.frame(id = id, x1 = x1, x2 = x2, y = 0)
  reterms <- lme4::glFormula(fmla, data = df, family = stats::binomial)
  Z <- t(reterms$reTrms$Zt)
  eta <- as.numeric(Z %*% u + ff1(x1) + ff2(x2))

  # Response
  pp <- 1 / (1 + exp(-eta))
  df$y <- stats::rbinom(N, 1, pp)

  df
}


## Function to fit all three models ##
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

  S1 <- thegam$smooth[[1]]$S[[1]] # Scaled penalty matrix
  S2 <- thegam$smooth[[2]]$S[[1]] # Scaled penalty matrix
  S <- as.matrix(bdiag(lambda[1] * S1, lambda[2] * S2))
  penmat <- rbind(
    0, cbind(0, S)
  ) 

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
  
  interceptest <- opt$par[1]

  beta1idx <- 1 + 1:ncol(S1)
  beta2idx <- 1 + (1 + ncol(S1)):(ncol(S2) + ncol(S1))

  beta1est <- opt$par[beta1idx]
  beta2est <- opt$par[beta2idx]

  xx <- seq(0, 1, length.out = 1e03)

  Xpred1 <- mgcv::PredictMat(thegam$smooth[[1]], data.frame(x1 = xx))
  ffpred1 <- Xpred1 %*% beta1est
  ffpred1sd <- sqrt(diag(Xpred1 %*% Hinv[beta1idx, beta1idx] %*% t(Xpred1)))
  ffpred1lower <- ffpred1 - zquant * ffpred1sd
  ffpred1upper <- ffpred1 + zquant * ffpred1sd

  Xpred2 <- mgcv::PredictMat(thegam$smooth[[2]], data.frame(x2 = xx))
  ffpred2 <- Xpred2 %*% beta2est
  ffpred2sd <- sqrt(diag(Xpred2 %*% Hinv[beta2idx, beta2idx] %*% t(Xpred2)))
  ffpred2lower <- ffpred2 - zquant * ffpred2sd
  ffpred2upper <- ffpred2 + zquant * ffpred2sd

  # Return the results from the GAM too
  ffpred1gamobj <- predict(thegam, newdata = data.frame(x1 = xx, x2 = 0), se.fit = TRUE, type = "terms")
  ffpred1gam <- ffpred1gamobj$fit[ ,1]
  ffpred1sdgam <- ffpred1gamobj$se.fit[ ,1]
  ffpred1lowergam <- ffpred1gam - zquant * ffpred1sdgam
  ffpred1uppergam <- ffpred1gam + zquant * ffpred1sdgam

  ffpred2gamobj <- predict(thegam, newdata = data.frame(x1 = 0, x2 = xx), se.fit = TRUE, type = "terms")
  ffpred2gam <- ffpred2gamobj$fit[ ,2]
  ffpred2sdgam <- ffpred2gamobj$se.fit[ ,2]
  ffpred2lowergam <- ffpred2gam - zquant * ffpred2sdgam
  ffpred2uppergam <- ffpred2gam + zquant * ffpred2sdgam

  # Return the results from the GAMM too
  ffpredgamm1obj <- predict(thegamm$gam, newdata = data.frame(x1 = xx, x2 = 0), se.fit = TRUE, type = "terms")
  ffpred1gamm <- ffpredgamm1obj$fit[ ,1]
  ffpred1sdgamm <- ffpredgamm1obj$se.fit[ ,1]
  ffpred1lowergamm <- ffpred1gamm - zquant * ffpred1sdgamm
  ffpred1uppergamm <- ffpred1gamm + zquant * ffpred1sdgamm

  ffpredgamm2obj <- predict(thegamm$gam, newdata = data.frame(x1 = 0, x2 = xx), se.fit = TRUE, type = "terms")
  ffpred2gamm <- ffpredgamm2obj$fit[ ,2]
  ffpred2sdgamm <- ffpredgamm2obj$se.fit[ ,2]
  ffpred2lowergamm <- ffpred2gamm - zquant * ffpred2sdgamm
  ffpred2uppergamm <- ffpred2gamm + zquant * ffpred2sdgamm

  list(
    x = xx,
    sigmainterval = c("estimate" = sigmaest),
    sigmagamm = sqrt(summary(thegamm$mer)$varcor[1]$id[1]),
    ffpred1new = data.frame("lower" = ffpred1lower, "estimate" = ffpred1, "upper" = ffpred1upper),
    ffpred1gam = data.frame("lower" = ffpred1lowergam, "estimate" = ffpred1gam, "upper" = ffpred1uppergam),
    ffpred1gamm = data.frame("lower" = ffpred1lowergamm, "estimate" = ffpred1gamm, "upper" = ffpred1uppergamm),
    ffpred2new = data.frame("lower" = ffpred2lower, "estimate" = ffpred2, "upper" = ffpred2upper),
    ffpred2gam = data.frame("lower" = ffpred2lowergam, "estimate" = ffpred2gam, "upper" = ffpred2uppergam),
    ffpred2gamm = data.frame("lower" = ffpred2lowergamm, "estimate" = ffpred2gamm, "upper" = ffpred2uppergamm),
    lambdagam = thegam$sp, lambdagamm = thegamm$lambdaest,
    pmlopt = opt # Return the optimization object for diagnostics
  )
}

## Simulation study ##
# numsims <- 200 # Set via command line
m <- c(1000, 2000, 5000, 10000)
n <- c(3) # Average group size
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

options(mc.cores = parallel::detectCores())

RNGkind("L'Ecuyer-CMRG") # For reproducibility with parallel

dosim <- function(lst) {
  cat(simname, ": simulation ", lst$idx, " of ", length(simlist), ", m = ", lst$m, ", n = ", lst$n, ", k = ", lst$k, ", sigma = ", lst$sigma, "\n", sep = "")
  # Average sample size = n
  nvec <- sample(2:(2 * (lst$n - 1)), size = lst$m, replace = TRUE)
  df <- with(lst, simulate_data(m, nvec, sigma))
  res <- tryCatch(fitnewmodel(
    smoothformula = y ~ s(x1, bs = "bs", k = lst$knots) + s(x2, bs = "bs", k = lst$knots),
    glmmformula = y ~ (1 | id),
    data = df,
    k = lst$k
  ), error = function(e) e, warning = function(w) w)
  if (inherits(res, "condition")) {
    cat("ERROR: Simulation ", lst$idx, " of ", length(simlist), ", m = ", lst$m, ", n = ", lst$n, ", k = ", lst$k, "\n", sep = "")
    cat("Error message: ", res$message, "\n")
    return(NULL)
  }

  truef1 <- ff1(res$x)
  bias1new <- mean(res$ffpred1new[["estimate"]] - truef1)
  bias1gam <- mean(res$ffpred1gam[["estimate"]] - truef1)
  bias1gamm <- mean(res$ffpred1gamm[["estimate"]] - truef1)
  
  covr1new <- with(res, mean(ffpred1new[["lower"]] <= truef1 & truef1 <= ffpred1new[["upper"]]))
  covr1gam <- with(res, mean(ffpred1gam[["lower"]] <= truef1 & truef1 <= ffpred1gam[["upper"]]))
  covr1gamm <- with(res, mean(ffpred1gamm[["lower"]] <= truef1 & truef1 <= ffpred1gamm[["upper"]]))

  truef2 <- ff2(res$x)
  bias2new <- mean(res$ffpred2new[["estimate"]] - truef2)
  bias2gam <- mean(res$ffpred2gam[["estimate"]] - truef2)
  bias2gamm <- mean(res$ffpred2gamm[["estimate"]] - truef2)
  
  covr2new <- with(res, mean(ffpred2new[["lower"]] <= truef2 & truef2 <= ffpred2new[["upper"]]))
  covr2gam <- with(res, mean(ffpred2gam[["lower"]] <= truef2 & truef2 <= ffpred2gam[["upper"]]))
  covr2gamm <- with(res, mean(ffpred2gamm[["lower"]] <= truef2 & truef2 <= ffpred2gamm[["upper"]]))

  biassigmagamm <- res$sigmagamm - lst$sigma
  biassigmanew <- res$sigmainterval["estimate"] - lst$sigma
  
  ressum <- c(
    "bias1new" = bias1new, "bias1gam" = bias1gam, "bias1gamm" = bias1gamm, 
    "covr1new" = covr1new, "covr1gam" = covr1gam, "covr1gamm" = covr1gamm,
    "bias2new" = bias2new, "bias2gam" = bias2gam, "bias2gamm" = bias2gamm, 
    "covr2new" = covr2new, "covr2gam" = covr2gam, "covr2gamm" = covr2gamm,
    "biassigmanew" = biassigmanew, "biassigmagamm" = biassigmagamm
  )

  c(unlist(lst), ressum)
}

if (runinparallel) {
  options(mc.cores = parallel::detectCores())
  RNGkind("L'Ecuyer-CMRG") # For reproducibility with parallel
  mc.reset.stream() # Reproducbility in parallel
}
# Do the simulations
tm <- Sys.time()
if (runinparallel) {
  # simlist <- simlist[sample.int(length(simlist), replace = FALSE)] # Do in random order so the progress bar is roughly accurate
  # sims <- bettermc::mclapply(simlist, dosim, mc.retry = -3, mc.cores = 14, mc.force.fork = TRUE, mc.allow.fatal = TRUE, mc.allow.error = TRUE)
  # sims <- mclapply(simlist, dosim)
  plan(multisession)
  sims <- future_lapply(simlist, dosim, future.stdout = NA, future.seed = TRUE)
} else {
  sims <- lapply(simlist, dosim)
}
simtime <- as.numeric(difftime(Sys.time(), tm, units='secs'))
cat("Finished simulations, they took", simtime, "seconds.\n")
cat("Saving simulations...\n")
save(sims, file = file.path(resultspath, simresultsname))
cat("Saved simulations to file:",file.path(resultspath,simresultsname),"\n")

