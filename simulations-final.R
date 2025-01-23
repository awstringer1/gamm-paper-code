## Simulations for "Inference for generalized additive mixed models via penalized marginal likelihood"
## Alex Stringer
## January 2025

runinparallel <- .Platform$OS.type == "unix" # Run in parallel?

## Packages ##
pkgs <- c(
  "lme4",
  "mgcv",
  "parallel",
  "TMB",
  "dplyr",
  "ggplot2",
  "tidyr",
  "remotes"
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
# Must contain a subdirectory "code/tmb" with the "gammapprox.cpp" file
simname <- "sims-20241115-v12"
resultspath <- file.path(basepath, "results")
if (!dir.exists(resultspath)) dir.create(resultspath)
figurespath <- file.path(basepath, "figures")
if (!dir.exists(figurespath)) dir.create(figurespath)
simresultsname <- paste0(simname, ".RData")
simsprocessedname <- paste0(simname, ".csv")

## Load TMB code for the GAMM ##

tmbfile <- file.path(basepath, "code/tmb/gammapprox")

compile(paste0(tmbfile, ".cpp"))
dyn.load(dynlib(tmbfile))

prepare_data_gamm_tmb <- function(df, knots = 30) {
  # Prepare the data for input into TMB

  # Random effects design matrix
  fmla <- y ~ (1 | id)
  reterms <- lme4::glFormula(fmla, data = df, family = stats::binomial)
  Z <- t(reterms$reTrms$Zt)

  # Spline design matrix
  smooth <- mgcv::smoothCon(s(x, bs = "bs", k = knots), data = df, absorb.cons = TRUE)[[1]]
  X <- smooth$X
  S <- smooth$S[[1]]

  # Prediction
  xx <- seq(0, 1, length.out = 1e03)
  Xpred <- mgcv::PredictMat(smooth, data.frame(x = xx))

  list(
    y = df$y,
    X = X,
    Z = Z,
    S = S,
    xx = xx,
    Xpred = Xpred
  )
}

fit_gamm_tmb <- function(df, knots = 20) {
  tmbdat <- prepare_data_gamm_tmb(df, knots = knots)

  tmbparams <- list(
    logprec = 0,
    loglambda = 0,
    alpha = 0,
    beta = rep(0, ncol(tmbdat$X)),
    u = rep(0, ncol(tmbdat$Z))
  )

  template <- TMB::MakeADFun(
    tmbdat,
    tmbparams,
    random = c("alpha", "beta", "u"),
    DLL = "gammapprox",
    silent = TRUE
  )

  opt <- optim(c(0, 0), template$fn, template$gr, method = "BFGS")

  uidx <- which(names(template$env$last.par.best) == "u")
  betaidx <- which(names(template$env$last.par.best) == "beta")
  alphaidx <- which(names(template$env$last.par.best) == "alpha")
  uest <- template$env$last.par.best[uidx]
  betaest <- template$env$last.par.best[betaidx]
  alphaest <- template$env$last.par.best[alphaidx]

  ffpred <- as.numeric(tmbdat$Xpred %*% betaest)
  H <- template$env$spHess(template$env$last.par.best, random = TRUE)
  L <- Cholesky(H, LDL = FALSE)
  tmpzero <- new("dgCMatrix")
  tmpzero@Dim <- as.integer(c(nrow(tmbdat$Xpred), ncol(H)))
  tmpzero@p <- rep(0L, tmpzero@Dim[2]+1)
  tmpzero[ , betaidx] <- tmbdat$Xpred
  ffpredsd <- as.numeric(sqrt(diag(tmbdat$Xpred %*% solve(L, t(tmpzero))[betaidx, ])))

  list(
    ffpred = ffpred, 
    ffpredsd = ffpredsd, 
    sigmaest = exp(-opt$par[1] / 2), 
    lambdaest = exp(opt$par[2]), normsdcurve = sqrt(sum(ffpredsd^2))
  )
}


## Generate data
ff <- function(x) {
  n <- length(x)
  out <- numeric(n)
  for (i in 1:n) {
    # out[i] <- (1 / 10) * (6 * dbeta(x[i], 30, 17) + 4 * dbeta(x[i], 3, 11)) - 1
    out[i] <- sin(2 * pi * x[i])
  }
  out
}

simulate_data <- function(m, n, sigma) {

  N <- m * n

  # generate the random effects
  u <- stats::rnorm(m, 0, sigma)

  # covariate
  id <- rep(1:m, each = n)
  # x <- rnorm(N, .5, .5/3)
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

## Function to fit the penalized GLMM ##
fitnewmodel <- function(smoothformula, glmmformula, data, k = 5) {
  # theta: 

  ## GAM ##
  thegam <- mgcv::bam(
    smoothformula,
    data = data,
    family = binomial(),
    method = "fREML",
    discrete = TRUE, nthreads = 1
  )
  S <- thegam$smooth[[1]]$S[[1]] # Scaled penalty matrix
  lambda <- thegam$sp # Scaled smoothing parameter
  penmat <- rbind(
    0, cbind(0, lambda * S)
  ) # Add row/column of zeroes for the intercept
  # So now the penalty is beta^T penmat beta
  betastart <- coef(thegam) # Includes intercept

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
  # ustart <- rep(0,d*length(modeldata$id))
  
  control <- aghqmm::aghqmm_control(onlynllgrad = TRUE)
  
  thetastart <- c(betastart, 0)
  optfun <- function(theta) aghqmm:::optimizegammscalar(theta, yy, XX, penmat, nn, ww, control)$nll
  optgrad <- function(theta) aghqmm:::optimizegammscalar(theta, yy, XX, penmat, nn, ww, control)$grad
  opt <- optim(thetastart, optfun, optgrad, method = "BFGS")
  ## GAMM ##
  thegamm <- fit_gamm_tmb(data, 20) # Setting knots here manually.


  # Confidence intervals
  H <- numDeriv::jacobian(optgrad, opt$par)
  Hinv <- solve(H)
  paramsd <- sqrt(diag(Hinv))

  s <- length(opt$par)
  sigmaest <- exp(opt$par[s] / 2)
  sigmalower <- exp( (opt$par[s] - 2 * paramsd[s]) / 2)
  sigmaupper <- exp( (opt$par[s] + 2 * paramsd[s]) / 2)
  
  betaest <- opt$par[1:(s - 1)] # Includes intercept
  xx <- seq(0, 1, length.out = 1e03)
  Xpred <- cbind(1, mgcv::PredictMat(thegam$smooth[[1]], data.frame(x = xx)))
  ffpred <- Xpred %*% betaest
  ffpredsd <- sqrt(diag(Xpred %*% Hinv[1:(s-1), 1:(s-1)] %*% t(Xpred)))
  ffpredlower <- ffpred - 2 * ffpredsd
  ffpredupper <- ffpred + 2 * ffpredsd
  
  # Return the results from the GAM too
  ffpredgamobj <- predict(thegam, newdata = data.frame(x = xx), se.fit = TRUE)
  ffpredgam <- ffpredgamobj$fit
  ffpredsdgam <- ffpredgamobj$se.fit
  ffpredlowergam <- ffpredgam - 2 * ffpredsdgam
  ffpreduppergam <- ffpredgam + 2 * ffpredsdgam

  # Return the results from the GAMM too
  ffpredlowergamm <- thegamm$ffpred - 2 * thegamm$ffpredsd
  ffpreduppergamm <- thegamm$ffpred + 2 * thegamm$ffpredsd

  list(
    x = xx,
    sigmainterval = c("lower" = sigmalower, "estimate" = sigmaest, "upper" = sigmaupper),
    sigmagamm = thegamm$sigmaest,
    ffprednew = data.frame("lower" = ffpredlower, "estimate" = ffpred, "upper" = ffpredupper),
    ffpredgam = data.frame("lower" = ffpredlowergam, "estimate" = ffpredgam, "upper" = ffpreduppergam),
    ffpredgamm = data.frame("lower" = ffpredlowergamm, "estimate" = thegamm$ffpred, "upper" = ffpreduppergamm)
  )  
}


## Simulation study ##
numsims <- 2 # Number of simulations in each category
m <- c(1000, 2000, 5000, 10000)
n <- c(3, 9)
sigma <- c(1)
k <- c(9)

# Create simulation objects
simstodoframe <- expand.grid(
  sim = 1:numsims,
  n = n,
  m = m,
  sigma = sigma,
  k = k
)
simlist <- list()
idx <- 1
for (i in 1:nrow(simstodoframe)) {
  lst <- list(
    idx = idx,
    n = simstodoframe[i, "n"],
    m = simstodoframe[i, "m"],
    sigma = simstodoframe[i, "sigma"],
    k = simstodoframe[i, "k"]
  )

  simlist <- c(simlist, list(lst))
  idx <- idx + 1
}

options(mc.cores = parallel::detectCores())

RNGkind("L'Ecuyer-CMRG") # For reproducibility with parallel

dosim <- function(lst) {
  cat("Simulation ", lst$idx, " of ", length(simlist), ", m = ", lst$m, ", n = ", lst$n, ", k = ", lst$k, "\n", sep = "")
  df <- with(lst, simulate_data(m, n, sigma))
  res <- tryCatch(fitnewmodel(
    smoothformula = y ~ s(x, bs = "bs", k = 20),
    glmmformula = y ~ (1 | id),
    data = df,
    k = lst$k
  ), error = function(e) e, warning = function(w) w)
  if (inherits(res, "condition")) {
    cat("ERROR: Simulation ", lst$idx, " of ", length(simlist), ", m = ", lst$m, ", n = ", lst$n, ", k = ", lst$k, "\n", sep = "")
    cat("Error message: ", res$message)
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

set.seed(6431972)
if (runinparallel) mc.reset.stream() # Reproducbility in parallel
# Do the simulations
tm <- Sys.time()
if (runinparallel) {
  sims <- mclapply(simlist, dosim)
} else {
  sims <- lapply(simlist, dosim)
}
simtime <- as.numeric(difftime(Sys.time(), tm, units='secs'))
cat("Finished simulations, they took", simtime, "seconds.\n")
cat("Saving simulations...\n")
save(sims, file = file.path(resultspath, simresultsname))
cat("Saved simulations to file:",file.path(resultspath,simresultsname),"\n")
simresults <- as.data.frame(Reduce(rbind, sims))

## Load if not doing interactively
# e <- new.env()
# load(file.path(resultspath, simresultsname), envir = e)
# simresults <- as.data.frame(Reduce(rbind, e$sims))
# rm(e)

rownames(simresults) <- NULL

simresultslongbias <- simresults %>%
  tidyr::pivot_longer(contains("bias"), names_to = "type", values_to = "bias") %>%
  mutate(type = stringr::str_replace(type, "bias", ""))

simresultslongcovr <- simresults %>%
  tidyr::pivot_longer(contains("covr"), names_to = "type", values_to = "covr") %>%
  mutate(type = stringr::str_replace(type, "covr", "")) %>%
  filter(!stringr::str_detect(type, "sigma"))

PLOTTEXTSIZE <- 18
GREYLOW <- .5
GREYHIGH <- .9

n_labeller <- function(value) paste0("Group size = ", value)

biasplot <- simresultslongbias %>%
  filter(!stringr::str_detect(type, "sigma")) %>%
  ggplot(aes(x = as.factor(m), y = bias)) +
  theme_bw() +
  facet_wrap(~n, labeller = labeller(n = n_labeller)) +
  geom_boxplot(aes(fill = type)) +
  scale_fill_grey(
    start = GREYLOW, end = GREYHIGH,
    labels = c("gamm" = "GAMM", "gam" = "GAM", "new" = "PML (new)")) +
  geom_hline(aes(yintercept = 0)) +
  labs(
    title = "Empirical bias, f(x)",
    x = "Number of Groups",
    y = "Empirical bias",
    fill = "Parameter/Model"
  ) +
  theme(text = element_text(size = PLOTTEXTSIZE), legend.position = "bottom") +
  coord_cartesian(ylim = c(-.15, .15)) +
  scale_y_continuous(breaks = seq(-.15, .15, by = .025), labels = ~round(.x, 3))

ggsave(file.path(figurespath, "biasplot.pdf"), plot = biasplot, width = 7, height = 7)


covrplot <- simresultslongcovr %>%
  ggplot(aes(x = as.factor(m), y = covr)) +
  theme_bw() +
  facet_wrap(~n, labeller = labeller(n = n_labeller)) +
  geom_boxplot(aes(fill = type)) +
  geom_hline(aes(yintercept = 0.95)) +
  scale_fill_grey(start = GREYLOW, end = GREYHIGH, labels = c("gam" = "GAM", "gamm" = "GAMM", "new" = "PML (new)")) +
  labs(
    title = "Empirical coverage, f(x)",
    x = "Number of Groups",
    y = "Empirical coverage, %",
    fill = "Model"
  ) +
  scale_y_continuous(breaks = seq(0, 1, by = .1), labels = scales::percent_format()) +
  theme(text = element_text(size = PLOTTEXTSIZE), legend.position = "bottom")

ggsave(file.path(figurespath, "covrplot.pdf"), plot = covrplot, width = 7, height = 7)

# Check bias of sigma
biasplotsigma <- simresultslongbias %>%
  filter(stringr::str_detect(type, "sigma")) %>%
  ggplot(aes(x = as.factor(m), y = bias)) +
  theme_bw() +
  facet_wrap(~n, labeller = labeller(n = n_labeller)) +
  geom_boxplot(aes(fill = type)) +
  scale_fill_grey(
    start = GREYLOW, end = GREYHIGH,
    labels = c("sigmagamm" = "GAMM", "sigmanew.estimate" = "PML (new)")) +
  geom_hline(aes(yintercept = 0)) +
  labs(
    title = expression("Empirical bias,"~sigma),
    x = "Number of Groups",
    y = "Empirical bias",
    fill = "Parameter/Model"
  ) +
  theme(text = element_text(size = PLOTTEXTSIZE), legend.position = "bottom") +
  coord_cartesian(ylim = c(-.3, .3)) +
  scale_y_continuous(breaks = seq(-.3, .3, by = .05), labels = ~round(.x, 2))

ggsave(file.path(figurespath, "biasplotsigma.pdf"), plot = biasplotsigma, width = 7, height = 7)

