## Simulations for "Inference for generalized additive mixed models via penalized marginal likelihood"
## Poisson model with gamm4
## Alex Stringer
## January 2025
## File 03: Poisson model with gamm4

runinparallel <- .Platform$OS.type == "unix" # Run in parallel?

seed <- 455139 # Arbitrary

args <- commandArgs(TRUE) # Returns character(0) if interactive
if (length(args) > 0) {
  numsims <- args[1]
  simname <- args[2]
  seed <- seed * as.numeric(args[3]) # Use a different random seed for each version
} else {
  numsims <- 100
  simname <- "sims-poisson-20250228-v1"
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
  "numDeriv",
  "future.apply"
)
for (pkg in pkgs) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste0("Could not find package ", pkg, ", installing from CRAN.\n"))
    install.packages(pkg)
    require(pkg, character.only = TRUE, quietly = TRUE)
  }
}

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
    # out[i] <- (1 / 10) * (6 * dbeta(x[i], 30, 17) + 4 * dbeta(x[i], 3, 11)) - 1
    out[i] <- sin(2 * pi * x[i])
    # out[i] <- x[i]^2 + x[i]
  }
  out
}

simulate_data <- function(m, n, sigma, alpha = 0) {
  # alpha = intercept, controls Poisson mean

  if (length(n) == 1) n <- rep(n, m)
  stopifnot(length(n) == m)
  N <- sum(n)

  # generate the random effects
  u <- stats::rnorm(m, 0, sigma)

  # covariate
  # id <- rep(1:m, each = n)
  id <- Reduce(c, mapply(function(x, y) rep(x, each = y), 1:m, n, SIMPLIFY = FALSE))
  # x <- rnorm(N, .5, .5/3)
  x <- runif(N)

  # Linear predictor
  fmla <- y ~ (1 | id)
  df <- data.frame(id = id, x = x, y = 0)
  reterms <- lme4::glFormula(fmla, data = df, family = stats::binomial)
  Z <- t(reterms$reTrms$Zt)
  eta <- as.numeric(Z %*% u + ff(x))

  # Response
  lambda <- exp(eta + alpha)
  df$y <- stats::rpois(N, lambda)

  df
}

fitnewmodel <- function(smoothformula, glmmformula, data) {
  thegamm <- gamm4::gamm4(smoothformula, random = update(glmmformula, NULL ~.), data = data, family = poisson(), REML = TRUE)

  zquant <- stats::qnorm(.975)
  xx <- seq(0, 1, length.out = 1e03)
  ffpredgamm4obj <- predict(thegamm$gam, newdata = data.frame(x = xx), se.fit = TRUE, type = "terms")
  ffpredgamm4 <- ffpredgamm4obj$fit[ ,1]
  ffpredsdgamm4 <- ffpredgamm4obj$se.fit[ ,1]
  ffpredlowergamm4 <- ffpredgamm4 - zquant * ffpredsdgamm4
  ffpreduppergamm4 <- ffpredgamm4 + zquant * ffpredsdgamm4

  list(
    x = xx,
    ffpredgamm = data.frame("lower" = ffpredlowergamm4, "estimate" = ffpredgamm4, "upper" = ffpreduppergamm4),
    sigmaest = sqrt(as.numeric(summary(thegamm$mer)$varcor[1]$id))
  )
}

# df <- simulate_data(1000, 3, 1, alpha = 2)
# smoothformula = y ~ s(x, bs = "bs", k = 20)
# glmmformula = y ~ (1 | id)
# data = df
# system.time(mod <- fitnewmodel(smoothformula, glmmformula, data))

## Simulation study ##
# numsims <- 200 # Set via command line
m <- c(1000, 2000, 5000, 10000)
# n <- c(3, 9)
n <- c(3) # Average group size
sigma <- c(1)
knots <- 20
alpha <- c(-2, 0, 2)

# Create simulation objects
simstodoframe <- expand.grid(
  sim = 1:numsims,
  n = n,
  m = m,
  sigma = sigma,
  knots = knots,
  alpha = alpha
)
simlist <- list()
idx <- 1
for (i in 1:nrow(simstodoframe)) {
  lst <- list(
    idx = idx,
    n = simstodoframe[i, "n"],
    m = simstodoframe[i, "m"],
    sigma = simstodoframe[i, "sigma"],
    knots = simstodoframe[i, "knots"],
    alpha = simstodoframe[i, "alpha"]
  )

  simlist <- c(simlist, list(lst))
  idx <- idx + 1
}

dosim <- function(lst) {
  cat(simname, ": simulation ", lst$idx, " of ", length(simlist), ", m = ", lst$m, ", n = ", lst$n, ", sigma = ", lst$sigma, "\n", sep = "")
  # Average sample size = n
  nvec <- sample(2:(2 * (lst$n - 1)), size = lst$m, replace = TRUE)
  df <- with(lst, simulate_data(m, nvec, sigma, alpha))
  res <- tryCatch(fitnewmodel(
    smoothformula = y ~ s(x, bs = "bs", k = lst$knots),
    glmmformula = y ~ (1 | id),
    data = df
  ), error = function(e) e, warning = function(w) w)
  if (inherits(res, "condition")) {
    cat("ERROR: Simulation ", lst$idx, " of ", length(simlist), ", m = ", lst$m, ", n = ", lst$n, "\n", sep = "")
    cat("Error message: ", res$message, "\n")
    return(NULL)
  }

  truef <- ff(res$x)
  biasgamm <- mean(res$ffpredgamm[["estimate"]] - truef)
  
  covrgamm <- with(res, mean(ffpredgamm[["lower"]] <= truef & truef <= ffpredgamm[["upper"]]))

  biassigma <- res$sigmaest - lst$sigma


  ressum <- c("biasgamm" = biasgamm, "covrgamm" = covrgamm, "biassigma" = biassigma)

  c(unlist(lst), ressum)
}

## Simulation ##

if (runinparallel) {
  options(mc.cores = parallel::detectCores() - 1)
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
