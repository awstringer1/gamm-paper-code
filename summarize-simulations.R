## Simulations for "Inference for generalized additive mixed models via penalized marginal likelihood"
## Alex Stringer
## January 2025
## Collate and summarize simulation results, create output

studies <- c("complexfunction", "k", "main", "multiple", "poisson", "sigma", "small") # Don't change
args <- commandArgs(TRUE) # Returns character(0) if interactive
if (length(args) > 0) {
  simdate <- args[1]s
} else {
  simdate <- "20250306"
}

## Packages ##
pkgs <- c(
  "dplyr",
  "tibble",
  "ggplot2",
  "tidyr",
  "stringr"
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
resultspath <- file.path(basepath, "results") # Where to load the results from
if (!dir.exists(resultspath)) dir.create(resultspath)
figurespath <- file.path(basepath, "figures") # Where to save the figures to
if (!dir.exists(figurespath)) dir.create(figurespath)


## Load the simulations ##
# Load all the simulations with "simdate"
files <- list.files(resultspath)
files <- files[grep(simdate, files)]


e <- new.env()
simresultslist <- list()
groups <- character()
filesloaded <- character()
i <- 1
for (file in files) {
  fgroup <- file %>%
    stringr::str_remove("sims") %>%
    stringr::str_remove_all("-") %>%
    stringr::str_remove("v[0-9]") %>%
    stringr::str_remove(".RData") %>%
    stringr::str_remove(simdate)
  if (fgroup == "") next
  groups <- c(groups, fgroup)
  filesloaded <- c(filesloaded, file)

  load(file.path(resultspath, file), envir = e) # Loads object "sims" into environment e
  simresults <- Reduce(rbind, e$sims)
  simresultslist <- c(simresultslist, list(simresults))
  names(simresultslist)[i] <- fgroup
  i <- i + 1
  rm(list = ls(envir = e), envir = e) # Clear e
}
rm(e)

groups <- unique(groups)
stopifnot(setequal(names(simresultslist), studies))
stopifnot(setequal(names(simresultslist), groups))

versions <- stringr::str_extract(filesloaded, "v[0-9]")
# table(versions) # Which version files were loaded

# Collate them by group
for (study in studies) {
  sims <- simresultslist[which(names(simresultslist) == study)]
  simresults <- tibble::as_tibble(Reduce(rbind, sims))
  rownames(simresults) <- NULL
  assign(paste0("simresults_", study), simresults)
}
rm(sims, simresults)


## Summarize the simulations ##

PLOTTEXTSIZE <- 18
ERRORBARWIDTH <- .1
col_gam <- "#717171"
col_gamm <- "#ACACAC"
col_pml <- "#DADADA"

n_labeller <- function(value) paste0("Avg. group size = ", value)
k_labeller <- function(value) paste0("Quad. points = ", value)


## Main study ##

# Doesn't appear in paper; just for visual check
simresults_main_table <- simresults_main %>%
  group_by(n, m, sigma, k) %>%
  summarize(numsims = n())
  
simresultslongbias_main <- simresults_main %>%
  tidyr::pivot_longer(contains("bias"), names_to = "type", values_to = "bias") %>%
  mutate(type = stringr::str_replace(type, "bias", ""))

simresultslongcovr_main <- simresults_main %>%
  tidyr::pivot_longer(contains("covr"), names_to = "type", values_to = "covr") %>%
  mutate(type = stringr::str_replace(type, "covr", "")) %>%
  filter(!stringr::str_detect(type, "sigma"))


biasplot_main <- simresultslongbias_main %>%
  filter(!stringr::str_detect(type, "sigma")) %>%
  ggplot(aes(x = as.factor(m), y = bias)) +
  theme_bw() +
  facet_grid(k~n, labeller = labeller(k = k_labeller, n = n_labeller)) +
  geom_boxplot(aes(fill = type)) +
  scale_fill_manual(
    values = c("gamm" = col_gamm, "gam" = col_gam, "new" = col_pml),
    labels = c("gamm" = "GAMM", "gam" = "GAM", "new" = "PML (new)")) +
  geom_hline(aes(yintercept = 0)) +
  labs(
    title = "Empirical bias, f(x)",
    x = "Number of Groups",
    y = "Empirical bias",
    fill = "Parameter/Model"
  ) +
  theme(text = element_text(size = PLOTTEXTSIZE), legend.position = "bottom") +
  coord_cartesian(ylim = c(-.05, .05)) +
  scale_y_continuous(breaks = seq(-.05, .05, by = .01), labels = ~round(.x, 3))

biasplotsigma_main <- simresultslongbias_main %>%
  filter(stringr::str_detect(type, "sigma")) %>%
  ggplot(aes(x = as.factor(m), y = bias)) +
  theme_bw() +
  facet_grid(k~n, labeller = labeller(k = k_labeller, n = n_labeller)) +
  geom_boxplot(aes(fill = type)) +
  scale_fill_manual(
    values = c("sigmagamm" = col_gamm, "sigmanew.estimate" = col_pml),
    labels = c("sigmagamm" = "GAMM", "sigmanew.estimate" = "PML (new)")) +
  geom_hline(aes(yintercept = 0)) +
  labs(
    title = expression("Empirical bias,"~sigma),
    x = "Number of Groups",
    y = "Empirical bias",
    fill = "Parameter/Model"
  ) +
  theme(text = element_text(size = PLOTTEXTSIZE), legend.position = "bottom") +
  coord_cartesian(ylim = c(-.4, .25)) +
  scale_y_continuous(breaks = seq(-.4, .25, by = .05), labels = ~round(.x, 2))

covrlineplot_main <- simresultslongcovr_main %>%
  group_by(m, n, sigma, k, type) %>%
  summarize(covrmean = mean(covr), covrse = sd(covr) / sqrt(n())) %>%
  ggplot(aes(x = as.factor(m), y = covrmean, group = type, linetype = type)) +
  theme_bw() +
  facet_grid(k~n, labeller = labeller(k = k_labeller, n = n_labeller)) +
  geom_line() +
  geom_errorbar(aes(ymin = covrmean - 2 * covrse, ymax = covrmean + 2 * covrse), width = ERRORBARWIDTH) +
  geom_hline(aes(yintercept = 0.95)) +
  scale_linetype_discrete(labels = c("gam" = "GAM", "gamm" = "GAMM", "new" = "PML (new)")) +
  labs(
    title = "Empirical coverage, f(x)",
    x = "Number of Groups",
    y = "Empirical coverage, %",
    linetype = "Model"
  ) +
  scale_y_continuous(breaks = seq(0, 1, by = .1), labels = scales::percent_format()) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(text = element_text(size = PLOTTEXTSIZE), legend.position = "bottom")

## Varying sigma ##

# Construct the labeller for sigma
sigma_values <- as.character(unique(simresults_sigma$sigma))
sigma_labels <- paste0("sigma == ", sigma_values)
names(sigma_labels) <- sigma_values
s_labeller <- as_labeller(sigma_labels, default = label_parsed)

simresults_sigma_table <- simresults_sigma %>%
  group_by(n, m, sigma, k) %>%
  summarize(numsims = n())

simresultslongbias_sigma <- simresults_sigma %>%
  tidyr::pivot_longer(contains("bias"), names_to = "type", values_to = "bias") %>%
  mutate(type = stringr::str_replace(type, "bias", ""))

simresultslongcovr_sigma <- simresults_sigma %>%
  tidyr::pivot_longer(contains("covr"), names_to = "type", values_to = "covr") %>%
  mutate(type = stringr::str_replace(type, "covr", "")) %>%
  filter(!stringr::str_detect(type, "sigma"))


biasplot_sigma <- simresultslongbias_sigma %>%
  filter(!stringr::str_detect(type, "sigma")) %>%
  ggplot(aes(x = as.factor(m), y = bias)) +
  theme_bw() +
  facet_grid(sigma~n, labeller = labeller(sigma = s_labeller, n = n_labeller)) +
  geom_boxplot(aes(fill = type)) +
  scale_fill_manual(
    values = c("gamm" = col_gamm, "gam" = col_gam, "new" = col_pml),
    labels = c("gamm" = "GAMM", "gam" = "GAM", "new" = "PML (new)")) +
  geom_hline(aes(yintercept = 0)) +
  labs(
    title = "Empirical bias, f(x)",
    x = "Number of Groups",
    y = "Empirical bias",
    fill = "Parameter/Model"
  ) +
  theme(text = element_text(size = PLOTTEXTSIZE), legend.position = "bottom") +
  coord_cartesian(ylim = c(-.05, .05)) +
  scale_y_continuous(breaks = seq(-.05, .05, by = .025), labels = ~round(.x, 3))

biasplotsigma_sigma <- simresultslongbias_sigma %>%
  filter(stringr::str_detect(type, "sigma")) %>%
  ggplot(aes(x = as.factor(m), y = bias)) +
  theme_bw() +
  facet_grid(sigma ~ n, labeller = labeller(sigma = s_labeller, n = n_labeller)) +
  geom_boxplot(aes(fill = type)) +
  scale_fill_manual(
    values = c("sigmagamm" = col_gamm, "sigmanew.estimate" = col_pml),
    labels = c("sigmagamm" = "GAMM", "sigmanew.estimate" = "PML (new)")) +
  geom_hline(aes(yintercept = 0)) +
  labs(
    title = expression("Empirical bias,"~sigma),
    x = "Number of Groups",
    y = "Empirical bias",
    fill = "Parameter/Model"
  ) +
  theme(text = element_text(size = PLOTTEXTSIZE), legend.position = "bottom") +
  coord_cartesian(ylim = c(-.5, .2)) +
  scale_y_continuous(breaks = seq(-.5, .2, by = .1), labels = ~round(.x, 2))

covrlineplot_sigma <- simresultslongcovr_sigma %>%
  group_by(m, n, sigma, k, type) %>%
  summarize(covrmean = mean(covr), covrse = sd(covr) / sqrt(n())) %>%
  ggplot(aes(x = as.factor(m), y = covrmean, group = type, linetype = type)) +
  theme_bw() +
  facet_grid(sigma ~ n, labeller = labeller(sigma = s_labeller, n = n_labeller)) +
  geom_line() +
  geom_errorbar(aes(ymin = covrmean - 2 * covrse, ymax = covrmean + 2 * covrse), width = ERRORBARWIDTH) +
  geom_hline(aes(yintercept = 0.95)) +
  scale_linetype_discrete(labels = c("gam" = "GAM", "gamm" = "GAMM", "new" = "PML (new)")) +
  labs(
    title = "Empirical coverage, f(x)",
    x = "Number of Groups",
    y = "Empirical coverage, %",
    linetype = "Model"
  ) +
  scale_y_continuous(breaks = seq(0, 1, by = .2), labels = scales::percent_format()) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(text = element_text(size = PLOTTEXTSIZE), legend.position = "bottom")

## Varying k ##

simresults_k_table <- simresults_k %>%
  group_by(n, m, sigma, k) %>%
  summarize(numsims = n())
  
simresultslongbias_k <- simresults_k %>%
  tidyr::pivot_longer(contains("bias"), names_to = "type", values_to = "bias") %>%
  mutate(type = stringr::str_replace(type, "bias", ""))

simresultslongcovr_k <- simresults_k %>%
  tidyr::pivot_longer(contains("covr"), names_to = "type", values_to = "covr") %>%
  mutate(type = stringr::str_replace(type, "covr", "")) %>%
  filter(!stringr::str_detect(type, "sigma"))


biasplot_k <- simresultslongbias_k %>%
  filter(!stringr::str_detect(type, "sigma")) %>%
  ggplot(aes(x = as.factor(m), y = bias)) +
  theme_bw() +
  facet_grid(k~n, labeller = labeller(k = k_labeller, n = n_labeller)) +
  geom_boxplot(aes(fill = type)) +
  scale_fill_manual(
    values = c("gamm" = col_gamm, "gam" = col_gam, "new" = col_pml),
    labels = c("gamm" = "GAMM", "gam" = "GAM", "new" = "PML (new)")) +
  geom_hline(aes(yintercept = 0)) +
  labs(
    title = "Empirical bias, f(x)",
    x = "Number of Groups",
    y = "Empirical bias",
    fill = "Parameter/Model"
  ) +
  theme(text = element_text(size = PLOTTEXTSIZE), legend.position = "bottom") +
  coord_cartesian(ylim = c(-.04, .04)) +
  scale_y_continuous(breaks = seq(-.04, .04, by = .01), labels = ~round(.x, 3))

biasplotsigma_k <- simresultslongbias_k %>%
  filter(stringr::str_detect(type, "sigma")) %>%
  ggplot(aes(x = as.factor(m), y = bias)) +
  theme_bw() +
  facet_grid(k~n, labeller = labeller(k = k_labeller, n = n_labeller)) +
  geom_boxplot(aes(fill = type)) +
  scale_fill_manual(
    values = c("sigmagamm" = col_gamm, "sigmanew.estimate" = col_pml),
    labels = c("sigmagamm" = "GAMM", "sigmanew.estimate" = "PML (new)")) +
  geom_hline(aes(yintercept = 0)) +
  labs(
    title = expression("Empirical bias,"~sigma),
    x = "Number of Groups",
    y = "Empirical bias",
    fill = "Parameter/Model"
  ) +
  theme(text = element_text(size = PLOTTEXTSIZE), legend.position = "bottom") +
  coord_cartesian(ylim = c(-.4, .2)) +
  scale_y_continuous(breaks = seq(-.4, .2, by = .05), labels = ~round(.x, 2))

covrlineplot_k <- simresultslongcovr_k %>%
  group_by(m, n, sigma, k, type) %>%
  summarize(covrmean = mean(covr), covrse = sd(covr) / sqrt(n())) %>%
  ggplot(aes(x = as.factor(m), y = covrmean, group = type, linetype = type)) +
  theme_bw() +
  facet_grid(k~n, labeller = labeller(k = k_labeller, n = n_labeller)) +
  geom_line() +
  geom_errorbar(aes(ymin = covrmean - 2 * covrse, ymax = covrmean + 2 * covrse), width = ERRORBARWIDTH) +
  geom_hline(aes(yintercept = 0.95)) +
  scale_linetype_discrete(labels = c("gam" = "GAM", "gamm" = "GAMM", "new" = "PML (new)")) +
  labs(
    title = "Empirical coverage, f(x)",
    x = "Number of Groups",
    y = "Empirical coverage, %",
    linetype = "Model"
  ) +
  scale_y_continuous(breaks = seq(0, 1, by = .2), labels = scales::percent_format()) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(text = element_text(size = PLOTTEXTSIZE), legend.position = "bottom")

## More wiggly function ##

simresults_complexfunction_table <- simresults_complexfunction %>%
  group_by(n, m, sigma, k) %>%
  summarize(numsims = n())
  
simresultslongbias_complexfunction <- simresults_complexfunction %>%
  tidyr::pivot_longer(contains("bias"), names_to = "type", values_to = "bias") %>%
  mutate(type = stringr::str_replace(type, "bias", ""))

simresultslongcovr_complexfunction <- simresults_complexfunction %>%
  tidyr::pivot_longer(contains("covr"), names_to = "type", values_to = "covr") %>%
  mutate(type = stringr::str_replace(type, "covr", "")) %>%
  filter(!stringr::str_detect(type, "sigma"))


biasplot_complexfunction <- simresultslongbias_complexfunction %>%
  filter(!stringr::str_detect(type, "sigma")) %>%
  ggplot(aes(x = as.factor(m), y = bias)) +
  theme_bw() +
  facet_grid(k~n, labeller = labeller(k = k_labeller, n = n_labeller)) +
  geom_boxplot(aes(fill = type)) +
  scale_fill_manual(
    values = c("gamm" = col_gamm, "gam" = col_gam, "new" = col_pml),
    labels = c("gamm" = "GAMM", "gam" = "GAM", "new" = "PML (new)")) +
  geom_hline(aes(yintercept = 0)) +
  labs(
    title = "Empirical bias, f(x)",
    x = "Number of Groups",
    y = "Empirical bias",
    fill = "Parameter/Model"
  ) +
  theme(text = element_text(size = PLOTTEXTSIZE), legend.position = "bottom") +
  coord_cartesian(ylim = c(-.2, .1)) +
  scale_y_continuous(breaks = seq(-.2, .1, by = .025), labels = ~round(.x, 3))

biasplotsigma_complexfunction <- simresultslongbias_complexfunction %>%
  filter(stringr::str_detect(type, "sigma")) %>%
  ggplot(aes(x = as.factor(m), y = bias)) +
  theme_bw() +
  facet_grid(k~n, labeller = labeller(k = k_labeller, n = n_labeller)) +
  geom_boxplot(aes(fill = type)) +
  scale_fill_manual(
    values = c("sigmagamm" = col_gamm, "sigmanew.estimate" = col_pml),
    labels = c("sigmagamm" = "GAMM", "sigmanew.estimate" = "PML (new)")) +
  geom_hline(aes(yintercept = 0)) +
  labs(
    title = expression("Empirical bias,"~sigma),
    x = "Number of Groups",
    y = "Empirical bias",
    fill = "Parameter/Model"
  ) +
  theme(text = element_text(size = PLOTTEXTSIZE), legend.position = "bottom") +
  coord_cartesian(ylim = c(-.4, .2)) +
  scale_y_continuous(breaks = seq(-.4, .2, by = .05), labels = ~round(.x, 2))

covrlineplot_complexfunction <- simresultslongcovr_complexfunction %>%
  group_by(m, n, sigma, k, type) %>%
  summarize(covrmean = mean(covr), covrse = sd(covr) / sqrt(n())) %>%
  ggplot(aes(x = as.factor(m), y = covrmean, group = type, linetype = type)) +
  theme_bw() +
  facet_grid(k~n, labeller = labeller(k = k_labeller, n = n_labeller)) +
  geom_line() +
  geom_errorbar(aes(ymin = covrmean - 2 * covrse, ymax = covrmean + 2 * covrse), width = ERRORBARWIDTH) +
  geom_hline(aes(yintercept = 0.95)) +
  scale_linetype_discrete(labels = c("gam" = "GAM", "gamm" = "GAMM", "new" = "PML (new)")) +
  labs(
    title = "Empirical coverage, f(x)",
    x = "Number of Groups",
    y = "Empirical coverage, %",
    linetype = "Model"
  ) +
  scale_y_continuous(breaks = seq(0, 1, by = .1), labels = scales::percent_format()) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(text = element_text(size = PLOTTEXTSIZE), legend.position = "bottom")


## Poisson model with GAMM4 only ##

# Construct the labeller for alpha
alpha_values <- as.character(unique(simresults_poisson$alpha))
alpha_labels <- paste0("alpha == ", alpha_values)
names(alpha_labels) <- alpha_values
a_labeller <- as_labeller(alpha_labels, default = label_parsed)

simresults_poisson_table <- simresults_poisson %>%
  group_by(n, m, sigma, alpha) %>%
  summarize(numsims = n())
  
simresultslongbias_poisson <- simresults_poisson %>%
  tidyr::pivot_longer(contains("bias"), names_to = "type", values_to = "bias") %>%
  mutate(type = stringr::str_replace(type, "bias", ""))

simresultslongcovr_poisson <- simresults_poisson %>%
  tidyr::pivot_longer(contains("covr"), names_to = "type", values_to = "covr") %>%
  mutate(type = stringr::str_replace(type, "covr", "")) %>%
  filter(!stringr::str_detect(type, "sigma"))

biasplot_poisson <- simresultslongbias_poisson %>%
  filter(stringr::str_detect(type, "gamm")) %>%
  ggplot(aes(x = as.factor(m), y = bias)) +
  theme_bw() +
  facet_grid(alpha ~ n, labeller = labeller(alpha = a_labeller, n = n_labeller)) +
  geom_boxplot(aes(fill = type)) +
  scale_fill_manual(
    values = c("gamm" = col_gamm),
    labels = c("gamm" = "GAMM")) +
  geom_hline(aes(yintercept = 0)) +
  labs(
    title = "Empirical bias, f(x)",
    x = "Number of Groups",
    y = "Empirical bias",
    fill = "Parameter/Model"
  ) +
  theme(text = element_text(size = PLOTTEXTSIZE), legend.position = "bottom") +
  coord_cartesian(ylim = c(-.05, .05)) +
  scale_y_continuous(breaks = seq(-.05, .05, by = .02), labels = ~round(.x, 3))

biasplotsigma_poisson <- simresultslongbias_poisson %>%
  filter(stringr::str_detect(type, "sigma")) %>%
  ggplot(aes(x = as.factor(m), y = bias)) +
  theme_bw() +
  facet_grid(alpha ~ n, labeller = labeller(alpha = a_labeller, n = n_labeller)) +
  geom_boxplot(aes(fill = type)) +
  scale_fill_manual(
    values = c("sigma" = col_gamm),
    labels = c("sigma" = "GAMM")) +
  geom_hline(aes(yintercept = 0)) +
  labs(
    title = expression("Empirical bias,"~sigma),
    x = "Number of Groups",
    y = "Empirical bias",
    fill = "Parameter/Model"
  ) +
  theme(text = element_text(size = PLOTTEXTSIZE), legend.position = "bottom") +
  coord_cartesian(ylim = c(-.15, .2)) +
  scale_y_continuous(breaks = seq(-.15, .2, by = .05), labels = ~round(.x, 2))


covrlineplot_poisson <- simresultslongcovr_poisson %>%
  group_by(m, n, sigma, type, alpha) %>%
  summarize(covrmean = mean(covr), covrse = sd(covr) / sqrt(n())) %>%
  ggplot(aes(x = as.factor(m), y = covrmean, group = type, linetype = type)) +
  theme_bw() +
  facet_grid(alpha ~ n, labeller = labeller(alpha = a_labeller, n = n_labeller)) +
  geom_line() +
  geom_errorbar(aes(ymin = covrmean - 2 * covrse, ymax = covrmean + 2 * covrse), width = ERRORBARWIDTH) +
  geom_hline(aes(yintercept = 0.95)) +
  scale_linetype_discrete(labels = c("gam" = "GAM", "gamm" = "GAMM", "new" = "PML (new)")) +
  labs(
    title = "Empirical coverage, f(x)",
    x = "Number of Groups",
    y = "Empirical coverage, %",
    linetype = "Model"
  ) +
  scale_y_continuous(breaks = seq(0, 1, by = .2), labels = scales::percent_format()) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(text = element_text(size = PLOTTEXTSIZE), legend.position = "bottom")


## Multiple smooth functions ##

simresults_multiple_table <- simresults_multiple %>%
  group_by(n, m, sigma, k) %>%
  summarize(numsims = n())
  
  
simresultslongbias_multiple <- simresults_multiple %>%
  tidyr::pivot_longer(contains("bias"), names_to = "type", values_to = "bias") %>%
  mutate(
    fun = stringr::str_extract(type, "[0-9]"),
    type = stringr::str_replace(type, "bias[0-9]?", "")
  )

simresultslongcovr_multiple <- simresults_multiple %>%
  tidyr::pivot_longer(contains("covr"), names_to = "type", values_to = "covr") %>%
  mutate(
    fun = stringr::str_extract(type, "[0-9]"),
    type = stringr::str_replace(type, "bias[0-9]?", "")
  ) %>%
  filter(!stringr::str_detect(type, "sigma"))


biasplot_multiple_f1 <- simresultslongbias_multiple %>%
  filter(fun == "1") %>%
  ggplot(aes(x = as.factor(m), y = bias)) +
  theme_bw() +
  facet_grid(k~n, labeller = labeller(k = k_labeller, n = n_labeller)) +
  geom_boxplot(aes(fill = type)) +
  scale_fill_manual(
    values = c("gamm" = col_gamm, "gam" = col_gam, "new" = col_pml),
    labels = c("gamm" = "GAMM", "gam" = "GAM", "new" = "PML (new)")) +
  geom_hline(aes(yintercept = 0)) +
  labs(
    title = expression("Empirical bias, f"[1]*"(x)"),
    x = "Number of Groups",
    y = "Empirical bias",
    fill = "Parameter/Model"
  ) +
  theme(text = element_text(size = PLOTTEXTSIZE), legend.position = "bottom") +
  coord_cartesian(ylim = c(-.05, .03)) +
  scale_y_continuous(breaks = seq(-.05, .03, by = .005), labels = ~round(.x, 3))

biasplot_multiple_f2 <- simresultslongbias_multiple %>%
  filter(fun == "2") %>%
  ggplot(aes(x = as.factor(m), y = bias)) +
  theme_bw() +
  facet_grid(k~n, labeller = labeller(k = k_labeller, n = n_labeller)) +
  geom_boxplot(aes(fill = type)) +
  scale_fill_manual(
    values = c("gamm" = col_gamm, "gam" = col_gam, "new" = col_pml),
    labels = c("gamm" = "GAMM", "gam" = "GAM", "new" = "PML (new)")) +
  geom_hline(aes(yintercept = 0)) +
  labs(
    title = expression("Empirical bias, f"[2]*"(x)"),
    x = "Number of Groups",
    y = "Empirical bias",
    fill = "Parameter/Model"
  ) +
  theme(text = element_text(size = PLOTTEXTSIZE), legend.position = "bottom") +
  coord_cartesian(ylim = c(-.05, .03)) +
  scale_y_continuous(breaks = seq(-.05, .03, by = .005), labels = ~round(.x, 3))

biasplotsigma_multiple <- simresultslongbias_multiple %>%
  filter(stringr::str_detect(type, "sigma")) %>%
  ggplot(aes(x = as.factor(m), y = bias)) +
  theme_bw() +
  facet_grid(k~n, labeller = labeller(k = k_labeller, n = n_labeller)) +
  geom_boxplot(aes(fill = type)) +
  scale_fill_manual(
    values = c("sigmagamm" = col_gamm, "sigmanew.estimate" = col_pml),
    labels = c("sigmagamm" = "GAMM", "sigmanew.estimate" = "PML (new)")) +
  geom_hline(aes(yintercept = 0)) +
  labs(
    title = expression("Empirical bias,"~sigma),
    x = "Number of Groups",
    y = "Empirical bias",
    fill = "Parameter/Model"
  ) +
  theme(text = element_text(size = PLOTTEXTSIZE), legend.position = "bottom") +
  coord_cartesian(ylim = c(-.4, .2)) +
  scale_y_continuous(breaks = seq(-.4, .2, by = .05), labels = ~round(.x, 2))

covrlineplot_multiple_f1 <- simresultslongcovr_multiple %>%
  filter(fun == "1") %>%
  group_by(m, n, sigma, k, type) %>%
  summarize(covrmean = mean(covr), covrse = sd(covr) / sqrt(n())) %>%
  ggplot(aes(x = as.factor(m), y = covrmean, group = type, linetype = type)) +
  theme_bw() +
  facet_grid(k~n, labeller = labeller(k = k_labeller, n = n_labeller)) +
  geom_line() +
  geom_errorbar(aes(ymin = covrmean - 2 * covrse, ymax = covrmean + 2 * covrse), width = ERRORBARWIDTH) +
  geom_hline(aes(yintercept = 0.95)) +
  scale_linetype_discrete(labels = c("covr1gam" = "GAM", "covr1gamm" = "GAMM", "covr1new" = "PML (new)")) +
  labs(
    title = expression("Empirical coverage, f"[1]*"(x)"),
    x = "Number of Groups",
    y = "Empirical coverage, %",
    linetype = "Model"
  ) +
  scale_y_continuous(breaks = seq(0, 1, by = .1), labels = scales::percent_format()) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(text = element_text(size = PLOTTEXTSIZE), legend.position = "bottom")

covrlineplot_multiple_f2 <- simresultslongcovr_multiple %>%
  filter(fun == "2") %>%
  group_by(m, n, sigma, k, type) %>%
  summarize(covrmean = mean(covr), covrse = sd(covr) / sqrt(n())) %>%
  ggplot(aes(x = as.factor(m), y = covrmean, group = type, linetype = type)) +
  theme_bw() +
  facet_grid(k~n, labeller = labeller(k = k_labeller, n = n_labeller)) +
  geom_line() +
  geom_errorbar(aes(ymin = covrmean - 2 * covrse, ymax = covrmean + 2 * covrse), width = ERRORBARWIDTH) +
  geom_hline(aes(yintercept = 0.95)) +
  scale_linetype_discrete(labels = c("covr2gam" = "GAM", "covr2gamm" = "GAMM", "covr2new" = "PML (new)")) +
  labs(
    title = expression("Empirical coverage, f"[2]*"(x)"),
    x = "Number of Groups",
    y = "Empirical coverage, %",
    linetype = "Model"
  ) +
  scale_y_continuous(breaks = seq(0, 1, by = .1), labels = scales::percent_format()) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(text = element_text(size = PLOTTEXTSIZE), legend.position = "bottom")
  theme(text = element_text(size = PLOTTEXTSIZE), legend.position = "bottom")

## Small m ##

simresults_small_table <- simresults_small %>%
  group_by(n, m, sigma, k) %>%
  summarize(numsims = n())
  
simresultslongbias_small <- simresults_small %>%
  tidyr::pivot_longer(contains("bias"), names_to = "type", values_to = "bias") %>%
  mutate(type = stringr::str_replace(type, "bias", ""))

simresultslongcovr_small <- simresults_small %>%
  tidyr::pivot_longer(contains("covr"), names_to = "type", values_to = "covr") %>%
  mutate(type = stringr::str_replace(type, "covr", "")) %>%
  filter(!stringr::str_detect(type, "sigma"))


biasplot_small <- simresultslongbias_small %>%
  filter(!stringr::str_detect(type, "sigma")) %>%
  ggplot(aes(x = as.factor(m), y = bias)) +
  theme_bw() +
  facet_grid(k~n, labeller = labeller(k = k_labeller, n = n_labeller)) +
  geom_boxplot(aes(fill = type)) +
  scale_fill_manual(
    values = c("gamm" = col_gamm, "gam" = col_gam, "new" = col_pml),
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

biasplotsigma_small <- simresultslongbias_small %>%
  filter(stringr::str_detect(type, "sigma")) %>%
  ggplot(aes(x = as.factor(m), y = bias)) +
  theme_bw() +
  facet_grid(k~n, labeller = labeller(k = k_labeller, n = n_labeller)) +
  geom_boxplot(aes(fill = type)) +
  scale_fill_manual(
    values = c("sigmagamm" = col_gamm, "sigmanew.estimate" = col_pml),
    labels = c("sigmagamm" = "GAMM", "sigmanew.estimate" = "PML (new)")) +
  geom_hline(aes(yintercept = 0)) +
  labs(
    title = expression("Empirical bias,"~sigma),
    x = "Number of Groups",
    y = "Empirical bias",
    fill = "Parameter/Model"
  ) +
  theme(text = element_text(size = PLOTTEXTSIZE), legend.position = "bottom") +
  coord_cartesian(ylim = c(-.7, .7)) +
  scale_y_continuous(breaks = seq(-.7, .7, by = .1), labels = ~round(.x, 2))

covrlineplot_small <- simresultslongcovr_small %>%
  group_by(m, n, sigma, k, type) %>%
  summarize(covrmean = mean(covr), covrse = sd(covr) / sqrt(n())) %>%
  ggplot(aes(x = as.factor(m), y = covrmean, group = type, linetype = type)) +
  theme_bw() +
  facet_grid(k~n, labeller = labeller(k = k_labeller, n = n_labeller)) +
  geom_line() +
  geom_errorbar(aes(ymin = covrmean - 2 * covrse, ymax = covrmean + 2 * covrse), width = ERRORBARWIDTH) +
  geom_hline(aes(yintercept = 0.95)) +
  scale_linetype_discrete(labels = c("gam" = "GAM", "gamm" = "GAMM", "new" = "PML (new)")) +
  labs(
    title = "Empirical coverage, f(x)",
    x = "Number of Groups",
    y = "Empirical coverage, %",
    linetype = "Model"
  ) +
  scale_y_continuous(breaks = seq(0, 1, by = .1), labels = scales::percent_format()) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(text = element_text(size = PLOTTEXTSIZE), legend.position = "bottom")


## Save figures ##
FIGUREWIDTH = 7
FIGUREHEIGHT = 7

# Main
ggsave(filename = file.path(figurespath, paste0("figure-", simdate, "-main-bias.pdf")), plot = biasplot_main, width = FIGUREWIDTH, height = FIGUREHEIGHT)
ggsave(filename = file.path(figurespath, paste0("figure-", simdate, "-main-sigma.pdf")), plot = biasplotsigma_main, width = FIGUREWIDTH, height = FIGUREHEIGHT)
ggsave(filename = file.path(figurespath, paste0("figure-", simdate, "-main-covr.pdf")), plot = covrlineplot_main, width = FIGUREWIDTH, height = FIGUREHEIGHT)

# Sigma
ggsave(filename = file.path(figurespath, paste0("figure-", simdate, "-sigma-bias.pdf")), plot = biasplot_sigma, width = FIGUREWIDTH, height = FIGUREHEIGHT)
ggsave(filename = file.path(figurespath, paste0("figure-", simdate, "-sigma-sigma.pdf")), plot = biasplotsigma_sigma, width = FIGUREWIDTH, height = FIGUREHEIGHT)
ggsave(filename = file.path(figurespath, paste0("figure-", simdate, "-sigma-covr.pdf")), plot = covrlineplot_sigma, width = FIGUREWIDTH, height = FIGUREHEIGHT)

# k
ggsave(filename = file.path(figurespath, paste0("figure-", simdate, "-k-bias.pdf")), plot = biasplot_k, width = FIGUREWIDTH, height = FIGUREHEIGHT)
ggsave(filename = file.path(figurespath, paste0("figure-", simdate, "-k-sigma.pdf")), plot = biasplotsigma_k, width = FIGUREWIDTH, height = FIGUREHEIGHT)
ggsave(filename = file.path(figurespath, paste0("figure-", simdate, "-k-covr.pdf")), plot = covrlineplot_k, width = FIGUREWIDTH, height = FIGUREHEIGHT)

# Wiggly
ggsave(filename = file.path(figurespath, paste0("figure-", simdate, "-complexfunction-bias.pdf")), plot = biasplot_complexfunction, width = FIGUREWIDTH, height = FIGUREHEIGHT)
ggsave(filename = file.path(figurespath, paste0("figure-", simdate, "-complexfunction-sigma.pdf")), plot = biasplotsigma_complexfunction, width = FIGUREWIDTH, height = FIGUREHEIGHT)
ggsave(filename = file.path(figurespath, paste0("figure-", simdate, "-complexfunction-covr.pdf")), plot = covrlineplot_complexfunction, width = FIGUREWIDTH, height = FIGUREHEIGHT)

# Poisson
ggsave(filename = file.path(figurespath, paste0("figure-", simdate, "-poisson-bias.pdf")), plot = biasplot_poisson, width = FIGUREWIDTH, height = FIGUREHEIGHT)
ggsave(filename = file.path(figurespath, paste0("figure-", simdate, "-poisson-sigma.pdf")), plot = biasplotsigma_poisson, width = FIGUREWIDTH, height = FIGUREHEIGHT)
ggsave(filename = file.path(figurespath, paste0("figure-", simdate, "-poisson-covr.pdf")), plot = covrlineplot_poisson, width = FIGUREWIDTH, height = FIGUREHEIGHT)

# Multiple functions
ggsave(filename = file.path(figurespath, paste0("figure-", simdate, "-multiple-f1-bias.pdf")), plot = biasplot_multiple_f1, width = FIGUREWIDTH, height = FIGUREHEIGHT)
ggsave(filename = file.path(figurespath, paste0("figure-", simdate, "-multiple-f2-bias.pdf")), plot = biasplot_multiple_f2, width = FIGUREWIDTH, height = FIGUREHEIGHT)
ggsave(filename = file.path(figurespath, paste0("figure-", simdate, "-multiple-sigma.pdf")), plot = biasplotsigma_multiple, width = FIGUREWIDTH, height = FIGUREHEIGHT)
ggsave(filename = file.path(figurespath, paste0("figure-", simdate, "-multiple-f1-covr.pdf")), plot = covrlineplot_multiple_f1, width = FIGUREWIDTH, height = FIGUREHEIGHT)
ggsave(filename = file.path(figurespath, paste0("figure-", simdate, "-multiple-f2-covr.pdf")), plot = covrlineplot_multiple_f2, width = FIGUREWIDTH, height = FIGUREHEIGHT)

# Small m
ggsave(filename = file.path(figurespath, paste0("figure-", simdate, "-small-bias.pdf")), plot = biasplot_small, width = FIGUREWIDTH, height = FIGUREHEIGHT)
ggsave(filename = file.path(figurespath, paste0("figure-", simdate, "-small-sigma.pdf")), plot = biasplotsigma_small, width = FIGUREWIDTH, height = FIGUREHEIGHT)
ggsave(filename = file.path(figurespath, paste0("figure-", simdate, "-small-covr.pdf")), plot = covrlineplot_small, width = FIGUREWIDTH, height = FIGUREHEIGHT)
