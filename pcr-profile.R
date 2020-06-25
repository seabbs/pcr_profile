library(data.table)

# SIMULATE DATA
# Number of patients
n_patients <- 20

# Generate patients with random number of tests performed
df <- data.table(pID = rep(1:n_patients, rpois(n_patients, 4)))

# Generate random actual exposure time for each patient
df[, actual_te := rnbinom(1, mu = 15, size = 15), by = "pID"]

# Parameters for fake incubation period
incub_alpha <- 2

incub_beta <- 0.75

plot(dgamma(1:10, 2, 0.75))

# Generate onset based on exposure + incubation period
df[, onset := mean(actual_te) + round(rgamma(1, shape = incub_alpha, rate = incub_beta)), by = "pID"]

# Generate times when tests were performed 
df[, dayoftest := sample((mean(onset) - 4):(mean(onset)+12), replace = FALSE, size = .N), by = "pID"]

df[order(pID, dayoftest)]

# If test after onset, symptomatic (keeping it simple for now)
df$symptomatic <- ifelse(df$dayoftest >= df$onset, TRUE, FALSE)

# Probability that test is positive based on test date and exposure date
df[, tp := fifelse(dayoftest <= actual_te, 0, 1 - pgamma(dayoftest - actual_te, 2, 0.2))]

# TRUE/FALSE return for tests
df$test_result <- purrr::rbernoulli(n = nrow(df), p = df$tp)

# Re-order
df <- df[order(pID, dayoftest)]

# Filter out patients that were only symptomatic or never had symptoms (for now)
df <- df[, any_symp := any(symptomatic == TRUE),by = "pID"][any_symp == TRUE]
df <- df[, any_asym := any(symptomatic == FALSE), by = "pID"][any_asym == TRUE]

# Get times of last asymptomatic test and first symptomatic test
first_symptomatic <- df[symptomatic == TRUE, min(dayoftest), by = "pID"]
last_asymptomatic <- df[symptomatic == FALSE, max(dayoftest), by = "pID"]

# Assign new ID from 1 to number of patients left
df[, new_id := .GRP, by = "pID"]

# Put in a list
dat <- list()
dat$N <- nrow(df)
dat$P <- length(unique(df$pID))
dat$day_of_test <- df$dayoftest
dat$symptomatic <- df$symptomatic
dat$test_result <- df$test_result
dat$incub_alpha <- incub_alpha
dat$incub_beta <- incub_beta
dat$patient_ID <- df$new_id
dat$time_first_symptom <- first_symptomatic$V1
dat$time_last_asym <- last_asymptomatic$V1
dat$T_e <- df[, actual_te[1], by = "pID"][,V1]
dat$time_first_test <- df[, min(dayoftest), by = "pID"][,V1]

options(mc.cores = parallel::detectCores())
mod <- rstan::stan_model("~/repos/pcr-profile/pcr_model.stan")

fit <- rstan::sampling(mod, chains = 1, 
                       iter = 1000, 
                       data = dat,
                       control = list(max_treedepth = 10))
                       # algorithm = "Fixed_param")#,
res <- rstan::extract(fit)

p1 <- data.frame(med = apply(res$T_e, 2, median),
           LQ = apply(res$T_e, 2, quantile, prob = 0.025),
           UQ = apply(res$T_e, 2, quantile, prob = 0.975),
           actual = dat$T_e) %>%
 ggplot2::ggplot(ggplot2::aes(x = actual, y = med, ymin = LQ, ymax = UQ)) +
 ggplot2::geom_point() + 
 ggplot2::geom_errorbar() +
 ggplot2::geom_abline(intercept = 0, slope = 1, lty = 2) +
 ggplot2::labs(x = "True value", y = "Estimated value") +
 ggplot2::theme_bw()

p1

