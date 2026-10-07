# 1. Import the original file
commodities <- read.csv2(
  file.choose(),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Retain dates and the three price series
commodities <- commodities[, c("Date", "Gold", "Silver", "Platinum")]
commodities$Date <- as.Date(commodities$Date, format = "%Y/%m/%d")
commodities <- commodities[order(commodities$Date), ]

price_names <- c("Gold", "Silver", "Platinum")

stopifnot(
  !anyNA(commodities),
  !anyDuplicated(commodities$Date),
  all(as.matrix(commodities[, price_names]) > 0)
)

# 2. Identify consecutive calendar months
month_number <- as.integer(format(commodities$Date, "%Y")) * 12 +
  as.integer(format(commodities$Date, "%m"))

consecutive <- diff(month_number) == 1

# Show intervals spanning missing months
gap_intervals <- data.frame(
  Previous_date = head(commodities$Date, -1),
  Current_date = tail(commodities$Date, -1),
  Months_apart = diff(month_number)
)

print(gap_intervals[!consecutive, ])

# 3. Calculate gross growth factors, dated at the ending month
prices <- as.matrix(commodities[, price_names])

growth <- data.frame(
  Date = commodities$Date[-1],
  prices[-1, ] / prices[-nrow(prices), ],
  check.names = FALSE
)

# Keep only genuine one-month intervals for analysis
growth_monthly <- growth[consecutive, ]
rownames(growth_monthly) <- NULL

cat("Monthly observations:", nrow(growth_monthly), "\n")

# 4. Descriptive statistics
describe_growth <- function(x) {
  c(
    N = length(x),
    Mean = mean(x),
    SD = sd(x),
    Minimum = min(x),
    Q1 = unname(quantile(x, 0.25)),
    Median = median(x),
    Q3 = unname(quantile(x, 0.75)),
    Maximum = max(x)
  )
}

descriptive_table <- data.frame(
  Commodity = price_names,
  t(vapply(
    growth_monthly[, price_names],
    describe_growth,
    numeric(8)
  )),
  row.names = NULL
)

print(descriptive_table, digits = 6)

# 5. Time plots
# Use NA at gaps so lines do not connect across missing months
growth_plot <- growth
growth_plot[!consecutive, price_names] <- NA_real_

old_par <- par(no.readonly = TRUE)
par(mfrow = c(3, 1), mar = c(3, 4, 2, 1))

plot_colours <- c("darkgoldenrod3", "steelblue4", "purple4")

for (j in seq_along(price_names)) {
  commodity <- price_names[j]
  
  plot(
    growth_plot$Date,
    growth_plot[[commodity]],
    type = "l",
    col = plot_colours[j],
    xlab = "Date",
    ylab = "Growth factor",
    main = paste(commodity, "monthly growth factors")
  )
  
  abline(h = 1, lty = 2, col = "grey60")
}

par(old_par)







# 6. Install and load the marginal-fitting package
if (!requireNamespace("fitdistrplus", quietly = TRUE)) {
  install.packages("fitdistrplus")
}

library(fitdistrplus)

# Confirm the prepared data
stopifnot(
  nrow(growth_monthly) == 298,
  all(vapply(growth_monthly[price_names], is.numeric, logical(1))),
  !anyNA(growth_monthly[price_names])
)

# 7. Fit candidate marginal distributions
candidate_distributions <- c("norm", "lnorm", "gamma", "weibull")

marginal_fits <- setNames(vector("list", length(price_names)), price_names)
comparison_rows <- list()

for (commodity in price_names) {
  
  x <- growth_monthly[[commodity]]
  commodity_fits <- list()
  
  for (distribution in candidate_distributions) {
    
    fitted_model <- tryCatch(
      fitdist(x, distribution, method = "mle"),
      error = function(e) {
        message(
          commodity, " / ", distribution, ": ",
          conditionMessage(e)
        )
        NULL
      }
    )
    
    if (!is.null(fitted_model)) {
      
      commodity_fits[[distribution]] <- fitted_model
      
      comparison_rows[[length(comparison_rows) + 1L]] <- data.frame(
        Commodity = commodity,
        Distribution = distribution,
        LogLik = fitted_model$loglik,
        AIC = fitted_model$aic,
        BIC = fitted_model$bic,
        Convergence = fitted_model$convergence
      )
    }
  }
  
  marginal_fits[[commodity]] <- commodity_fits
}

marginal_comparison <- do.call(rbind, comparison_rows)
rownames(marginal_comparison) <- NULL

# Print each commodity's comparison, ordered by AIC
for (commodity in price_names) {
  
  cat("\n---", commodity, ": marginal comparison ---\n")
  
  comparison <- subset(
    marginal_comparison,
    Commodity == commodity
  )
  
  comparison <- comparison[
    order(comparison$AIC),
    names(comparison),
    drop = FALSE
  ]
  
  print(comparison, row.names = FALSE, digits = 6)
}

# 8. Select the lowest-AIC converged fit
best_marginals <- setNames(
  vector("list", length(price_names)),
  price_names
)

selection_rows <- list()
parameter_rows <- list()

for (commodity in price_names) {
  
  comparison <- subset(
    marginal_comparison,
    Commodity == commodity & Convergence == 0
  )
  
  if (nrow(comparison) == 0) {
    stop("No converged marginal fit for ", commodity)
  }
  
  aic_choice <- comparison$Distribution[which.min(comparison$AIC)]
  bic_choice <- comparison$Distribution[which.min(comparison$BIC)]
  
  best_fit <- marginal_fits[[commodity]][[aic_choice]]
  best_marginals[[commodity]] <- best_fit
  
  selection_rows[[commodity]] <- data.frame(
    Commodity = commodity,
    Best_AIC = aic_choice,
    Best_BIC = bic_choice
  )
  
  standard_errors <- best_fit$sd
  
  if (is.null(standard_errors)) {
    standard_errors <- rep(NA_real_, length(best_fit$estimate))
  }
  
  parameter_rows[[commodity]] <- data.frame(
    Commodity = commodity,
    Distribution = aic_choice,
    Parameter = names(best_fit$estimate),
    Estimate = unname(best_fit$estimate),
    Standard_Error = unname(standard_errors)
  )
}

marginal_selection <- do.call(rbind, selection_rows)
marginal_parameters <- do.call(rbind, parameter_rows)

rownames(marginal_selection) <- NULL
rownames(marginal_parameters) <- NULL

cat("\n--- Selected marginal distributions ---\n")
print(marginal_selection, row.names = FALSE)

cat("\n--- Selected marginal parameter estimates ---\n")
print(marginal_parameters, row.names = FALSE, digits = 6)

# 9. Four diagnostic plots for each selected marginal
# Run separately to view each commodity's four-panel figure.

plot(best_marginals[["Gold"]])
title("Gold: selected marginal diagnostics", outer = TRUE, line = -1)

plot(best_marginals[["Silver"]])
title("Silver: selected marginal diagnostics", outer = TRUE, line = -1)

plot(best_marginals[["Platinum"]])
title("Platinum: selected marginal diagnostics", outer = TRUE, line = -1)







# 10. Transform growth factors using the selected marginal CDFs

uniform_data <- data.frame(Date = growth_monthly$Date)

for (commodity in price_names) {
  
  fitted_model <- best_marginals[[commodity]]
  
  cdf_function <- get(
    paste0("p", fitted_model$distname),
    mode = "function"
  )
  
  uniform_data[[commodity]] <- do.call(
    cdf_function,
    c(
      list(q = growth_monthly[[commodity]]),
      as.list(fitted_model$estimate)
    )
  )
}

# Check the transformed values
uniform_checks <- do.call(
  rbind,
  lapply(price_names, function(commodity) {
    
    u <- uniform_data[[commodity]]
    
    data.frame(
      Commodity = commodity,
      N = length(u),
      Minimum = min(u),
      Maximum = max(u),
      Missing = sum(is.na(u)),
      Boundary_values = sum(u <= 0 | u >= 1)
    )
  })
)

cat("\n--- Fitted CDF checks ---\n")
print(uniform_checks, row.names = FALSE, digits = 8)

stopifnot(
  all(vapply(
    uniform_data[price_names],
    function(u) all(is.finite(u) & u > 0 & u < 1),
    logical(1)
  ))
)

# 11. Define the three commodity pairs

commodity_pairs <- list(
  Gold_Silver = c("Gold", "Silver"),
  Gold_Platinum = c("Gold", "Platinum"),
  Silver_Platinum = c("Silver", "Platinum")
)

# Store the paired CDF values for subsequent copula fitting
copula_data <- lapply(
  commodity_pairs,
  function(pair) as.matrix(uniform_data[pair])
)

# 12. Scatterplots of the transformed variables
# Run each call separately to view each figure.

plot_uniform_pair <- function(pair_name) {
  
  pair <- commodity_pairs[[pair_name]]
  uv <- copula_data[[pair_name]]
  
  plot(
    x = uv[, 1L],
    y = uv[, 2L],
    pch = 16,
    cex = 0.7,
    col = adjustcolor("steelblue4", alpha.f = 0.6),
    xlim = c(0, 1),
    ylim = c(0, 1),
    asp = 1,
    xlab = paste0("u1: ", pair[1L], " fitted CDF"),
    ylab = paste0("u2: ", pair[2L], " fitted CDF"),
    main = paste(pair, collapse = "–")
  )
}

plot_uniform_pair("Gold_Silver")
plot_uniform_pair("Gold_Platinum")
plot_uniform_pair("Silver_Platinum")

# 13. Empirical Kendall's tau for each pair

kendall_results <- do.call(
  rbind,
  lapply(names(commodity_pairs), function(pair_name) {
    
    pair <- commodity_pairs[[pair_name]]
    
    data.frame(
      Pair = paste(pair, collapse = "–"),
      N = nrow(growth_monthly),
      Kendall_tau = cor(
        growth_monthly[[pair[1L]]],
        growth_monthly[[pair[2L]]],
        method = "kendall"
      )
    )
  })
)

rownames(kendall_results) <- NULL

cat("\n--- Empirical Kendall's tau ---\n")
print(kendall_results, row.names = FALSE, digits = 6)








# 14. Load the copula package
if (!requireNamespace("copula", quietly = TRUE)) {
  install.packages("copula")
}

library(copula)

# Adapted from the supplied copula examples and package documentation:
# compare all required families for each commodity pair,
# retaining convergence checks and fitted parameters.

copula_fits <- list()
comparison_rows <- list()

# 15. Fit candidate copulas
for (pair_name in names(copula_data)) {
  
  cat("\nFitting:", pair_name, "\n")
  
  uv <- copula_data[[pair_name]]
  n <- nrow(uv)
  
  empirical_tau <- cor(
    uv[, 1L], uv[, 2L],
    method = "kendall"
  )
  
  rho_start <- sin(pi * empirical_tau / 2)
  
  specifications <- list(
    Clayton = list(
      model = claytonCopula(dim = 2),
      start = 2 * empirical_tau / (1 - empirical_tau),
      lower = 1e-6,
      upper = Inf
    ),
    Gumbel = list(
      model = gumbelCopula(dim = 2),
      start = 1 / (1 - empirical_tau),
      lower = 1 + 1e-6,
      upper = Inf
    ),
    Frank = list(
      model = frankCopula(dim = 2),
      start = iTau(frankCopula(), empirical_tau),
      lower = 1e-6,
      upper = Inf
    ),
    Gaussian = list(
      model = normalCopula(dim = 2, dispstr = "un"),
      start = rho_start,
      lower = -0.999,
      upper = 0.999
    ),
    Student_t = list(
      model = tCopula(
        dim = 2, dispstr = "un",
        df = 4, df.fixed = FALSE
      ),
      start = c(rho_start, 4),
      lower = c(-0.999, 0.1),
      upper = c(0.999, 200)
    )
  )
  
  # Independence has no estimated copula parameters
  pair_fits <- list(Independence = indepCopula(dim = 2))
  
  comparison_rows[[length(comparison_rows) + 1L]] <- data.frame(
    Pair = pair_name,
    Family = "Independence",
    Parameters = 0L,
    LogLik = 0,
    AIC = 0,
    BIC = 0,
    Convergence = 0L
  )
  
  for (family in names(specifications)) {
    
    spec <- specifications[[family]]
    
    fitted_model <- tryCatch(
      fitCopula(
        copula = spec$model,
        data = uv,
        method = "ml",
        start = spec$start,
        lower = spec$lower,
        upper = spec$upper,
        optim.method = "L-BFGS-B",
        optim.control = list(maxit = 2000),
        estimate.variance = FALSE
      ),
      error = function(e) {
        message(pair_name, " / ", family, ": ", conditionMessage(e))
        NULL
      }
    )
    
    if (is.null(fitted_model)) next
    
    pair_fits[[family]] <- fitted_model
    
    log_likelihood <- as.numeric(logLik(fitted_model))
    k <- length(coef(fitted_model))
    convergence <- fitted_model@fitting.stats$convergence
    
    if (is.null(convergence)) convergence <- NA_integer_
    
    comparison_rows[[length(comparison_rows) + 1L]] <- data.frame(
      Pair = pair_name,
      Family = family,
      Parameters = k,
      LogLik = log_likelihood,
      AIC = -2 * log_likelihood + 2 * k,
      BIC = -2 * log_likelihood + log(n) * k,
      Convergence = convergence
    )
  }
  
  copula_fits[[pair_name]] <- pair_fits
}

copula_comparison <- do.call(rbind, comparison_rows)
rownames(copula_comparison) <- NULL

# 16. Print comparisons ordered by AIC
for (pair_name in names(copula_data)) {
  
  cat("\n---", pair_name, ": copula comparison ---\n")
  
  comparison <- subset(
    copula_comparison,
    Pair == pair_name
  )
  
  comparison <- comparison[
    order(comparison$AIC),
    names(comparison),
    drop = FALSE
  ]
  
  print(comparison, row.names = FALSE, digits = 6)
}

# 17. Select converged models and retain their copula objects
best_copulas <- list()
selection_rows <- list()

for (pair_name in names(copula_data)) {
  
  comparison <- subset(
    copula_comparison,
    Pair == pair_name &
      Convergence == 0 &
      is.finite(LogLik)
  )
  
  aic_choice <- comparison$Family[which.min(comparison$AIC)]
  bic_choice <- comparison$Family[which.min(comparison$BIC)]
  
  selected_fit <- copula_fits[[pair_name]][[aic_choice]]
  
  if (aic_choice == "Independence") {
    best_copulas[[pair_name]] <- selected_fit
  } else {
    best_copulas[[pair_name]] <- selected_fit@copula
  }
  
  selection_rows[[pair_name]] <- data.frame(
    Pair = pair_name,
    Best_AIC = aic_choice,
    Best_BIC = bic_choice
  )
}

copula_selection <- do.call(rbind, selection_rows)
rownames(copula_selection) <- NULL

cat("\n--- Selected copulas ---\n")
print(copula_selection, row.names = FALSE)

# 18. Print selected parameter estimates
for (pair_name in names(best_copulas)) {
  
  family <- copula_selection$Best_AIC[
    copula_selection$Pair == pair_name
  ]
  
  cat("\n---", pair_name, ":", family, "parameters ---\n")
  
  if (family == "Independence") {
    cat("No estimated parameters.\n")
  } else {
    print(coef(copula_fits[[pair_name]][[family]]))
  }
}








# 19. Plot fitted copula CDF and density surfaces

plot_copula_surfaces <- function(pair_name) {
  
  fitted_copula <- best_copulas[[pair_name]]
  pair <- commodity_pairs[[pair_name]]
  
  grid_values <- seq(0.01, 0.99, length.out = 60)
  
  grid_points <- as.matrix(expand.grid(
    u1 = grid_values,
    u2 = grid_values
  ))
  
  cdf_values <- matrix(
    pCopula(grid_points, fitted_copula),
    nrow = length(grid_values)
  )
  
  density_values <- matrix(
    dCopula(grid_points, fitted_copula),
    nrow = length(grid_values)
  )
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  par(mfrow = c(1, 2), mar = c(3, 3, 3, 1))
  
  persp(
    grid_values, grid_values, cdf_values,
    theta = 35, phi = 25,
    col = "lightblue",
    border = NA,
    ticktype = "detailed",
    xlab = "u1",
    ylab = "u2",
    zlab = "C(u1, u2)",
    main = paste(paste(pair, collapse = "–"), "CDF")
  )
  
  persp(
    grid_values, grid_values, density_values,
    theta = 35, phi = 25,
    col = "wheat",
    border = NA,
    ticktype = "detailed",
    xlab = "u1",
    ylab = "u2",
    zlab = "c(u1, u2)",
    main = paste(paste(pair, collapse = "–"), "density")
  )
}

# Run separately to view each figure
plot_copula_surfaces("Gold_Silver")
plot_copula_surfaces("Gold_Platinum")
plot_copula_surfaces("Silver_Platinum")


# 20. Simulate 10,000 observations from each fitted copula
# Uses the estimated parameters stored in best_copulas.

set.seed(6823)

simulated_uniform <- lapply(
  best_copulas,
  function(fitted_copula) rCopula(10000, fitted_copula)
)


# 21. Transform simulations back to commodity growth factors

marginal_quantile <- function(u, commodity) {
  
  fitted_marginal <- best_marginals[[commodity]]
  
  quantile_function <- get(
    paste0("q", fitted_marginal$distname),
    mode = "function"
  )
  
  do.call(
    quantile_function,
    c(
      list(p = u),
      as.list(fitted_marginal$estimate)
    )
  )
}

simulated_growth <- lapply(
  names(commodity_pairs),
  function(pair_name) {
    
    pair <- commodity_pairs[[pair_name]]
    uv <- simulated_uniform[[pair_name]]
    
    simulated_pair <- data.frame(
      marginal_quantile(uv[, 1L], pair[1L]),
      marginal_quantile(uv[, 2L], pair[2L])
    )
    
    names(simulated_pair) <- pair
    simulated_pair
  }
)

names(simulated_growth) <- names(commodity_pairs)


# 22. Compare observed and simulated scatterplots

plot_simulation_comparison <- function(pair_name) {
  
  pair <- commodity_pairs[[pair_name]]
  observed_uv <- copula_data[[pair_name]]
  simulated_uv <- simulated_uniform[[pair_name]]
  simulated_pair <- simulated_growth[[pair_name]]
  
  x_observed <- growth_monthly[[pair[1L]]]
  y_observed <- growth_monthly[[pair[2L]]]
  
  x_limits <- range(x_observed, simulated_pair[[pair[1L]]])
  y_limits <- range(y_observed, simulated_pair[[pair[2L]]])
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  
  plot(
    observed_uv[, 1L], observed_uv[, 2L],
    pch = 16, cex = 0.6,
    col = adjustcolor("steelblue4", alpha.f = 0.6),
    xlim = c(0, 1), ylim = c(0, 1),
    xlab = paste(pair[1L], "fitted CDF"),
    ylab = paste(pair[2L], "fitted CDF"),
    main = "Observed: uniform scale"
  )
  
  plot(
    simulated_uv[, 1L], simulated_uv[, 2L],
    pch = 16, cex = 0.3,
    col = adjustcolor("firebrick3", alpha.f = 0.15),
    xlim = c(0, 1), ylim = c(0, 1),
    xlab = paste(pair[1L], "uniform value"),
    ylab = paste(pair[2L], "uniform value"),
    main = "Simulated: 10,000 observations"
  )
  
  plot(
    x_observed, y_observed,
    pch = 16, cex = 0.6,
    col = adjustcolor("steelblue4", alpha.f = 0.6),
    xlim = x_limits, ylim = y_limits,
    xlab = paste(pair[1L], "growth factor"),
    ylab = paste(pair[2L], "growth factor"),
    main = "Observed: growth scale"
  )
  
  plot(
    simulated_pair[[pair[1L]]],
    simulated_pair[[pair[2L]]],
    pch = 16, cex = 0.3,
    col = adjustcolor("firebrick3", alpha.f = 0.15),
    xlim = x_limits, ylim = y_limits,
    xlab = paste(pair[1L], "growth factor"),
    ylab = paste(pair[2L], "growth factor"),
    main = "Simulated: growth scale"
  )
}

# Run separately to view each figure
plot_simulation_comparison("Gold_Silver")
plot_simulation_comparison("Gold_Platinum")
plot_simulation_comparison("Silver_Platinum")


# 23. Check simulated dependence against fitted dependence

simulation_checks <- do.call(
  rbind,
  lapply(names(best_copulas), function(pair_name) {
    
    uv <- simulated_uniform[[pair_name]]
    observed_uv <- copula_data[[pair_name]]
    
    data.frame(
      Pair = pair_name,
      Simulations = nrow(uv),
      Observed_tau = cor(
        observed_uv[, 1L], observed_uv[, 2L],
        method = "kendall"
      ),
      Fitted_tau = as.numeric(tau(best_copulas[[pair_name]])),
      Simulated_tau = cor(
        uv[, 1L], uv[, 2L],
        method = "kendall"
      ),
      Mean_u1 = mean(uv[, 1L]),
      Mean_u2 = mean(uv[, 2L])
    )
  })
)

rownames(simulation_checks) <- NULL

cat("\n--- Simulation checks ---\n")
print(simulation_checks, row.names = FALSE, digits = 6)




# 24. Helper: selected marginal CDF
marginal_cdf <- function(x, commodity) {
  
  fitted_marginal <- best_marginals[[commodity]]
  
  cdf_function <- get(
    paste0("p", fitted_marginal$distname),
    mode = "function"
  )
  
  do.call(
    cdf_function,
    c(
      list(q = x),
      as.list(fitted_marginal$estimate)
    )
  )
}


# 25. Joint probabilities at growth-factor thresholds
# 0.95 = a 5% decrease; 1.00 = no change; 1.05 = a 5% increase.
#
# Both below: C(F1(x), F2(y))
# Both above: 1 - F1(x) - F2(y) + C(F1(x), F2(y))

threshold_scenarios <- data.frame(
  Scenario = c(
    "Both decrease by at least 5%",
    "Both have non-positive growth",
    "Both increase by more than 5%"
  ),
  Threshold_1 = c(0.95, 1.00, 1.05),
  Threshold_2 = c(0.95, 1.00, 1.05),
  Direction = c("Below", "Below", "Above")
)

probability_rows <- list()

for (pair_name in names(best_copulas)) {
  
  pair <- commodity_pairs[[pair_name]]
  fitted_copula <- best_copulas[[pair_name]]
  simulated_pair <- simulated_growth[[pair_name]]
  
  x_observed <- growth_monthly[[pair[1L]]]
  y_observed <- growth_monthly[[pair[2L]]]
  
  for (j in seq_len(nrow(threshold_scenarios))) {
    
    scenario <- threshold_scenarios[j, ]
    x <- scenario$Threshold_1
    y <- scenario$Threshold_2
    
    u <- marginal_cdf(x, pair[1L])
    v <- marginal_cdf(y, pair[2L])
    
    joint_cdf <- as.numeric(
      pCopula(matrix(c(u, v), nrow = 1), fitted_copula)
    )
    
    if (scenario$Direction == "Below") {
      
      model_probability <- joint_cdf
      independence_probability <- u * v
      
      empirical_probability <- mean(
        x_observed <= x & y_observed <= y
      )
      
      simulated_probability <- mean(
        simulated_pair[[pair[1L]]] <= x &
          simulated_pair[[pair[2L]]] <= y
      )
      
    } else {
      
      model_probability <- 1 - u - v + joint_cdf
      independence_probability <- (1 - u) * (1 - v)
      
      empirical_probability <- mean(
        x_observed > x & y_observed > y
      )
      
      simulated_probability <- mean(
        simulated_pair[[pair[1L]]] > x &
          simulated_pair[[pair[2L]]] > y
      )
    }
    
    probability_rows[[length(probability_rows) + 1L]] <- data.frame(
      Pair = pair_name,
      Scenario = scenario$Scenario,
      Threshold_1 = x,
      Threshold_2 = y,
      Marginal_CDF_1 = u,
      Marginal_CDF_2 = v,
      Joint_CDF = joint_cdf,
      Model_Probability = model_probability,
      Independence_Probability = independence_probability,
      Empirical_Probability = empirical_probability,
      Simulated_Probability = simulated_probability
    )
  }
}

joint_probabilities <- do.call(rbind, probability_rows)
rownames(joint_probabilities) <- NULL

for (pair_name in names(best_copulas)) {
  
  cat("\n---", pair_name, ": joint probabilities ---\n")
  
  print(
    subset(joint_probabilities, Pair == pair_name),
    row.names = FALSE,
    digits = 6
  )
}


# 26. Exact asymptotic tail-dependence coefficients
# lambda() returns lower and upper coefficients.

tail_dependence <- do.call(
  rbind,
  lapply(names(best_copulas), function(pair_name) {
    
    coefficients <- lambda(best_copulas[[pair_name]])
    
    data.frame(
      Pair = pair_name,
      Family = copula_selection$Best_AIC[
        copula_selection$Pair == pair_name
      ],
      Lower_Tail = unname(coefficients["lower"]),
      Upper_Tail = unname(coefficients["upper"])
    )
  })
)

rownames(tail_dependence) <- NULL

cat("\n--- Asymptotic tail dependence ---\n")
print(tail_dependence, row.names = FALSE, digits = 6)


# 27. Finite-threshold tail probabilities
# These are conditional probabilities at chosen quantiles,
# separate from the limiting coefficients above.
#
# Lower: P(U1 <= q | U2 <= q) = C(q, q) / q
# Upper: P(U1 > 1-q | U2 > 1-q)
#        = [2q - 1 + C(1-q, 1-q)] / q

tail_levels <- c(0.10, 0.05, 0.01)
finite_tail_rows <- list()

for (pair_name in names(best_copulas)) {
  
  fitted_copula <- best_copulas[[pair_name]]
  
  for (q in tail_levels) {
    
    lower_joint <- as.numeric(
      pCopula(matrix(c(q, q), nrow = 1), fitted_copula)
    )
    
    upper_joint <- 2 * q - 1 + as.numeric(
      pCopula(
        matrix(c(1 - q, 1 - q), nrow = 1),
        fitted_copula
      )
    )
    
    finite_tail_rows[[length(finite_tail_rows) + 1L]] <-
      data.frame(
        Pair = pair_name,
        Tail_Probability = q,
        Lower_Joint = lower_joint,
        Upper_Joint = upper_joint,
        Lower_Conditional = lower_joint / q,
        Upper_Conditional = upper_joint / q,
        Independence_Conditional = q
      )
  }
}

finite_tail_probabilities <- do.call(rbind, finite_tail_rows)
rownames(finite_tail_probabilities) <- NULL

cat("\n--- Finite-threshold tail probabilities ---\n")
print(finite_tail_probabilities, row.names = FALSE, digits = 6)

# Check probability bounds
stopifnot(
  all(is.finite(joint_probabilities$Model_Probability)),
  all(joint_probabilities$Model_Probability >= -1e-10),
  all(joint_probabilities$Model_Probability <= 1 + 1e-10),
  !anyNA(tail_dependence),
  all(finite_tail_probabilities$Lower_Conditional >= -1e-10),
  all(finite_tail_probabilities$Lower_Conditional <= 1 + 1e-10),
  all(finite_tail_probabilities$Upper_Conditional >= -1e-10),
  all(finite_tail_probabilities$Upper_Conditional <= 1 + 1e-10)
)





# 28. Create an output folder
output_dir <- file.path(getwd(), "Assignment_5_Output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 29. Export tables
tables_to_save <- list(
  "01_descriptive_statistics" = descriptive_table,
  "02_marginal_comparison" = marginal_comparison,
  "03_marginal_selection" = marginal_selection,
  "04_marginal_parameters" = marginal_parameters,
  "05_uniform_checks" = uniform_checks,
  "06_kendall_tau" = kendall_results,
  "07_copula_comparison" = copula_comparison,
  "08_copula_selection" = copula_selection,
  "09_simulation_checks" = simulation_checks,
  "10_joint_probabilities" = joint_probabilities,
  "11_tail_dependence" = tail_dependence,
  "12_finite_tail_probabilities" = finite_tail_probabilities
)

for (table_name in names(tables_to_save)) {
  write.csv(
    tables_to_save[[table_name]],
    file.path(output_dir, paste0(table_name, ".csv")),
    row.names = FALSE
  )
}

# 30. Save fitted models and data for reuse
saveRDS(
  list(
    growth_monthly = growth_monthly,
    best_marginals = best_marginals,
    marginal_fits = marginal_fits,
    copula_data = copula_data,
    copula_fits = copula_fits,
    best_copulas = best_copulas,
    simulated_uniform = simulated_uniform,
    simulated_growth = simulated_growth,
    tables = tables_to_save
  ),
  file.path(output_dir, "Assignment_5_results.rds")
)

# 31. Save package versions and citation information
capture.output(
  sessionInfo(),
  file = file.path(output_dir, "session_info.txt")
)

capture.output(
  citation("copula"),
  citation("fitdistrplus"),
  file = file.path(output_dir, "package_references.txt")
)

cat("\nOutput folder:\n", normalizePath(output_dir), "\n")
print(list.files(output_dir))


# 32. Create a figures folder
figure_dir <- file.path(output_dir, "Figures")
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

# Helper: save a plot and close its graphics device
save_figure <- function(filename, draw, width = 8, height = 6) {
  
  png(
    filename = file.path(figure_dir, filename),
    width = width,
    height = height,
    units = "in",
    res = 300,
    bg = "white"
  )
  
  on.exit(dev.off(), add = TRUE)
  draw()
}


# 33. Save time plots
save_figure(
  "Figure_01_Commodity_Growth_Time_Plots.png",
  draw = function() {
    
    par(mfrow = c(3, 1), mar = c(4, 4, 3, 1))
    
    plot_colours <- c("darkgoldenrod3", "steelblue4", "purple4")
    
    for (j in seq_along(price_names)) {
      
      commodity <- price_names[j]
      
      plot(
        growth_plot$Date,
        growth_plot[[commodity]],
        type = "l",
        col = plot_colours[j],
        xlab = "Date",
        ylab = "Growth factor",
        main = paste(commodity, "monthly growth factors")
      )
      
      abline(h = 1, lty = 2, col = "grey60")
    }
  },
  width = 8,
  height = 9
)


# 34. Save marginal diagnostic plots
for (j in seq_along(price_names)) {
  
  commodity <- price_names[j]
  
  save_figure(
    sprintf(
      "Figure_%02d_%s_Marginal_Diagnostics.png",
      j + 1L, commodity
    ),
    draw = function() {
      plot(best_marginals[[commodity]])
    },
    width = 8,
    height = 7
  )
}


# 35. Save transformed scatterplots
pair_names <- names(commodity_pairs)

for (j in seq_along(pair_names)) {
  
  pair_name <- pair_names[j]
  
  save_figure(
    sprintf(
      "Figure_%02d_%s_Transformed_Scatterplot.png",
      j + 4L, pair_name
    ),
    draw = function() {
      par(mar = c(5, 5, 4, 2))
      plot_uniform_pair(pair_name)
    },
    width = 7,
    height = 7
  )
}


# 36. Save copula CDF and density surfaces
for (j in seq_along(pair_names)) {
  
  pair_name <- pair_names[j]
  
  save_figure(
    sprintf(
      "Figure_%02d_%s_Copula_Surfaces.png",
      j + 7L, pair_name
    ),
    draw = function() {
      plot_copula_surfaces(pair_name)
    },
    width = 10,
    height = 5
  )
}


# 37. Save observed versus simulated comparisons
for (j in seq_along(pair_names)) {
  
  pair_name <- pair_names[j]
  
  save_figure(
    sprintf(
      "Figure_%02d_%s_Simulation_Comparison.png",
      j + 10L, pair_name
    ),
    draw = function() {
      plot_simulation_comparison(pair_name)
    },
    width = 9,
    height = 8
  )
}


# 38. Check the saved figures
figure_files <- list.files(
  figure_dir,
  pattern = "\\.png$",
  full.names = TRUE
)

stopifnot(
  length(figure_files) == 13L,
  all(file.info(figure_files)$size > 0)
)

cat("\nSaved figures:", length(figure_files), "\n")
cat("Figure folder:\n", normalizePath(figure_dir), "\n")
print(basename(figure_files))


