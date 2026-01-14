# Replication.R  (Python -> R translation of Replication.py)
# -----------------------------------------------------------
# Requirements:
#   - country_cleaned.dta
#   - country_sector_cleaned.dta
#
# Edit the DATA_DIR path below before running.
#
# Outputs:
#   - Several .csv tables and .png figures (mirroring the Python script names).
#
# Notes:
#   - This is a faithful translation of the logic in Replication.py.
#   - Some model tables differ slightly in formatting from Python/linearmodels.
#   - Robust/clustered SEs are computed where applicable.

# ---------------------------
# 0) Packages
# ---------------------------
suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(patchwork)
  library(moments)     # skewness/kurtosis
  library(plm)         # panel models
  library(sandwich)    # robust vcov
  library(lmtest)      # coeftest
  library(AER)         # ivreg
  library(tseries)     # adf.test
})

# ---------------------------
# 1) Paths (EDIT THIS)
# ---------------------------
DATA_DIR <- "/Users/raph/Library/Mobile Documents/com~apple~CloudDocs/M2/Panel Data/Project/Data replication/Data"
path_sector  <- file.path(DATA_DIR, "country_sector_cleaned.dta")
path_country <- file.path(DATA_DIR, "country_cleaned.dta")

# ---------------------------
# 2) Load data
# ---------------------------
df_sector  <- read_dta(path_sector)  %>% as.data.frame()
df_country <- read_dta(path_country) %>% as.data.frame()

# Helper: histogram + KDE + normal curve
plot_dist <- function(x, title, fill = "grey80") {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(ggplot() + ggtitle(paste(title, "(insufficient data)")))
  mu <- mean(x); sig <- sd(x)
  sk <- moments::skewness(x); ku <- moments::kurtosis(x)
  ggplot(data.frame(x = x), aes(x)) +
    geom_histogram(aes(y = after_stat(density)), bins = 40, fill = fill, color = NA, alpha = 0.6) +
    geom_density(linewidth = 1) +
    stat_function(fun = dnorm, args = list(mean = mu, sd = sig), linetype = "dashed") +
    labs(title = sprintf("%s\nSkew: %.2f, Kurtosis: %.2f", title, sk, ku),
         x = NULL, y = "Density")
}

# Helper: clustered SE coeftest for plm
plm_coeftest_cluster <- function(model, cluster = "group") {
  vc <- plm::vcovHC(model, type = "HC1", cluster = cluster)
  lmtest::coeftest(model, vcov = vc)
}

# Helper: Fisher-type ADF combination (country-by-country)
fisher_adf <- function(df, id, value, min_n = 11) {
  ids <- unique(df[[id]])
  pvals <- c()
  for (g in ids) {
    x <- df %>% filter(.data[[id]] == g) %>% pull(.data[[value]])
    x <- x[is.finite(x)]
    if (length(x) >= min_n) {
      # adf.test can fail on constant series; tryCatch
      pv <- tryCatch(tseries::adf.test(x)$p.value, error = function(e) NA_real_)
      if (is.finite(pv)) pvals <- c(pvals, pv)
    }
  }
  if (length(pvals) == 0) {
    return(list(n = 0, fisher = NA_real_, df = NA_integer_, p = NA_real_))
  }
  fisher <- -2 * sum(log(pvals))
  dfree  <- 2 * length(pvals)
  p      <- stats::pchisq(fisher, df = dfree, lower.tail = FALSE)
  list(n = length(pvals), fisher = fisher, df = dfree, p = p)
}

# ===========================================================
# Question 1
# ===========================================================
cat("\n====================\nQuestion 1\n====================\n")
n_countries <- n_distinct(df_country$country)
n_sectors   <- n_distinct(df_sector$industry)

df_sector <- df_sector %>%
  mutate(country_sector = paste0(as.character(country), "_", as.character(industry)))

n_individuals <- n_distinct(df_sector$country_sector)

t_min <- min(df_country$year, na.rm = TRUE)
t_max <- max(df_country$year, na.rm = TRUE)
max_t_span <- as.integer(t_max - t_min + 1)

obs_counts <- df_sector %>% count(country_sector, name = "n")
n_single_obs <- sum(obs_counts$n == 1)

df_sorted <- df_sector %>% arrange(country_sector, year) %>%
  group_by(country_sector) %>% mutate(year_diff = year - lag(year)) %>% ungroup()

consecutive_obs <- df_sorted %>%
  filter(year_diff == 1) %>% summarise(n = n_distinct(country_sector)) %>% pull(n)

cat(sprintf("Number of Countries: %d\n", n_countries))
cat(sprintf("Number of Sectors: %d\n", n_sectors))
cat(sprintf("Total Country-Sector individuals: %d\n", n_individuals))
cat(sprintf("Period: %d to %d (Max T Span = %d)\n", t_min, t_max, max_t_span))
cat(sprintf("Individuals with only 1 observation: %d\n", n_single_obs))
cat(sprintf("Individuals with ≥2 consecutive obs: %d\n", consecutive_obs))

# ===========================================================
# Question 4: Variance decomposition + between/within columns
# ===========================================================
cat("\n====================\nQuestion 4\n====================\n")
df <- df_country %>% as_tibble()

id_col <- "country"
numeric_cols <- df %>%
  select(where(is.numeric)) %>% names()
numeric_cols <- setdiff(numeric_cols, "year")

variance_results <- list()

for (col in numeric_cols) {
  temp <- df %>% select(all_of(c(id_col, col))) %>% filter(!is.na(.data[[col]]))
  if (nrow(temp) < 2) next

  pooled_var <- var(temp[[col]])
  if (!is.finite(pooled_var) || pooled_var == 0) next

  country_means <- temp %>% group_by(.data[[id_col]]) %>%
    mutate(mean_i = mean(.data[[col]], na.rm = TRUE)) %>% ungroup()

  between_var <- var(country_means$mean_i)
  within_val  <- country_means[[col]] - country_means$mean_i
  within_var  <- var(within_val)

  share_between <- 100 * between_var / pooled_var
  share_within  <- 100 * within_var  / pooled_var

  variance_results[[col]] <- data.frame(
    Variable = col,
    Pooled_Var = pooled_var,
    Between_Var = between_var,
    Within_Var = within_var,
    Between_Share_pct = share_between,
    Within_Share_pct = share_within
  )

  # add transformations to main df
  df[[paste0(col, "_between")]] <- df %>% group_by(.data[[id_col]]) %>%
    mutate(tmp = mean(.data[[col]], na.rm = TRUE)) %>% pull(tmp)

  df[[paste0(col, "_within")]] <- df[[col]] - df[[paste0(col, "_between")]]
}

summary_table <- bind_rows(variance_results)
print(summary_table %>% mutate(across(where(is.numeric), ~round(.x, 4))))
write.csv(summary_table, "Question4_Variance_Decomposition.csv", row.names = FALSE)

# ===========================================================
# Question 5: Distribution analysis (Within vs Between)
# ===========================================================
cat("\n====================\nQuestion 5\n====================\n")
dep_var <- "cotwo_total_per_cap"
exp_var <- "log_gdp_per_cap"

plot_data <- list()

for (v in c(dep_var, exp_var)) {
  temp <- df_country %>% select(country, year, all_of(v)) %>% filter(!is.na(.data[[v]]))
  b_mean <- temp %>% group_by(country) %>% mutate(b = mean(.data[[v]], na.rm = TRUE)) %>% pull(b)
  w_val  <- temp[[v]] - b_mean
  plot_data[[paste0(v, "_within")]]  <- w_val
  plot_data[[paste0(v, "_between")]] <- b_mean
}

p1 <- plot_dist(plot_data[[paste0(dep_var, "_within")]],  "Within-Variation: CO2 per Cap")
p2 <- plot_dist(plot_data[[paste0(dep_var, "_between")]], "Between-Variation: CO2 per Cap")
p3 <- plot_dist(plot_data[[paste0(exp_var, "_within")]],  "Within-Variation: Log GDP per Cap")
p4 <- plot_dist(plot_data[[paste0(exp_var, "_between")]], "Between-Variation: Log GDP per Cap")

g_q5 <- (p1 | p2) / (p3 | p4)
ggsave("Question5_Distribution_Within_Between.png", g_q5, width = 16, height = 12, dpi = 150)

# ===========================================================
# Question 6: First differences for all numeric vars
# ===========================================================
cat("\n====================\nQuestion 6\n====================\n")
df6 <- df_country %>% arrange(country, year)

numeric_vars <- df6 %>% select(where(is.numeric)) %>% names()
numeric_vars <- setdiff(numeric_vars, "year")

df6 <- df6 %>%
  group_by(country) %>%
  mutate(across(all_of(numeric_vars), ~ .x - lag(.x), .names = "{.col}_fd")) %>%
  ungroup()

# verification: show first obs of each country (should have NA diffs)
verification_cols <- c("country", "year", "log_gdp_per_cap", "log_gdp_per_cap_fd")
country_changes <- df6$country != dplyr::lag(df6$country)
cat("\nVerification of First Difference Logic (first obs per country):\n")
print(df6 %>% filter(country_changes) %>% select(any_of(verification_cols)) %>% head(10))

write.csv(df6, "Question_6_country_first_differences.csv", row.names = FALSE)

# ===========================================================
# Question 7: Within vs FD for y=cotwo_total_per_gdp, x=fin_str2
# ===========================================================
cat("\n====================\nQuestion 7\n====================\n")
dep_var <- "cotwo_total_per_gdp"
exp_var <- "fin_str2"

df7 <- df_country %>% arrange(country, year) %>%
  group_by(country) %>%
  mutate(dep_within = .data[[dep_var]] - mean(.data[[dep_var]], na.rm = TRUE),
         exp_within = .data[[exp_var]] - mean(.data[[exp_var]], na.rm = TRUE),
         dep_fd     = .data[[dep_var]] - lag(.data[[dep_var]]),
         exp_fd     = .data[[exp_var]] - lag(.data[[exp_var]])) %>%
  ungroup()

cfg <- list(
  list("dep_within", "Within: CO2/GDP (Dependent)", "grey80"),
  list("dep_fd",     "First-Difference: Δ(CO2/GDP)", "grey80"),
  list("exp_within", "Within: Equity Share (Explanatory)", "grey80"),
  list("exp_fd",     "First-Difference: Δ Equity Share", "grey80")
)

p_q7 <- map(cfg, \(z) plot_dist(df7[[z[[1]]]], z[[2]], fill = z[[3]]))
g_q7 <- (p_q7[[1]] | p_q7[[2]]) / (p_q7[[3]] | p_q7[[4]])
ggsave("Question7_Comparative_Distributions.png", g_q7, width = 16, height = 12, dpi = 150)

# ===========================================================
# Question 8: Joint relationship in FD space (scatter + lm)
# ===========================================================
cat("\n====================\nQuestion 8\n====================\n")
plot_df <- df7 %>% select(dep_fd, exp_fd) %>% filter(complete.cases(.))
corr_q8 <- cor(plot_df$exp_fd, plot_df$dep_fd)

p_q8 <- ggplot(plot_df, aes(x = exp_fd, y = dep_fd)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "lm", se = FALSE) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = sprintf("Pearson r = %.4f", corr_q8)) +
  labs(x = "Δ Equity Share (fin_str2)", y = "Δ CO2/GDP", title = "Question 8: FD Correlation") +
  theme_minimal()

ggsave("Question8_FD_Correlation.png", p_q8, width = 10, height = 8, dpi = 150)

# ===========================================================
# Question 9: Balanced panel + time component + TWFE transform
# ===========================================================
cat("\n====================\nQuestion 9\n====================\n")
year_counts <- df_country %>%
  group_by(year) %>%
  summarise(n_countries = n_distinct(country), .groups = "drop")

max_countries <- max(year_counts$n_countries)
balanced_years <- year_counts %>% filter(n_countries == max_countries) %>% pull(year)

country_counts <- df_country %>%
  filter(year %in% balanced_years) %>%
  group_by(country) %>%
  summarise(n = n(), .groups = "drop")

balanced_countries <- country_counts %>% filter(n == length(balanced_years)) %>% pull(country)

df_balanced <- df_country %>%
  filter(country %in% balanced_countries, year %in% balanced_years) %>%
  arrange(country, year)

dep_var <- "cotwo_total_per_gdp"
grand_mean <- mean(df_balanced[[dep_var]], na.rm = TRUE)

df_balanced <- df_balanced %>%
  group_by(country) %>% mutate(mu_i = mean(.data[[dep_var]], na.rm = TRUE)) %>% ungroup() %>%
  group_by(year)    %>% mutate(mu_t = mean(.data[[dep_var]], na.rm = TRUE)) %>% ungroup() %>%
  mutate(dep_twfe = .data[[dep_var]] - mu_i - mu_t + grand_mean)

time_component <- df_balanced %>%
  group_by(year) %>%
  summarise(mu_t = mean(mu_t, na.rm = TRUE), .groups = "drop") %>%
  mutate(time_component = -mu_t + grand_mean)

p_time <- ggplot(time_component, aes(x = year, y = time_component)) +
  geom_line() + geom_point() +
  labs(title = "Time-Specific Component: -mean_t + grand mean", x = "Year", y = "Value") +
  theme_minimal()

ggsave("Question9_Time_Component.png", p_time, width = 10, height = 6, dpi = 150)

cat("\nFirst 5 observations of the TWFE transformed variable:\n")
print(df_balanced %>% select(country, year, all_of(dep_var), dep_twfe) %>% head(5))

# ===========================================================
# Question 10: Boxplots of CO2/GDP by country ordered by within variance
# ===========================================================
cat("\n====================\nQuestion 10\n====================\n")
country_variances <- df_balanced %>%
  group_by(country) %>%
  summarise(within_var = var(cotwo_total_per_gdp, na.rm = TRUE), .groups = "drop") %>%
  arrange(within_var)

ordered_countries <- country_variances$country

p_q10 <- ggplot(df_balanced, aes(x = factor(country, levels = ordered_countries),
                                y = cotwo_total_per_gdp)) +
  geom_boxplot(outlier.alpha = 0.2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(title = "Question 10: CO2/GDP by Country (ordered by within-country variance)",
       x = "Country", y = "CO2 emissions per GDP")

ggsave("Question10_Boxplots_by_Variance.png", p_q10, width = 16, height = 8, dpi = 150)

# ===========================================================
# Question 11: Country-specific TWFE correlations + 'reg coef' definition
# ===========================================================
cat("\n====================\nQuestion 11\n====================\n")
y_col <- "cotwo_total_per_gdp"
x_col <- "log_gdp_per_cap"

df11 <- df_balanced

# TWFE transform both
for (pair in list(c(y_col, "y_twfe"), c(x_col, "x_twfe"))) {
  col <- pair[[1]]; target <- pair[[2]]
  mu_grand <- mean(df11[[col]], na.rm = TRUE)
  mu_i <- df11 %>% group_by(country) %>% mutate(mi = mean(.data[[col]], na.rm = TRUE)) %>% pull(mi)
  mu_t <- df11 %>% group_by(year) %>% mutate(mt = mean(.data[[col]], na.rm = TRUE)) %>% pull(mt)
  df11[[target]] <- df11[[col]] - mu_i - mu_t + mu_grand
}

summary_q11 <- df11 %>%
  group_by(country) %>%
  summarise(
    rho = cor(x_twfe, y_twfe, use = "complete.obs"),
    std_x = sd(x_twfe, na.rm = TRUE),
    std_y = sd(y_twfe, na.rm = TRUE),
    reg_coef = ifelse(std_y == 0, 0, rho * (std_x / std_y)),
    .groups = "drop"
  ) %>%
  arrange(desc(rho))

print(summary_q11 %>% mutate(across(where(is.numeric), ~round(.x, 4))))
write.csv(summary_q11, "Question11_TWFE_Correlations.csv", row.names = FALSE)

# ===========================================================
# Question 12: Unbalanced TWFE residuals via within + year dummies
# ===========================================================
cat("\n====================\nQuestion 12\n====================\n")
obs_counts <- df_country %>% count(country, name = "n")
countries_to_keep <- obs_counts %>% filter(n > 1) %>% pull(country)

df_unbal <- df_country %>% filter(country %in% countries_to_keep)

y_col <- "cotwo_total_per_gdp"
x_col <- "log_gdp_per_cap"

twfe_unbal <- function(data, col) {
  tmp <- data %>% select(country, year, all_of(col)) %>% filter(!is.na(.data[[col]]))
  tmp <- tmp %>% group_by(country) %>% mutate(within = .data[[col]] - mean(.data[[col]], na.rm = TRUE)) %>% ungroup()
  # year dummies regression (no intercept to match Python intent)
  fit <- lm(within ~ 0 + factor(year), data = tmp)
  tmp$resid <- resid(fit)
  tmp %>% select(country, year, resid)
}

y_res <- twfe_unbal(df_unbal, y_col) %>% rename(y_twfe_unbal = resid)
x_res <- twfe_unbal(df_unbal, x_col) %>% rename(x_twfe_unbal = resid)

df_unbal_final <- df_unbal %>%
  left_join(y_res, by = c("country", "year")) %>%
  left_join(x_res, by = c("country", "year"))

cat("\nQuestion 12: Descriptive stats for unbalanced TWFE residuals\n")
print(summary(df_unbal_final %>% select(y_twfe_unbal, x_twfe_unbal)))

# ===========================================================
# Question 13: Unbalanced between/within/TWFE distribution comparison (2x3)
# ===========================================================
cat("\n====================\nQuestion 13\n====================\n")
df_unbal2 <- df_unbal %>% arrange(country, year)

y_col <- "cotwo_total_per_gdp"
x_col <- "log_gdp_per_cap"

df_unbal2 <- df_unbal2 %>%
  group_by(country) %>%
  mutate(
    y_between = mean(.data[[y_col]], na.rm = TRUE),
    y_within  = .data[[y_col]] - y_between,
    x_between = mean(.data[[x_col]], na.rm = TRUE),
    x_within  = .data[[x_col]] - x_between
  ) %>% ungroup()

# TWFE residuals (unbalanced): within then year dummies
y_twfe <- twfe_unbal(df_unbal2, y_col) %>% rename(y_twfe = resid)
x_twfe <- twfe_unbal(df_unbal2, x_col) %>% rename(x_twfe = resid)

df_unbal2 <- df_unbal2 %>%
  left_join(y_twfe, by = c("country", "year")) %>%
  left_join(x_twfe, by = c("country", "year"))

p13 <- list(
  plot_dist(df_unbal2$y_between, "Between: CO2/GDP", "grey80"),
  plot_dist(df_unbal2$y_within,  "Within: CO2/GDP",  "grey80"),
  plot_dist(df_unbal2$y_twfe,    "TWFE: CO2/GDP",    "grey80"),
  plot_dist(df_unbal2$x_between, "Between: Log GDP per Cap", "grey80"),
  plot_dist(df_unbal2$x_within,  "Within: Log GDP per Cap",  "grey80"),
  plot_dist(df_unbal2$x_twfe,    "TWFE: Log GDP per Cap",    "grey80")
)

g13 <- (p13[[1]] | p13[[2]] | p13[[3]]) / (p13[[4]] | p13[[5]] | p13[[6]])
ggsave("Question13_Distribution_Comparison.png", g13, width = 20, height = 12, dpi = 150)

# ===========================================================
# Question 14: Boxplots of transformations (select 20 countries)
# ===========================================================
cat("\n====================\nQuestion 14\n====================\n")
df_work <- df_country %>%
  filter(country %in% countries_to_keep) %>%
  arrange(country, year)

y_col <- "cotwo_total_per_gdp"
x_col <- "fin_str2"

df_plot <- df_work %>%
  group_by(country) %>%
  mutate(
    y_between = mean(.data[[y_col]], na.rm = TRUE),
    y_within  = .data[[y_col]] - y_between,
    y_fd      = .data[[y_col]] - lag(.data[[y_col]]),
    x_between = mean(.data[[x_col]], na.rm = TRUE),
    x_within  = .data[[x_col]] - x_between,
    x_fd      = .data[[x_col]] - lag(.data[[x_col]])
  ) %>% ungroup()

# TWFE residuals for each series
y_twfe2 <- twfe_unbal(df_plot, y_col) %>% rename(y_twfe = resid)
x_twfe2 <- twfe_unbal(df_plot, x_col) %>% rename(x_twfe = resid)

df_plot <- df_plot %>%
  left_join(y_twfe2, by = c("country", "year")) %>%
  left_join(x_twfe2, by = c("country", "year"))

selected_countries <- sort(unique(df_plot$country))[1:min(20, length(unique(df_plot$country)))]
df_sub <- df_plot %>% filter(country %in% selected_countries)

plot_q14 <- function(prefix, label) {
  between_data <- df_plot %>%
    group_by(country) %>%
    summarise(val = first(.data[[paste0(prefix, "_between")]]), .groups = "drop")

  pA <- ggplot(between_data, aes(y = val)) + geom_boxplot() +
    theme_minimal() + labs(title = paste("Global Between:", label), x = NULL, y = NULL)

  pB <- ggplot(df_sub, aes(x = country, y = .data[[paste0(prefix, "_within")]])) +
    geom_boxplot(outlier.alpha = 0.2) +
    theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = paste("Within by Country:", label), x = NULL, y = NULL)

  pC <- ggplot(df_sub, aes(x = country, y = .data[[paste0(prefix, "_twfe")]])) +
    geom_boxplot(outlier.alpha = 0.2) +
    theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = paste("TWFE by Country:", label), x = NULL, y = NULL)

  pD <- ggplot(df_sub, aes(x = country, y = .data[[paste0(prefix, "_fd")]])) +
    geom_boxplot(outlier.alpha = 0.2) +
    theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = paste("First Differences by Country:", label), x = NULL, y = NULL)

  (pA | pB) / (pC | pD)
}

g14y <- plot_q14("y", "CO2 per GDP")
ggsave("Question14_y_boxplots.png", g14y, width = 22, height = 16, dpi = 150)

g14x <- plot_q14("x", "Equity Share")
ggsave("Question14_x_boxplots.png", g14x, width = 22, height = 16, dpi = 150)

# ===========================================================
# Question 15: Univariate statistics with standardized extremes
# ===========================================================
cat("\n====================\nQuestion 15\n====================\n")

target_vars <- c("y", "x")
transformations <- c("between", "within", "twfe", "fd")

get_stats <- function(x) {
  # if column missing or NULL
  if (is.null(x)) {
    return(tibble::tibble(
      Count = 0,
      Mean = NA_real_, Std_Dev = NA_real_,
      Std_Min_Z = NA_real_, Q1 = NA_real_, Median = NA_real_, Q3 = NA_real_,
      Std_Max_Z = NA_real_, Mean_Median_Diff = NA_real_
    ))
  }
  
  x <- x[is.finite(x)]
  n <- length(x)
  
  if (n == 0) {
    return(tibble::tibble(
      Count = 0,
      Mean = NA_real_, Std_Dev = NA_real_,
      Std_Min_Z = NA_real_, Q1 = NA_real_, Median = NA_real_, Q3 = NA_real_,
      Std_Max_Z = NA_real_, Mean_Median_Diff = NA_real_
    ))
  }
  
  mu <- mean(x)
  s  <- sd(x)
  
  zmin <- if (is.finite(s) && s > 0) (min(x) - mu) / s else NA_real_
  zmax <- if (is.finite(s) && s > 0) (max(x) - mu) / s else NA_real_
  
  tibble::tibble(
    Count = n,
    Mean = mu,
    Std_Dev = s,
    Std_Min_Z = zmin,
    Q1 = as.numeric(stats::quantile(x, 0.25, na.rm = TRUE)),
    Median = stats::median(x, na.rm = TRUE),
    Q3 = as.numeric(stats::quantile(x, 0.75, na.rm = TRUE)),
    Std_Max_Z = zmax,
    Mean_Median_Diff = mu - stats::median(x, na.rm = TRUE)
  )
}

q15_grid <- expand.grid(var = target_vars, trans = transformations, stringsAsFactors = FALSE) %>%
  dplyr::mutate(col = paste0(var, "_", trans))

q15 <- q15_grid %>%
  dplyr::mutate(stats = purrr::map(col, \(cn) {
    if (!cn %in% names(df_plot)) return(get_stats(NULL))
    get_stats(df_plot[[cn]])
  })) %>%
  tidyr::unnest(stats) %>%
  dplyr::select(
    Variable = var,
    Transformation = trans,
    dplyr::everything()
  )

print(q15 %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(.x, 4))))
write.csv(q15, "Question15_Descriptive_Stats.csv", row.names = FALSE)

# ===========================================================
# Question 16: Between vs Within correlation matrices (balanced panel)
# ===========================================================
cat("\n====================\nQuestion 16\n====================\n")

df_bal <- df_balanced %>%
  mutate(trend = year - min(year, na.rm = TRUE) + 1)

vars_list <- c(
  "cotwo_total_per_gdp", "log_gdp_per_cap", "fin_str2",
  "credit", "stock", "pop_mil", "trend",
  "log_gdp_per_cap_l1", "fin_str2_l1", "credit_l1", "stock_l1"
)
vars_list <- vars_list[vars_list %in% names(df_bal)]

# Between: replace each variable with its country mean
df_between_corr <- df_bal %>%
  group_by(country) %>%
  mutate(across(all_of(vars_list), ~mean(.x, na.rm = TRUE))) %>%
  ungroup() %>%
  select(all_of(vars_list))

# Within: deviation from the between component
df_within_corr <- df_bal %>%
  select(all_of(vars_list)) %>%
  mutate(across(everything(), as.numeric)) - df_between_corr

corr_between <- cor(df_between_corr, use = "pairwise.complete.obs")
corr_within  <- cor(df_within_corr,  use = "pairwise.complete.obs")

comparison <- data.frame(
  Variable = rownames(corr_between),
  Between_Correlation = corr_between[, "cotwo_total_per_gdp"],
  Within_Correlation  = corr_within[,  "cotwo_total_per_gdp"],
  row.names = NULL
)

cat("\nQuestion 16: Correlation with CO2/GDP\n")
print(comparison %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(.x, 4))))

cat("\nBetween correlation matrix (rounded):\n")
print(round(as.matrix(corr_between), 4))
cat("\nWithin correlation matrix (rounded):\n")
print(round(as.matrix(corr_within), 4))

find_high_corr <- function(M, threshold = 0.8) {
  idx <- which(abs(M) > threshold & upper.tri(M), arr.ind = TRUE)
  if (nrow(idx) == 0) return(data.frame())
  data.frame(
    var1 = rownames(M)[idx[,1]],
    var2 = colnames(M)[idx[,2]],
    corr = M[idx],
    row.names = NULL
  ) %>% arrange(desc(abs(corr)))
}

cat("\nHigh Correlations (> 0.8) in Within domain:\n")
print(find_high_corr(corr_within, 0.8))
cat("\nHigh Correlations (> 0.8) in Between domain:\n")
print(find_high_corr(corr_between, 0.8))

# ===========================================================
# Question 17: Autocorrelation and trend correlation (balanced)
# ===========================================================
cat("\n====================\nQuestion 17\n====================\n")
df_bal <- df_bal %>% arrange(country, year)

analysis_vars <- c("cotwo_total_per_gdp", "fin_str2", "log_gdp_per_cap", "credit")

auto_trend <- map_dfr(analysis_vars, \(v) {
  lagv <- df_bal %>% group_by(country) %>% mutate(lagv = lag(.data[[v]])) %>% ungroup() %>% pull(lagv)
  auto_corr <- cor(df_bal[[v]], lagv, use = "complete.obs")
  n_auto <- sum(is.finite(df_bal[[v]]) & is.finite(lagv))
  trend_corr <- cor(df_bal[[v]], df_bal$trend, use = "complete.obs")
  n_trend <- sum(is.finite(df_bal[[v]]))
  data.frame(
    Variable = v,
    Autocorrelation = auto_corr,
    N_Auto = n_auto,
    Trend_Correlation = trend_corr,
    N_Trend = n_trend
  )
})

print(auto_trend %>% mutate(across(where(is.numeric), ~round(.x, 4))))
write.csv(auto_trend, "Question17_Autocorr_Trend.csv", row.names = FALSE)

# ===========================================================
# Question 18: TWFE and FD correlation matrices (unbalanced)
# ===========================================================
cat("\n====================\nQuestion 18\n====================\n")

dep_var <- "cotwo_total_per_gdp"
indep_vars <- c("fin_str2", "log_gdp_per_cap", "credit", "stock", "pop_mil")
all_vars <- c(dep_var, indep_vars)

df_q18 <- df_country %>%
  filter(country %in% countries_to_keep) %>%
  arrange(country, year)

# Master key (this keeps row alignment consistent)
base <- df_q18 %>% select(country, year)

# Build TWFE residual columns (aligned by country-year) + FD columns on the full panel
for (v in all_vars) {
  
  # ---- TWFE residuals: within-demean by country, then partial out year dummies
  tmp <- df_q18 %>%
    select(country, year, value = all_of(v)) %>%
    filter(!is.na(value)) %>%
    group_by(country) %>%
    mutate(w = value - mean(value, na.rm = TRUE)) %>%
    ungroup()
  
  fit <- lm(w ~ 0 + factor(year), data = tmp)
  tmp$resid <- resid(fit)
  
  twfe_name <- paste0(v, "_twfe")
  base <- base %>%
    left_join(tmp %>% select(country, year, !!twfe_name := resid),
              by = c("country", "year"))
  
  # ---- First differences + lagged first differences (on full data, aligned)
  fd_name <- paste0(v, "_fd")
  fd_lag_name <- paste0(v, "_fd_lag")
  
  df_q18 <- df_q18 %>%
    group_by(country) %>%
    mutate(
      !!fd_name := .data[[v]] - lag(.data[[v]])
    ) %>%
    ungroup()
  
  df_q18 <- df_q18 %>%
    group_by(country) %>%
    mutate(
      !!fd_lag_name := lag(.data[[fd_name]])
    ) %>%
    ungroup()
}

# Attach FD columns to the same key frame
fd_cols <- c(paste0(all_vars, "_fd"), paste0(all_vars, "_fd_lag"))
fd_cols <- intersect(fd_cols, names(df_q18))
base <- base %>%
  left_join(df_q18 %>% select(country, year, all_of(fd_cols)),
            by = c("country", "year"))

# Correlation matrices
twfe_cols <- intersect(paste0(all_vars, "_twfe"), names(base))
fd_cols2  <- intersect(c(paste0(all_vars, "_fd"), paste0(all_vars, "_fd_lag")), names(base))

corr_twfe <- cor(base[, twfe_cols, drop = FALSE], use = "pairwise.complete.obs")
corr_fd   <- cor(base[, fd_cols2,  drop = FALSE], use = "pairwise.complete.obs")

cat("\n--- Question 18: TWFE correlation matrix ---\n")
print(round(as.matrix(corr_twfe), 4))

cat("\n--- FD correlation matrix (with lags) ---\n")
print(round(as.matrix(corr_fd), 4))

# Lists of weak correlations with dep var
dep_twfe <- paste0(dep_var, "_twfe")
dep_fd   <- paste0(dep_var, "_fd")

if (dep_twfe %in% colnames(corr_twfe)) {
  low_twfe <- names(which(abs(corr_twfe[, dep_twfe]) < 0.1))
  cat("\nPoor correlations with Dep (<0.1) in TWFE:\n"); print(low_twfe)
}

if (dep_fd %in% colnames(corr_fd)) {
  low_fd <- names(which(abs(corr_fd[, dep_fd]) < 0.1))
  cat("\nPoor correlations with Dep (<0.1) in FD:\n"); print(low_fd)
}

cat("\n--- Verification: first 30 rows of FD and lagged FD ---\n")
ver_cols <- intersect(c("country", "year", dep_fd, paste0(dep_var, "_fd_lag")), names(base))
print(head(base[, ver_cols, drop = FALSE], 30))

# ===========================================================
# Question 19: Bivariate analysis in four spaces (linear/quadratic/lowess)
# ===========================================================
cat("\n====================\nQuestion 19\n====================\n")
y_raw <- "cotwo_total_per_gdp"
x_raw <- "fin_str2"

# ensure between columns exist
if (!paste0(y_raw, "_between") %in% names(df_plot)) {
  df_plot <- df_plot %>% group_by(country) %>%
    mutate(!!paste0(y_raw, "_between") := mean(.data[[y_raw]], na.rm = TRUE)) %>% ungroup()
}
if (!paste0(x_raw, "_between") %in% names(df_plot)) {
  df_plot <- df_plot %>% group_by(country) %>%
    mutate(!!paste0(x_raw, "_between") := mean(.data[[x_raw]], na.rm = TRUE)) %>% ungroup()
}

scenarios <- list(
  list("x_within",              "y_within",              "1. One-Way Within Transformed"),
  list(paste0(x_raw, "_between"), paste0(y_raw, "_between"), "2. Between Transformed (Levels)"),
  list("x_fd",                  "y_fd",                  "3. First Differences (Growth/Changes)"),
  list("x_twfe",                "y_twfe",                "4. Two-Way Fixed Effects (TWFE)")
)

plot_bivar <- function(data, x, y, title) {
  sub <- data %>% select(all_of(c(x, y))) %>% filter(complete.cases(.))
  ggplot(sub, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(alpha = 0.25, size = 1) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
    geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE, linetype = "dotted", linewidth = 1) +
    geom_smooth(method = "loess", se = FALSE, linetype = "dashed", linewidth = 1) +
    theme_minimal() +
    labs(title = title, x = paste("Explanatory:", x), y = paste("Dependent:", y))
}

p19 <- map(scenarios, \(s) plot_bivar(df_plot, s[[1]], s[[2]], s[[3]]))
g19 <- (p19[[1]] | p19[[2]]) / (p19[[3]] | p19[[4]])
ggsave("Question19_Bivariate_Analysis.png", g19, width = 22, height = 18, dpi = 150)

# ===========================================================
# Question 20: Model comparison (Between, FE, Mundlak, TWFE, FD)
# ===========================================================
cat("\n====================\nQuestion 20\n====================\n")

if (!requireNamespace("fixest", quietly = TRUE)) install.packages("fixest")
library(fixest)
library(dplyr)

dep_var <- "cotwo_total_per_gdp"
exog_vars <- c("fin_str2", "log_gdp_per_cap", "credit", "stock", "pop_mil")

# Align sample (levels)
df_reg <- df_country %>%
  select(country, year, all_of(dep_var), all_of(exog_vars)) %>%
  filter(complete.cases(.)) %>%
  arrange(country, year)

# -------------------------
# Between (BE): country means
# -------------------------
df_be <- df_reg %>%
  group_by(country) %>%
  summarise(across(c(all_of(dep_var), all_of(exog_vars)), ~mean(.x, na.rm = TRUE)),
            .groups = "drop")

# IMPORTANT: vcov uses ONLY variables in df_be (no df_be$country)
m_be <- feols(
  as.formula(paste(dep_var, "~", paste(exog_vars, collapse = " + "))),
  data = df_be,
  vcov = ~country
)

# -------------------------
# One-way FE (country)
# -------------------------
m_fe <- feols(
  as.formula(paste(dep_var, "~", paste(exog_vars, collapse = " + "), "| country")),
  data = df_reg,
  vcov = ~country
)

# -------------------------
# Mundlak: levels + country means of X (drop identical mean vars)
# -------------------------
df_m <- df_reg %>%
  group_by(country) %>%
  mutate(across(all_of(exog_vars), ~mean(.x, na.rm = TRUE), .names = "{.col}_mean")) %>%
  ungroup()

mean_vars <- paste0(exog_vars, "_mean")

# Drop mean vars that are numerically identical to the level var (time-invariant within country)
keep_mean <- mean_vars[!vapply(seq_along(exog_vars), function(i) {
  v  <- exog_vars[i]
  vm <- mean_vars[i]
  all(is.finite(df_m[[v]] - df_m[[vm]]) & abs(df_m[[v]] - df_m[[vm]]) < 1e-12)
}, logical(1))]

# IMPORTANT: vcov uses ONLY variables in df_m (no df_m$country)
m_mundlak <- feols(
  as.formula(paste(dep_var, "~", paste(c(exog_vars, keep_mean), collapse = " + "))),
  data = df_m,
  vcov = ~country
)

# -------------------------
# TWFE (country + year)
# -------------------------
m_twfe <- feols(
  as.formula(paste(dep_var, "~", paste(exog_vars, collapse = " + "), "| country + year")),
  data = df_reg,
  vcov = ~country
)

# -------------------------
# First differences (FD)
#   Skip if dy is constant / empty
# -------------------------
df_fd <- df_reg %>%
  group_by(country) %>%
  arrange(year) %>%
  mutate(dy = .data[[dep_var]] - lag(.data[[dep_var]])) %>%
  ungroup()

for (v in exog_vars) {
  df_fd <- df_fd %>%
    group_by(country) %>%
    arrange(year) %>%
    mutate(!!paste0("d_", v) := .data[[v]] - lag(.data[[v]])) %>%
    ungroup()
}

dvars <- paste0("d_", exog_vars)

df_fd_use <- df_fd %>%
  select(country, year, dy, all_of(dvars)) %>%
  filter(is.finite(dy)) %>%
  filter(complete.cases(.))

m_fd <- NULL
if (nrow(df_fd_use) == 0) {
  message("FD: no usable observations after differencing; skipping FD.")
} else if (dplyr::n_distinct(df_fd_use$dy) <= 1) {
  message("FD: dy is constant; skipping FD.")
} else {
  m_fd <- feols(
    as.formula(paste("dy ~", paste(dvars, collapse = " + "), "| year")),
    data = df_fd_use,
    vcov = ~country
  )
}

# -------------------------
# Output table
# IMPORTANT: do NOT add se='cluster' here, because vcov is already set in each model
# -------------------------
models <- list(BE = m_be, FE = m_fe, MUNDLAK = m_mundlak, TWFE = m_twfe)
if (!is.null(m_fd)) models$FD <- m_fd

tab <- etable(models)
print(tab)
writeLines(capture.output(tab), "Question20_Regression_Table.txt")

# ===========================================================
# Question 21: Anderson–Hsiao ARDL(1,1) in first differences + IV
# ===========================================================
cat("\n====================\nQuestion 21\n====================\n")
y_name <- "cotwo_total_per_gdp"
x_name <- "fin_str2"

df_ah <- df_plot %>% arrange(country, year)

# Required pieces (already in df_plot): y_fd, x_fd
df_ah <- df_ah %>%
  group_by(country) %>%
  mutate(
    dy_lag1 = lag(y_fd, 1),
    dx_lag1 = lag(x_fd, 1),
    y_level_lag2 = lag(.data[[y_name]], 2),
    x_level_lag2 = lag(.data[[x_name]], 2)
  ) %>%
  ungroup()

# 21.1 Univariate stats
dynamic_vars <- c("y_fd", "dy_lag1", "x_fd", "dx_lag1", "y_level_lag2", "x_level_lag2")
q21_1 <- map_dfr(dynamic_vars, \(v) {
  x <- df_ah[[v]]
  x <- x[is.finite(x)]
  data.frame(
    Variable = v,
    Mean = mean(x), Median = median(x),
    Std_Error = sd(x) / sqrt(length(x)),
    Std_Dev = sd(x),
    Skewness = moments::skewness(x),
    Kurtosis = moments::kurtosis(x),
    Min = min(x), Max = max(x), N = length(x)
  )
})
write.csv(q21_1, "Question21_1_ARDL_Stats.csv", row.names = FALSE)
cat("\nQuestion 21.1 (head):\n")
print(q21_1 %>% mutate(across(where(is.numeric), ~round(.x, 4))) %>% head(10))

# 21.2 Correlation matrix for instrument strength check
matrix_vars <- c("dy_lag1", "x_fd", "dx_lag1", "y_level_lag2", "x_level_lag2")
df_corr_ah <- df_ah %>% select(all_of(matrix_vars)) %>% filter(complete.cases(.))
corr_matrix_ah <- cor(df_corr_ah)
cat("\nQuestion 21.2: ARDL Correlation Matrix\n")
print(round(corr_matrix_ah, 4))
cat(sprintf("\nCorr(dy_lag1, y_level_lag2): %.4f\n", corr_matrix_ah["dy_lag1", "y_level_lag2"]))
cat(sprintf("Corr(dx_lag1, x_level_lag2): %.4f\n", corr_matrix_ah["dx_lag1", "x_level_lag2"]))

# 21.3 Fisher-type ADF on first differences
cat("\nQuestion 21.3: Fisher-type ADF (by country)\n")
adf_y <- fisher_adf(df_ah %>% select(country, y_fd), id = "country", value = "y_fd")
adf_x <- fisher_adf(df_ah %>% select(country, x_fd), id = "country", value = "x_fd")
cat(sprintf("ΔY: N=%d, Fisher=%.4f, df=%d, p=%.4f\n", adf_y$n, adf_y$fisher, adf_y$df, adf_y$p))
cat(sprintf("ΔX: N=%d, Fisher=%.4f, df=%d, p=%.4f\n", adf_x$n, adf_x$fisher, adf_x$df, adf_x$p))

# 21.4 OLS ARDL in differences (year FE via factor(year)), clustered SE by country
ctrl_raw_vars <- c("log_gdp_per_cap", "credit", "stock", "pop_mil")
for (v in ctrl_raw_vars) {
  df_ah <- df_ah %>% group_by(country) %>%
    mutate(!!paste0(v, "_fd") := .data[[v]] - lag(.data[[v]])) %>%
    ungroup()
}
ctrls_fd <- paste0(ctrl_raw_vars, "_fd")

df_ardl <- df_ah %>%
  select(country, year, y_fd, dy_lag1, x_fd, dx_lag1, all_of(ctrls_fd)) %>%
  filter(complete.cases(.))

ols_ardl <- lm(
  y_fd ~ dy_lag1 + x_fd + dx_lag1 + .,
  data = df_ardl %>% mutate(year = factor(year)) %>%
    select(-country) %>% # keep year + regressors
    mutate(country = df_ardl$country) # put country back at end
)

# Cluster-robust vcov by country
vc_ols <- sandwich::vcovCL(ols_ardl, cluster = df_ardl$country, type = "HC1")
res_ols_tab <- lmtest::coeftest(ols_ardl, vcov. = vc_ols)
cat("\nQuestion 21.4: OLS ARDL (clustered by country)\n")
print(res_ols_tab)
writeLines(paste(capture.output(res_ols_tab), collapse = "\n"), "Question21.4_Regression_Table.txt")

# 21.5 Anderson–Hsiao IV via ivreg
df_iv <- df_ah %>%
  select(country, year, y_fd, dy_lag1, x_fd, dx_lag1, y_level_lag2, x_level_lag2, all_of(ctrls_fd)) %>%
  filter(complete.cases(.)) %>%
  mutate(year = factor(year))

iv_form <- as.formula(
  paste0(
    "y_fd ~ dy_lag1 + x_fd + dx_lag1 + ", paste(ctrls_fd, collapse = " + "), " + year",
    " | dx_lag1 + ", paste(ctrls_fd, collapse = " + "), " + year + y_level_lag2 + x_level_lag2"
  )
)

iv_ah <- AER::ivreg(iv_form, data = df_iv)

vc_iv <- sandwich::vcovCL(iv_ah, cluster = df_iv$country, type = "HC1")
iv_tab <- lmtest::coeftest(iv_ah, vcov. = vc_iv)

cat("\nQuestion 21.5: Anderson–Hsiao IV (clustered by country)\n")
print(iv_tab)
writeLines(paste(capture.output(iv_tab), collapse = "\n"), "Question21.5_Regression_Table.txt")

# 21.6 First-stage diagnostics (roughly mirroring Python):
# Run separate first-stage regressions with excluded instruments and compute joint F test.
endogenous_vars <- c("dy_lag1", "x_fd")
instruments <- c("y_level_lag2", "x_level_lag2")
exog_fs <- c("dx_lag1", ctrls_fd, "year")

first_stage_results <- map_dfr(endogenous_vars, \(endo) {
  fs_form <- as.formula(
    paste0(endo, " ~ ", paste(c(exog_fs, instruments), collapse = " + "))
  )
  fs <- lm(fs_form, data = df_iv)
  # Joint F-test: instruments excluded
  # Use linearHypothesis with robust vcov if installed; otherwise classic F
  # We'll do classic F for portability:
  # Compare with reduced model excluding instruments:
  fs_red_form <- as.formula(
    paste0(endo, " ~ ", paste(exog_fs, collapse = " + "))
  )
  fs_red <- lm(fs_red_form, data = df_iv)
  an <- anova(fs_red, fs)
  fstat <- an$F[2]
  pval  <- an$`Pr(>F)`[2]
  data.frame(
    Endogenous_Var = endo,
    R2 = summary(fs)$r.squared,
    Adj_R2 = summary(fs)$adj.r.squared,
    F_Instruments = fstat,
    p_value = pval
  )
})
cat("\nQuestion 21.6: First-stage diagnostics (classic F-test)\n")
print(first_stage_results %>% mutate(across(where(is.numeric), ~round(.x, 4))))
write.csv(first_stage_results, "Question21.6_FirstStage.csv", row.names = FALSE)

# 21.7 IRF (t=1..4) from IV coefficients
b <- coef(iv_ah)
beta_y <- unname(b["dy_lag1"])
beta_1 <- unname(b["x_fd"])
beta_2 <- unname(b["dx_lag1"])

irf_values <- c(
  beta_1,
  beta_y * beta_1 + beta_2,
  (beta_y^2) * beta_1 + beta_y * beta_2,
  (beta_y^3) * beta_1 + (beta_y^2) * beta_2
)
periods <- 1:4
irf_df <- data.frame(t = periods, irf = irf_values)

cat("\nQuestion 21.7: IRF values\n")
print(irf_df)

p_irf <- ggplot(irf_df, aes(x = t, y = irf)) +
  geom_line() + geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(title = "IRF: response of Δ CO2/GDP to 1-unit shock in Δ Equity Share",
       x = "Periods (years after shock)", y = "Impact on emissions growth")
ggsave("Question21_7_IRF_Plot.png", p_irf, width = 10, height = 6, dpi = 150)

# 21.8 Long-run propensity
beta_lt <- (beta_1 + beta_2) / (1 - beta_y)
cat("\n==================================================\n")
cat("Question 21.8: Long-Run Propensity (LRP)\n")
cat("==================================================\n")
cat(sprintf("beta_1 (short-run contemporaneous): %.6f\n", beta_1))
cat(sprintf("beta_2 (short-run lagged):          %.6f\n", beta_2))
cat(sprintf("beta_y (persistence):              %.6f\n", beta_y))
cat("--------------------------------------------------\n")
cat(sprintf("beta_LT (long-run coefficient):    %.6f\n", beta_lt))
cat("==================================================\n")
sum_short_run <- beta_1 + beta_2
multiplier <- 1 / (1 - beta_y)
cat(sprintf("\nDiagnostic breakdown:\nSum short-run: %.6f\nDynamic multiplier: %.4f\n",
            sum_short_run, multiplier))

cat("\n\n--- DONE ---\n")
