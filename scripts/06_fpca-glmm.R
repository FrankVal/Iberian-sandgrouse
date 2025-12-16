## ============================
## Code 06 — FPCA-GLMM framework
# FPCA-GLMM (functional logistic regression; scalar-on-function)
# Persistence (1) vs Extinction (0) using 2005–2019 full time series
# of landscape metrics per UTMCODE.
#
# Notes:
# - Only Status ∈ {Persistence, Extinction} are retained (Colonization removed)
# - Species is ignored (optional filter if you want)
# - Random effect: (1 | RegSect). No Year random effect (Year enters via FPCA)
# - Metrics: AREA_MN, ED, LPI, NP, CA  (CA)
# - Variance target: 0.80; max 2 PCs per metric
# ============================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(glmmTMB)
  library(pROC)
  library(performance)
})

# ----------------------------
# 0) read & prep
# ----------------------------
data_path <- "E:/path/to/your/status_metrics_file.csv"  # <- set your path (csv/tsv)
variance_target <- 0.80
max_pc <- 2
years_vec <- 2005:2019

dat_raw <- read.csv(data_path, sep = ";", stringsAsFactors = FALSE, check.names = FALSE)

# required columns (species is not required; UTMCODE is only an id, not a predictor)
req_cols <- c("UTMCODE","Status","Year","RegSect","AREA_MN","ED","LPI","NP","CA")
stopifnot(all(req_cols %in% names(dat_raw)))

dat_long <- dat_raw %>%
  mutate(
    Year    = as.integer(Year),
    Status  = as.character(Status),
    RegSect = as.character(RegSect)
  ) %>%
  filter(Year %in% years_vec) %>%
  select(all_of(req_cols))

# final (binary) status per UTMCODE, keeping ONLY Persistence / Extinction
status_final <- dat_long %>%
  group_by(UTMCODE) %>%
  slice_max(Year, with_ties = FALSE) %>%
  ungroup() %>%
  filter(Status %in% c("Persistence","Extinction")) %>%
  transmute(UTMCODE, status01 = as.integer(Status == "Persistence"))

# keep only modeled UTMCODEs (colonization removed here)
dat_long <- dat_long %>%
  semi_join(status_final, by = "UTMCODE")

metrics <- c("AREA_MN","ED","LPI","NP","CA")

# ----------------------------
# helpers
# ----------------------------
interp_row <- function(x){
  tt <- seq_along(x)
  ok <- is.finite(x)
  if (!any(ok)) return(rep(0, length(x)))
  approx(tt[ok], x[ok], xout = tt, rule = 2)$y
}

make_wide <- function(long_df, metric, years = years_vec){
  vals <- long_df %>%
    mutate(value = .data[[metric]]) %>%
    select(UTMCODE, Year, value)
  
  grid <- tidyr::expand_grid(
    UTMCODE = sort(unique(long_df$UTMCODE)),
    Year    = years
  )
  
  grid %>%
    left_join(vals, by = c("UTMCODE","Year")) %>%
    pivot_wider(names_from = Year, values_from = value) %>%
    arrange(UTMCODE)
}

scale_to_int <- function(w){
  w <- as.numeric(w)
  w[!is.finite(w) | w < 0] <- 0
  mp <- suppressWarnings(min(w[w > 0]))
  if (!is.finite(mp)) return(integer(length(w)))
  as.integer(round(w / mp))  # smallest positive weight -> 1
}

best_threshold <- function(y, p){
  roc_obj <- pROC::roc(response = y, predictor = p, quiet = TRUE, direction = "<")
  co <- pROC::coords(roc_obj, "best", ret = "threshold", best.method = "youden", transpose = FALSE)
  thr <- suppressWarnings(as.numeric(if (is.list(co)) co$threshold else co))
  if (length(thr) == 0 || all(is.na(thr))) return(0.5)
  if (length(thr) > 1) thr <- mean(thr, na.rm = TRUE)
  thr
}

safe_r2 <- function(fit, has_re){
  has_pos_var <- FALSE
  vc <- try(VarCorr(fit)$cond, silent = TRUE)
  if (!inherits(vc, "try-error") && !is.null(vc)) {
    vv <- unlist(lapply(vc, function(x) if (is.numeric(x)) x[1] else NA_real_))
    has_pos_var <- any(is.finite(vv) & vv > 0)
  }
  
  if (!has_re || !has_pos_var) {
    nk <- try(suppressWarnings(suppressMessages(performance::r2_nagelkerke(fit))), silent = TRUE)
    if (!inherits(nk, "try-error") && !is.null(nk$R2_Nagelkerke)) {
      r2m <- as.numeric(nk$R2_Nagelkerke)
      return(c(r2_m = r2m, r2_c = r2m))
    }
    return(c(r2_m = NA_real_, r2_c = NA_real_))
  }
  
  out <- try(suppressWarnings(suppressMessages(performance::r2_nakagawa(fit, tolerance = 0))), silent = TRUE)
  if (!inherits(out, "try-error") && !is.null(out$r2_marginal)) {
    r2m <- as.numeric(out$r2_marginal)
    r2c <- if (!is.null(out$r2_conditional)) as.numeric(out$r2_conditional) else r2m
    return(c(r2_m = r2m, r2_c = r2c))
  }
  
  nk <- try(suppressWarnings(suppressMessages(performance::r2_nagelkerke(fit))), silent = TRUE)
  if (!inherits(nk, "try-error") && !is.null(nk$R2_Nagelkerke)) {
    r2m <- as.numeric(nk$R2_Nagelkerke)
    return(c(r2_m = r2m, r2_c = r2m))
  }
  c(r2_m = NA_real_, r2_c = NA_real_)
}

# ----------------------------
# FPCA per metric (train -> project to test)
# ----------------------------
train_project_metric_pca <- function(train_long, test_long, metric,
                                     years = years_vec,
                                     variance_cut = variance_target,
                                     max_pc = max_pc){
  
  wide_tr <- make_wide(train_long, metric, years)
  m_tr <- as.matrix(wide_tr[, -1, drop = FALSE])
  if (anyNA(m_tr)) m_tr <- t(apply(m_tr, 1, interp_row))
  
  pr <- prcomp(m_tr, center = TRUE, scale. = FALSE)
  varexp <- cumsum(pr$sdev^2) / sum(pr$sdev^2)
  
  k_var <- which(varexp >= variance_cut)[1]
  if (is.na(k_var)) k_var <- ncol(m_tr)
  k <- max(1, min(k_var, max_pc))  # <= 2 PCs
  
  s_tr <- pr$x[, 1:k, drop = FALSE]
  mu_s <- apply(s_tr, 2, mean)
  sd_s <- apply(s_tr, 2, sd); sd_s[sd_s == 0] <- 1
  s_tr_sc <- scale(s_tr, center = mu_s, scale = sd_s)
  colnames(s_tr_sc) <- paste0(metric, "_pc", seq_len(k))
  
  wide_te <- make_wide(test_long, metric, years)
  m_te <- as.matrix(wide_te[, -1, drop = FALSE])
  if (anyNA(m_te)) m_te <- t(apply(m_te, 1, interp_row))
  
  m_te_ctr <- sweep(m_te, 2, pr$center, "-")
  s_te <- m_te_ctr %*% pr$rotation[, 1:k, drop = FALSE]
  s_te_sc <- sweep(s_te, 2, mu_s, "-")
  s_te_sc <- sweep(s_te_sc, 2, sd_s, "/")
  colnames(s_te_sc) <- paste0(metric, "_pc", seq_len(k))
  
  list(
    train = tibble(UTMCODE = wide_tr$UTMCODE) %>% bind_cols(as.data.frame(s_tr_sc)),
    test  = tibble(UTMCODE = wide_te$UTMCODE) %>% bind_cols(as.data.frame(s_te_sc))
  )
}

# ----------------------------
# build train/test frames
# ----------------------------
build_frames <- function(train_long, test_long){
  pcs_tr_list <- list(UTMCODE = tibble(UTMCODE = sort(unique(train_long$UTMCODE))))
  pcs_te_list <- list(UTMCODE = tibble(UTMCODE = sort(unique(test_long$UTMCODE))))
  
  for (m in metrics){
    pp <- train_project_metric_pca(train_long, test_long, metric = m)
    pcs_tr_list[[m]] <- pp$train
    pcs_te_list[[m]] <- pp$test
  }
  
  pcs_tr <- reduce(pcs_tr_list, full_join, by = "UTMCODE")
  pcs_te <- reduce(pcs_te_list, full_join, by = "UTMCODE")
  
  # RegSect per UTMCODE (majority vote if inconsistent)
  meta_df <- dat_long %>%
    group_by(UTMCODE) %>%
    summarise(
      RegSect = names(sort(table(RegSect), decreasing = TRUE))[1],
      .groups = "drop"
    ) %>%
    mutate(RegSect = factor(RegSect))
  
  train_df <- pcs_tr %>%
    left_join(status_final, by = "UTMCODE") %>%
    left_join(meta_df, by = "UTMCODE")
  
  test_df <- pcs_te %>%
    left_join(status_final, by = "UTMCODE") %>%
    left_join(meta_df, by = "UTMCODE")
  
  # class weights (train only)
  prev <- train_df %>% count(status01) %>% pivot_wider(names_from = status01, values_from = n, values_fill = 0)
  n0 <- as.numeric(prev$`0`); n1 <- as.numeric(prev$`1`)
  n_tot <- n0 + n1
  
  if (isTRUE(n_tot > 0 && n0 > 0 && n1 > 0)) {
    w0 <- n_tot/(2*n0); w1 <- n_tot/(2*n1)
    train_df <- train_df %>% mutate(w = ifelse(status01 == 1L, w1, w0))
  } else {
    train_df$w <- 1
  }
  test_df$w <- 1
  
  train_df$w <- scale_to_int(train_df$w)
  
  list(train = train_df, test = test_df)
}

# ----------------------------
# fit + evaluate one split
# ----------------------------
fit_eval <- function(train_df, test_df){
  pc_cols <- grep("_pc\\d+$", names(train_df), value = TRUE)
  if (length(pc_cols)) {
    sds <- sapply(train_df[pc_cols], sd, na.rm = TRUE)
    pc_cols <- pc_cols[sds > 0]
  }
  if (!length(pc_cols)) stop("No pc predictors to fit.")
  
  train_df <- train_df %>% mutate(RegSect = factor(RegSect))
  test_df  <- test_df  %>% mutate(RegSect = factor(RegSect, levels = levels(train_df$RegSect)))
  
  re_terms <- character(0)
  if (nlevels(droplevels(train_df$RegSect)) >= 2) re_terms <- c(re_terms, "(1|RegSect)")
  
  rhs <- paste(c(pc_cols, re_terms), collapse = " + ")
  fm <- as.formula(paste("status01 ~", rhs))
  
  null_rhs <- if (length(re_terms)) paste("1 +", paste(re_terms, collapse = " + ")) else "1"
  null_fm <- as.formula(paste("status01 ~", null_rhs))
  
  fit_tr  <- suppressWarnings(glmmTMB(fm, data = train_df, family = binomial(), weights = w))
  null_tr <- suppressWarnings(glmmTMB(null_fm, data = train_df, family = binomial(), weights = w))
  
  pr_te <- predict(fit_tr, newdata = test_df, type = "response", allow.new.levels = TRUE)
  y_te  <- test_df$status01
  
  roc_te   <- pROC::roc(response = y_te, predictor = pr_te, quiet = TRUE, direction = "<")
  auc_te   <- as.numeric(pROC::auc(roc_te))
  tjur_te  <- mean(pr_te[y_te == 1], na.rm = TRUE) - mean(pr_te[y_te == 0], na.rm = TRUE)
  brier_te <- mean((pr_te - y_te)^2, na.rm = TRUE)
  
  pr_null_te <- predict(null_tr, newdata = test_df, type = "response", allow.new.levels = TRUE)
  ll_full_te <- sum(y_te*log(pmax(pr_te,1e-12)) + (1 - y_te)*log(pmax(1 - pr_te,1e-12)))
  ll_null_te <- sum(y_te*log(pmax(pr_null_te,1e-12)) + (1 - y_te)*log(pmax(1 - pr_null_te,1e-12)))
  mcfad_oos  <- 1 - (ll_full_te / ll_null_te)
  
  has_re <- length(re_terms) > 0
  r2_vals <- safe_r2(fit_tr, has_re)
  r2_m <- unname(r2_vals["r2_m"])
  r2_c <- unname(r2_vals["r2_c"])
  
  thr <- best_threshold(y_te, pr_te)
  pred_cls <- as.integer(pr_te >= thr)
  tab <- table(observed = y_te, predicted = pred_cls)
  acc <- sum(diag(tab)) / sum(tab)
  
  list(
    metrics = tibble(
      auc = auc_te,
      tjur = tjur_te,
      brier = brier_te,
      mcfadden_oos = mcfad_oos,
      accuracy = acc,
      threshold = thr,
      r2_marginal = r2_m,
      r2_conditional = r2_c
    ),
    table = tab,
    fit = fit_tr
  )
}

# ----------------------------
# k-fold CV (stratified)
# ----------------------------
cv_run <- function(k = 5, seed = 12){
  set.seed(seed)
  
  idx <- status_final %>%
    group_by(status01) %>%
    mutate(fold = sample(rep(1:k, length.out = n()))) %>%
    ungroup()
  
  out_list <- vector("list", k)
  
  for (i in 1:k){
    test_ids  <- idx$UTMCODE[idx$fold == i]
    train_ids <- setdiff(idx$UTMCODE, test_ids)
    
    train_long <- dat_long %>% filter(UTMCODE %in% train_ids)
    test_long  <- dat_long %>% filter(UTMCODE %in% test_ids)
    
    frames <- build_frames(train_long, test_long)
    fe <- fit_eval(frames$train, frames$test)
    
    out_list[[i]] <- fe$metrics %>% mutate(fold = i)
  }
  
  fold_metrics <- bind_rows(out_list)
  
  summary <- fold_metrics %>%
    summarise(
      auc_mean = mean(auc), auc_sd = sd(auc),
      tjur_mean = mean(tjur), tjur_sd = sd(tjur),
      brier_mean = mean(brier), brier_sd = sd(brier),
      acc_mean = mean(accuracy), acc_sd = sd(accuracy)
    )
  
  list(fold_metrics = fold_metrics, summary = summary)
}

# ---- run CV ----
cv_res <- cv_run(k = 5, seed = 12)
print(cv_res$summary)

# ============================
# full-data FPCA-GLMM + beta(t)
# ============================

fpca_full <- function(long_df, metric, years = years_vec,
                      variance_cut = variance_target, max_pc = max_pc){
  wide <- long_df %>%
    mutate(value = .data[[metric]]) %>%
    select(UTMCODE, Year, value) %>%
    complete(UTMCODE, Year = years, fill = list(value = NA_real_)) %>%
    pivot_wider(names_from = Year, values_from = value) %>%
    arrange(UTMCODE)
  
  m <- as.matrix(wide[, -1, drop = FALSE])
  if (anyNA(m)) m <- t(apply(m, 1, interp_row))
  
  pr <- prcomp(m, center = TRUE, scale. = FALSE)
  varexp <- cumsum(pr$sdev^2) / sum(pr$sdev^2)
  
  k_var <- which(varexp >= variance_cut)[1]
  if (is.na(k_var)) k_var <- ncol(m)
  k <- max(1, min(k_var, max_pc))
  
  s <- pr$x[, 1:k, drop = FALSE]
  mu_s <- colMeans(s)
  sd_s <- apply(s, 2, sd); sd_s[sd_s == 0] <- 1
  s_sc <- scale(s, center = mu_s, scale = sd_s)
  colnames(s_sc) <- paste0(metric, "_pc", seq_len(k))
  
  list(metric = metric, utm = wide$UTMCODE, scores = as.data.frame(s_sc),
       pr = pr, k = k, mu_s = mu_s, sd_s = sd_s, years = years)
}

build_full_frame <- function(dat_long){
  pcs <- lapply(metrics, fpca_full, long_df = dat_long)
  pc_df <- reduce(
    lapply(pcs, function(x) tibble(UTMCODE = x$utm) %>% bind_cols(x$scores)),
    full_join, by = "UTMCODE"
  )
  
  meta_df <- dat_long %>%
    group_by(UTMCODE) %>%
    summarise(RegSect = names(sort(table(RegSect), decreasing = TRUE))[1], .groups = "drop") %>%
    mutate(RegSect = factor(RegSect))
  
  full_df <- pc_df %>%
    left_join(status_final, by = "UTMCODE") %>%
    left_join(meta_df, by = "UTMCODE")
  
  prev <- full_df %>% count(status01) %>% pivot_wider(names_from = status01, values_from = n, values_fill = 0)
  n0 <- as.numeric(prev$`0`); n1 <- as.numeric(prev$`1`); n_tot <- n0 + n1
  w <- if (isTRUE(n_tot > 0 && n0 > 0 && n1 > 0)) ifelse(full_df$status01 == 1L, n_tot/(2*n1), n_tot/(2*n0)) else 1
  full_df$w <- w
  
  list(data = full_df, pcs = pcs)
}

fit_full <- function(full_df){
  pc_cols <- grep("_pc\\d+$", names(full_df), value = TRUE)
  sds <- sapply(full_df[pc_cols], sd, na.rm = TRUE)
  pc_cols <- pc_cols[sds > 0]
  stopifnot(length(pc_cols) > 0)
  
  full_df <- full_df %>% mutate(RegSect = factor(RegSect))
  form <- as.formula(paste0("status01 ~ ", paste(pc_cols, collapse = " + "), " + (1|RegSect)"))
  glmmTMB(form, data = full_df, family = binomial(), weights = w)
}

metric_year_sd <- function(dat_long, metric, years = years_vec){
  wide <- dat_long %>%
    mutate(value = .data[[metric]]) %>%
    select(UTMCODE, Year, value) %>%
    complete(UTMCODE, Year = years, fill = list(value = NA_real_)) %>%
    pivot_wider(names_from = Year, values_from = value) %>%
    arrange(UTMCODE)
  m <- as.matrix(wide[, -1, drop = FALSE])
  apply(m, 2, sd, na.rm = TRUE)
}

reconstruct_betas <- function(fit, pcs_list){
  b <- as.numeric(fixef(fit)$cond)
  bn <- names(fixef(fit)$cond)
  
  map_dfr(pcs_list, function(pp){
    metric <- pp$metric
    k <- pp$k
    years <- pp$years
    
    coefs <- sapply(seq_len(k), function(j){
      nm <- paste0(metric, "_pc", j)
      idx <- match(nm, bn)
      if (is.na(idx)) 0 else b[idx]
    })
    
    coefs_rescaled <- coefs / pp$sd_s[seq_len(k)]
    phi <- pp$pr$rotation[, seq_len(k), drop = FALSE]
    beta_t <- as.numeric(phi %*% coefs_rescaled)
    
    tibble(metric = metric, Year = years, beta_raw = beta_t)
  })
}

built <- build_full_frame(dat_long)
full_df <- built$data
pcs <- built$pcs

fit <- fit_full(full_df)
betas_raw <- reconstruct_betas(fit, pcs)

sds_df <- bind_rows(lapply(metrics, function(m){
  tibble(metric = m, Year = years_vec, sd_year = metric_year_sd(dat_long, m, years_vec))
}))

betas <- betas_raw %>%
  left_join(sds_df, by = c("metric","Year")) %>%
  mutate(beta_std = beta_raw * sd_year)  # +1 SD change at that year

importance_tbl <- betas %>%
  group_by(metric) %>%
  summarise(
    importance_sd = sqrt(mean(beta_std^2, na.rm = TRUE)),
    mean_all_sd   = mean(beta_std, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(importance_sd))

print(importance_tbl)
