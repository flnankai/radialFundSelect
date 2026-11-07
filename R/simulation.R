# The simulation retains only additive sufficient statistics in memory.
# No panel, repetition, block, checkpoint, or RNG state is written to disk.

scenario_inventory <- function() {
  rows <- list()
  for (distribution in c("normal", "t3", "mixture")) {
    for (T in c(60L, 120L)) {
      deltas <- if (T == 60L) seq_len(9L) * 0.06 else seq_len(9L) * 0.04
      for (delta in deltas) {
        rows[[length(rows) + 1L]] <- data.frame(
          case_id = .rfs_case_id(distribution, T, 0.5, delta),
          distribution = distribution, N = 200L, T = T, rho = 0.5,
          delta = delta, signal_grid = TRUE,
          rho_sensitivity = T == 60L && abs(delta - 0.48) < 1e-12,
          stringsAsFactors = FALSE
        )
      }
    }
    for (rho in (seq_len(9L) / 10)[-5L]) {
      rows[[length(rows) + 1L]] <- data.frame(
        case_id = .rfs_case_id(distribution, 60L, rho, 0.48),
        distribution = distribution, N = 200L, T = 60L, rho = rho,
        delta = 0.48, signal_grid = FALSE, rho_sensitivity = TRUE,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

.rfs_case_id <- function(distribution, T, rho, delta) {
  sprintf("%s_T%d_rho%.2f_delta%.2f", distribution, T, rho, delta)
}

.rfs_integer <- function(x, name, lower = 1, upper = .Machine$integer.max) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) ||
      x != floor(x) || x < lower || x > upper) {
    stop(name, " must be an integer between ", lower, " and ", upper, ".",
         call. = FALSE)
  }
  as.integer(x)
}

.rfs_rng_snapshot <- function() {
  list(kind = RNGkind(),
       seed = if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
         get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL)
}

.rfs_rng_restore <- function(state) {
  do.call(RNGkind, as.list(state$kind))
  if (is.null(state$seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  } else {
    assign(".Random.seed", state$seed, envir = .GlobalEnv)
  }
  invisible(NULL)
}

.rfs_seed_table <- function(seed, reps) {
  old <- .rfs_rng_snapshot()
  on.exit(.rfs_rng_restore(old))
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  set.seed(seed)
  stream <- get(".Random.seed", envir = .GlobalEnv)
  seeds <- vector("list", 8L)
  # Eight independent streams: (T=60,120) x (nuisance, normal, t3, mixture).
  # Repetition r uses substream r, independently of scheduling or block size.
  for (j in seq_len(8L)) {
    substream <- stream
    seeds[[j]] <- vector("list", reps)
    for (r in seq_len(reps)) {
      seeds[[j]][[r]] <- substream
      substream <- parallel::nextRNGSubStream(substream)
    }
    stream <- parallel::nextRNGStream(stream)
  }
  seeds
}

.rfs_draw_nuisance <- function(T, N, rng_seed) {
  assign(".Random.seed", rng_seed, envir = .GlobalEnv)
  mu_f <- stats::runif(3L, 0, 1)
  cov_f <- 0.5 ^ abs(outer(seq_len(3L), seq_len(3L), "-"))
  factors <- sweep(matrix(stats::rnorm(T * 3L), T, 3L) %*% chol(cov_f),
                   2L, mu_f, "+")
  beta <- cbind(stats::runif(N, 0.2, 2), stats::runif(N, -1, 1.5),
                stats::runif(N, -1.5, 1.5))
  order <- sample.int(N)
  n_positive <- as.integer(round(0.1 * N))
  positive <- order[seq_len(n_positive)]
  negative <- order[n_positive + seq_len(n_positive)]
  alpha <- numeric(N)
  alpha[positive] <- stats::runif(n_positive, 0.15, 0.30)
  alpha[negative] <- stats::runif(n_positive, -0.30, -0.20)
  list(factors = factors, beta = beta, mu_f = mu_f, alpha_base = alpha,
       positive = positive, negative = negative)
}

.rfs_draw_error_sources <- function(T, N, distribution, rng_seed) {
  assign(".Random.seed", rng_seed, envir = .GlobalEnv)
  components <- matrix(stats::rnorm(T * N), T, N)
  scale <- switch(distribution,
                  normal = rep(1, T),
                  t3 = 1 / sqrt(stats::rchisq(T, df = 3)),
                  mixture = ifelse(stats::runif(T) < 0.8, 1, 3) / sqrt(2.6))
  list(components = components, scale = scale)
}

.rfs_correlate_errors <- function(sources, rho) {
  gaussian <- sources$components
  innovation_sd <- sqrt(1 - rho^2)
  for (i in seq.int(2L, ncol(gaussian))) {
    gaussian[, i] <- rho * gaussian[, i - 1L] +
      innovation_sd * sources$components[, i]
  }
  gaussian * sources$scale
}

simulate_panel <- function(T = 60L, N = 200L, distribution = "normal",
                           rho = 0.5, delta = 0.06, seed = 20261105L,
                           replicate = 1L) {
  T <- .rfs_integer(T, "T", 5L)
  N <- .rfs_integer(N, "N", 10L)
  if (!T %in% c(60L, 120L))
    stop("The study uses T = 60 or 120.", call. = FALSE)
  distribution <- match.arg(distribution, c("normal", "t3", "mixture"))
  if (length(rho) != 1L || !is.finite(rho) || rho < 0 || rho >= 1)
    stop("rho must be in [0, 1).", call. = FALSE)
  if (length(delta) != 1L || !is.finite(delta) || delta < 0)
    stop("delta must be finite and nonnegative.", call. = FALSE)
  seed <- .rfs_integer(seed, "seed", 0L)
  replicate <- .rfs_integer(replicate, "replicate")
  old <- .rfs_rng_snapshot()
  on.exit(.rfs_rng_restore(old))
  seeds <- .rfs_seed_table(seed, replicate)
  offset <- if (T == 60L) 0L else 4L
  nuisance <- .rfs_draw_nuisance(T, N, seeds[[offset + 1L]][[replicate]])
  code <- match(distribution, c("normal", "t3", "mixture"))
  sources <- .rfs_draw_error_sources(T, N, distribution,
                                   seeds[[offset + code + 1L]][[replicate]])
  alpha <- nuisance$alpha_base
  alpha[nuisance$positive] <- alpha[nuisance$positive] + delta
  errors <- .rfs_correlate_errors(sources, rho)
  Y <- sweep(nuisance$factors %*% t(nuisance$beta) + errors, 2L, alpha, "+")
  list(Y = Y, factors = nuisance$factors, alpha = alpha, beta = nuisance$beta,
       errors = errors, positive = nuisance$positive, negative = nuisance$negative,
       mu_f = nuisance$mu_f, distribution = distribution, rho = rho,
       delta = delta, seed = seed, replicate = replicate)
}

.rfs_empty_totals <- function(k) {
  dn <- list(NULL, c("LS", "SS", "IN"),
             c("FDR", "TDP", "discoveries", "false_discoveries", "true_discoveries"))
  pn <- list(NULL, c("SS-LS", "IN-LS", "IN-SS"), dn[[3L]])
  list(n = integer(k), sum = array(0, c(k, 3L, 5L), dn),
       sumsq = array(0, c(k, 3L, 5L), dn),
       failed = matrix(0L, k, 3L, dimnames = dn[1:2]),
       pair_sum = array(0, c(k, 3L, 5L), pn),
       pair_sumsq = array(0, c(k, 3L, 5L), pn),
       pair_successful = matrix(0L, k, 3L, dimnames = pn[1:2]))
}

.rfs_study_block <- function(task, cfg) {
  cases <- cfg$scenarios[task$indices, , drop = FALSE]
  totals <- .rfs_empty_totals(nrow(cases))
  T <- cases$T[1L]
  N <- cases$N[1L]
  distribution <- cases$distribution[1L]
  for (r in seq_along(task$nuisance_seeds)) {
    nuisance <- .rfs_draw_nuisance(T, N, task$nuisance_seeds[[r]])
    projection <- observed_projection(T, nuisance$factors)
    sources <- .rfs_draw_error_sources(T, N, distribution, task$error_seeds[[r]])
    factor_component <- nuisance$factors %*% t(nuisance$beta)
    for (rho in unique(cases$rho)) {
      Z0 <- project_observations(factor_component + .rfs_correlate_errors(sources, rho),
                                 projection)
      scale <- tryCatch(pqr_scale_projected(Z0, projection), error = identity)
      pilot_radius <- if (cfg$smoothing > 0 && !inherits(scale, "error"))
        tryCatch(.rfs_smooth_pilot(Z0, projection, scale$d)$median_radius,
                 error = identity) else NULL
      for (i in which(cases$rho == rho)) {
        alpha <- nuisance$alpha_base
        alpha[nuisance$positive] <- alpha[nuisance$positive] + cases$delta[i]
        Z <- Z0 + tcrossprod(projection$h, alpha)
        fits <- fit_three_projected(Z, projection, q = cfg$q, tol = cfg$tol,
                                    max_iter = cfg$max_iter, capture_errors = TRUE,
                                    pqr_scale = scale, smoothing = cfg$smoothing,
                                    in_pilot_radius = pilot_radius)
        values <- matrix(0, 3L, 5L)
        successful <- rep(TRUE, 3L)
        for (m in seq_len(3L)) {
          fit <- fits[[c("LS", "SS", "IN")[m]]]
          if (inherits(fit, "error")) {
            successful[m] <- FALSE
          } else {
            selected <- fit$selected
            if (!is.logical(selected) || length(selected) != N || anyNA(selected))
              stop("A fit returned invalid selection indicators.")
            R <- sum(selected)
            V <- sum(selected & alpha <= 0)
            S <- sum(selected & alpha > 0)
            values[m, ] <- c(V / max(R, 1), S / max(sum(alpha > 0), 1), R, V, S)
          }
        }
        differences <- rbind(values[2L, ] - values[1L, ],
                             values[3L, ] - values[1L, ],
                             values[3L, ] - values[2L, ])
        totals$n[i] <- totals$n[i] + 1L
        totals$sum[i, , ] <- totals$sum[i, , ] + values
        totals$sumsq[i, , ] <- totals$sumsq[i, , ] + values^2
        totals$failed[i, ] <- totals$failed[i, ] + !successful
        totals$pair_sum[i, , ] <- totals$pair_sum[i, , ] + differences
        totals$pair_sumsq[i, , ] <- totals$pair_sumsq[i, , ] + differences^2
        totals$pair_successful[i, ] <- totals$pair_successful[i, ] +
          c(successful[2L] && successful[1L], successful[3L] && successful[1L],
            successful[3L] && successful[2L])
      }
    }
  }
  list(indices = task$indices, totals = totals)
}

.rfs_add_totals <- function(total, block) {
  idx <- block$indices
  add <- block$totals
  total$n[idx] <- total$n[idx] + add$n
  for (name in c("sum", "sumsq", "pair_sum", "pair_sumsq"))
    total[[name]][idx, , ] <- total[[name]][idx, , , drop = FALSE] + add[[name]]
  for (name in c("failed", "pair_successful"))
    total[[name]][idx, ] <- total[[name]][idx, , drop = FALSE] + add[[name]]
  total
}

.rfs_mcse <- function(total, square_total, n) {
  if (n < 2L) return(rep(NA_real_, length(total)))
  sqrt(pmax(0, square_total - total^2 / n) / (n * (n - 1)))
}

.rfs_finalize_totals <- function(totals, scenarios, q) {
  summaries <- pairs <- list()
  metrics <- dimnames(totals$sum)[[3L]]
  methods <- dimnames(totals$sum)[[2L]]
  comparisons <- dimnames(totals$pair_sum)[[2L]]
  for (i in seq_len(nrow(scenarios))) {
    n <- totals$n[i]
    for (m in seq_along(methods)) {
      row <- scenarios[i, , drop = FALSE]
      row$method <- methods[m]
      row$repetitions <- n
      row$successful <- n - totals$failed[i, m]
      row$failed <- totals$failed[i, m]
      means <- totals$sum[i, m, ] / n
      se <- .rfs_mcse(totals$sum[i, m, ], totals$sumsq[i, m, ], n)
      for (j in seq_along(metrics)) {
        metric <- metrics[j]
        row[[metric]] <- means[j]
        row[[paste0(metric, "_MCSE")]] <- se[j]
        row[[paste0(metric, "_valid_only")]] <- if (row$successful > 0)
          totals$sum[i, m, j] / row$successful else NA_real_
        if (j <= 2L) {
          row[[paste0(metric, "_CI_low")]] <- max(0, means[j] - 1.96 * se[j])
          row[[paste0(metric, "_CI_high")]] <- min(1, means[j] + 1.96 * se[j])
          row[[paste0(metric, "_failure_bound_low")]] <- means[j]
          row[[paste0(metric, "_failure_bound_high")]] <-
            (totals$sum[i, m, j] + row$failed) / n
        }
      }
      row$FDR_above_nominal <- row$FDR > q
      row$FDR_lower_CI_above_nominal <- row$FDR_CI_low > q
      summaries[[length(summaries) + 1L]] <- row
    }
    for (p in seq_along(comparisons)) {
      means <- totals$pair_sum[i, p, ] / n
      se <- .rfs_mcse(totals$pair_sum[i, p, ], totals$pair_sumsq[i, p, ], n)
      for (j in seq_along(metrics)) {
        row <- scenarios[i, , drop = FALSE]
        row$comparison <- comparisons[p]
        row$metric <- metrics[j]
        row$attempted_pairs <- row$operational_pairs <- n
        row$successful_pairs <- totals$pair_successful[i, p]
        row$failed_pairs <- n - row$successful_pairs
        row$difference <- means[j]
        row$MCSE <- se[j]
        row$CI_low <- means[j] - 1.96 * se[j]
        row$CI_high <- means[j] + 1.96 * se[j]
        pairs[[length(pairs) + 1L]] <- row
      }
    }
  }
  summary <- do.call(rbind, summaries)
  paired <- do.call(rbind, pairs)
  rownames(summary) <- rownames(paired) <- NULL
  list(summary = summary, paired_differences = paired)
}

.rfs_export_source_dependencies <- function(cl, fn) {
  env <- environment(fn)
  todo <- ".rfs_study_block"
  found <- character()
  while (length(todo)) {
    name <- todo[1L]
    todo <- todo[-1L]
    if (name %in% found || !exists(name, envir = env, inherits = FALSE)) next
    found <- c(found, name)
    object <- get(name, envir = env, inherits = FALSE)
    if (is.function(object)) {
      globals <- codetools::findGlobals(object, merge = TRUE)
      todo <- unique(c(todo, setdiff(globals, found)))
    }
  }
  parallel::clusterExport(cl, found, envir = env)
  invisible(NULL)
}

simulate_study <- function(reps = 2000L, workers = 25L, seed = 20261105L,
                           progress = TRUE, q = 0.10, tol = 1e-8,
                           max_iter = 1000L, block_size = 5L,
                           scenarios = scenario_inventory(), smoothing = 0.5) {
  smoothing <- .rfs_smoothing(smoothing)
  reps <- .rfs_integer(reps, "reps")
  workers <- .rfs_integer(workers, "workers", 1L, 25L)
  seed <- .rfs_integer(seed, "seed", 0L)
  max_iter <- .rfs_integer(max_iter, "max_iter")
  block_size <- .rfs_integer(block_size, "block_size")
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress))
    stop("progress must be TRUE or FALSE.", call. = FALSE)
  if (length(q) != 1L || !is.finite(q) || q <= 0 || q >= 1)
    stop("q must be in (0, 1).", call. = FALSE)
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0)
    stop("tol must be positive and finite.", call. = FALSE)
  inventory <- scenario_inventory()
  if (!is.data.frame(scenarios) || !"case_id" %in% names(scenarios) ||
      nrow(scenarios) == 0L || anyDuplicated(scenarios$case_id) ||
      anyNA(scenarios$case_id) || !all(scenarios$case_id %in% inventory$case_id))
    stop("scenarios must identify a nonempty, unique subset of scenario_inventory().",
         call. = FALSE)
  # Use the canonical design for each ID; caller-supplied numerical columns
  # cannot silently change the experiment while retaining its name.
  scenarios <- inventory[match(scenarios$case_id, inventory$case_id), , drop = FALSE]
  rownames(scenarios) <- NULL
  old_rng <- .rfs_rng_snapshot()
  on.exit(.rfs_rng_restore(old_rng), add = TRUE)
  seeds <- .rfs_seed_table(seed, reps)
  groups <- unique(scenarios[c("distribution", "T")])
  tasks <- list()
  for (g in seq_len(nrow(groups))) {
    idx <- which(scenarios$distribution == groups$distribution[g] &
                   scenarios$T == groups$T[g])
    offset <- if (groups$T[g] == 60L) 0L else 4L
    code <- match(groups$distribution[g], c("normal", "t3", "mixture"))
    for (first in seq.int(1L, reps, by = block_size)) {
      rr <- seq.int(first, min(first + block_size - 1L, reps))
      tasks[[length(tasks) + 1L]] <- list(
        indices = idx, nuisance_seeds = seeds[[offset + 1L]][rr],
        error_seeds = seeds[[offset + code + 1L]][rr]
      )
    }
  }
  rm(seeds)
  actual_workers <- min(workers, length(tasks))
  cfg <- list(scenarios = scenarios, q = q, tol = tol, max_iter = max_iter,
              smoothing = smoothing)
  total <- .rfs_empty_totals(nrow(scenarios))
  start <- proc.time()[["elapsed"]]
  # Set limits before spawning so numerical libraries see them at startup.
  env_names <- c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                 "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")
  previous_env <- Sys.getenv(env_names, unset = NA_character_)
  do.call(Sys.setenv, as.list(stats::setNames(rep("1", length(env_names)), env_names)))
  on.exit({
    for (j in seq_along(env_names)) {
      if (is.na(previous_env[j])) Sys.unsetenv(env_names[j])
      else do.call(Sys.setenv, stats::setNames(list(previous_env[j]), env_names[j]))
    }
  }, add = TRUE)
  cl <- NULL
  if (actual_workers > 1L) {
    cl <- parallel::makePSOCKcluster(actual_workers, outfile = "")
    on.exit(parallel::stopCluster(cl), add = TRUE)
    if (isNamespace(environment(simulate_study))) {
      pkg <- getNamespaceName(environment(simulate_study))
      parallel::clusterCall(cl, function(paths, package) {
        .libPaths(paths)
        loadNamespace(package)
        NULL
      }, .libPaths(), pkg)
    } else {
      .rfs_export_source_dependencies(cl, .rfs_study_block)
    }
  }
  if (progress) {
    message(sprintf("%d scenarios x %d repetitions; %d %s; summaries only.",
                    nrow(scenarios), reps, actual_workers,
                    if (actual_workers > 1L) "PSOCK workers" else "serial worker"))
  }
  batch_size <- max(1L, 2L * actual_workers)
  for (first in seq.int(1L, length(tasks), by = batch_size)) {
    last <- min(first + batch_size - 1L, length(tasks))
    batch <- tasks[seq.int(first, last)]
    blocks <- if (is.null(cl)) lapply(batch, .rfs_study_block, cfg = cfg) else
      parallel::clusterApplyLB(cl, batch, .rfs_study_block, cfg = cfg)
    # Merge in task order, never arrival order. No replicate records survive.
    for (block in blocks) total <- .rfs_add_totals(total, block)
    rm(blocks, block)
    if (progress) message(sprintf("%d/%d blocks; %d/%d panels; %.1f seconds.",
                                  last, length(tasks), sum(total$n),
                                  nrow(scenarios) * reps,
                                  proc.time()[["elapsed"]] - start))
  }
  if (any(total$n != reps)) stop("Internal repetition-count audit failed.")
  answer <- .rfs_finalize_totals(total, scenarios, q)
  answer$scenarios <- scenarios
  answer$config <- list(
    repetitions = reps, seed = seed, workers_requested = workers,
    workers_used = actual_workers, backend = if (is.null(cl)) "serial" else "PSOCK",
    block_size = block_size, q = q, tol = tol, max_iter = max_iter,
    IN_smoothing = smoothing,
    IN_implementation = if (smoothing == 0) "legacy_inverse_norm" else "smoothed_inverse_norm",
    IN_lambda = "smoothing * median(PQR-standardized unrestricted OLS residual norms); fixed during fitting",
    IN_studentizer = if (smoothing == 0) "Original inverse-norm identity" else
      "General radial A=mean(h^2/(r^2+lambda^2)); B=mean(r^2*(h/(r^2+lambda^2)-P_F(h/(r^2+lambda^2)))^2)",
    N = 200L, p = 3L, unique_scenarios = nrow(scenarios),
    panels = nrow(scenarios) * reps, method_evaluations = nrow(scenarios) * reps * 3L,
    rng = "L'Ecuyer-CMRG / Inversion / Rejection",
    seed_rule = paste("Independent streams in order T=60 then 120, each with",
                      "nuisance, normal, t3, mixture; repetition uses its numbered substream."),
    reproduction = paste("R distributional reproduction of the Python design, not bitwise",
                         "reproduction of NumPy SeedSequence/PCG64 draws or historical results."),
    paired = paste("Identical panels across methods; common nuisance draws across distributions;",
                   "common errors across rho and delta; overlapping scenario simulated once."),
    failure_policy = paste("Failed fit makes no rejections. All attempted repetitions enter",
                           "operational means and paired differences; no replacement draws."),
    persistence = "No panels, repetition records, block files, checkpoints, or RNG states saved.",
    IN_calibration = "Residual degrees-of-freedom multiplier sqrt((T-rank([1,F]))/T)",
    SS_calibration = "Original joint spatial-sign statistic; no df correction",
    LS_calibration = "Original OLS statistic; residual variance divisor T",
    elapsed_seconds = proc.time()[["elapsed"]] - start,
    R_version = R.version.string
  )
  class(answer) <- c("radial_study", "list")
  answer
}
