# Paired smooth-IN experiment. All panels, seeds and replication/block
# accumulators remain in memory. Only final summaries leave this function.

.sin_methods <- function(c_grid) c("IN", paste0("sIN_c", format(c_grid, trim = TRUE, scientific = FALSE)))

.sin_groups <- function(cases) {
  names <- c("Overall", paste0("distribution=", unique(cases$distribution)),
             paste0("T=", sort(unique(cases$T))),
             paste0("distribution=", unique(cases[c("distribution", "T")])$distribution,
                    ";T=", unique(cases[c("distribution", "T")])$T))
  membership <- matrix(FALSE, nrow(cases), length(names), dimnames = list(NULL, names))
  membership[, 1L] <- TRUE
  for (g in seq_along(names)[-1L]) {
    name <- names[g]
    if (grepl(";T=", name, fixed = TRUE)) {
      parts <- strsplit(name, ";T=", fixed = TRUE)[[1L]]
      membership[, g] <- cases$distribution == sub("distribution=", "", parts[1L], fixed = TRUE) &
        cases$T == as.integer(parts[2L])
    } else if (startsWith(name, "distribution=")) {
      membership[, g] <- cases$distribution == sub("distribution=", "", name, fixed = TRUE)
    } else membership[, g] <- cases$T == as.integer(sub("T=", "", name, fixed = TRUE))
  }
  list(names = names, membership = membership, counts = colSums(membership))
}

.sin_empty <- function(k, methods) {
  m <- length(methods)
  metrics <- c("failure_rate", "FDR", "TDP")
  dn <- list(NULL, methods, metrics)
  pn <- list(NULL, methods[-1L], metrics)
  reasons <- c("radius_collapse", "iteration_limit", "scale_failure", "invalid_numeric", "other_fit_failure")
  list(n = integer(k), sum = array(0, c(k, m, 3L), dn),
       sumsq = array(0, c(k, m, 3L), dn),
       pair_sum = array(0, c(k, m - 1L, 3L), pn),
       pair_sumsq = array(0, c(k, m - 1L, 3L), pn),
       both_sum = array(0, c(k, m - 1L, 3L), pn),
       both_sumsq = array(0, c(k, m - 1L, 3L), pn),
       rescued_sum = array(0, c(k, m - 1L, 3L), pn),
       both_n = matrix(0L, k, m - 1L), rescued = matrix(0L, k, m - 1L),
       regressed = matrix(0L, k, m - 1L),
       reasons = array(0L, c(k, m, length(reasons)), list(NULL, methods, reasons)),
       successful = matrix(0L, k, m), high_weight = matrix(0L, k, m),
       near_center = matrix(0L, k, m), weight_sum = matrix(0, k, m),
       lambda_sum = matrix(0, k, m), iterations_sum = matrix(0, k, m),
       iterations_max = matrix(0, k, m), score_max = matrix(0, k, m),
       radius_ratio_min = matrix(Inf, k, m))
}

.sin_reason <- function(error) {
  message <- conditionMessage(error)
  if (grepl("PQR|scale|quantile", message, ignore.case = TRUE)) return("scale_failure")
  if (grepl("collaps|radius", message, ignore.case = TRUE)) return("radius_collapse")
  if (grepl("converg|iteration", message, ignore.case = TRUE)) return("iteration_limit")
  if (grepl("nonfinite|nonpositive|curvature|variance|numeric", message, ignore.case = TRUE)) return("invalid_numeric")
  "other_fit_failure"
}

.sin_capture <- function(expr) tryCatch(force(expr), radial_fit_error = identity)

.sin_block <- function(task, cfg) {
  cases <- cfg$scenarios[task$indices, , drop = FALSE]
  methods <- cfg$methods
  m <- length(methods)
  total <- .sin_empty(nrow(cases), methods)
  rep_groups <- array(0, c(length(task$replicates), m, length(cfg$groups$names), 3L))
  for (j in seq_along(task$replicates)) {
    T <- cases$T[1L]; N <- cases$N[1L]; distribution <- cases$distribution[1L]
    nuisance <- .rfs_draw_nuisance(T, N, task$nuisance_seeds[[j]])
    projection <- observed_projection(T, nuisance$factors)
    sources <- .rfs_draw_error_sources(T, N, distribution, task$error_seeds[[j]])
    factor_component <- nuisance$factors %*% t(nuisance$beta)
    for (rho in unique(cases$rho)) {
      Z0 <- project_observations(factor_component + .rfs_correlate_errors(sources, rho), projection)
      scale <- .sin_capture(pqr_scale_projected(Z0, projection))
      pilot <- if (inherits(scale, "error")) scale else .sin_capture(.rfs_smooth_pilot(Z0, projection, scale$d))
      for (i in which(cases$rho == rho)) {
        alpha <- nuisance$alpha_base
        alpha[nuisance$positive] <- alpha[nuisance$positive] + cases$delta[i]
        Z <- Z0 + tcrossprod(projection$h, alpha)
        fits <- vector("list", m)
        lambdas <- numeric(m)
        if (inherits(scale, "error")) {
          fits[] <- list(scale)
        } else {
          fits[[1L]] <- .sin_capture(fit_in_projected(Z, projection, scale$d,
                    q = cfg$q, tol = cfg$tol, max_iter = cfg$max_iter))
          for (s in seq_along(cfg$c_grid)) {
            if (inherits(pilot, "error")) fits[[s + 1L]] <- pilot else {
              lambdas[s + 1L] <- cfg$c_grid[s] * pilot$median_radius
              fits[[s + 1L]] <- .sin_capture(.rfs_fit_smooth_in(Z, projection, scale$d,
                lambda = lambdas[s + 1L], q = cfg$q, tol = cfg$tol, max_iter = cfg$max_iter))
            }
          }
        }
        values <- matrix(0, m, 3L)
        successful <- logical(m)
        for (s in seq_len(m)) {
          fit <- fits[[s]]
          if (inherits(fit, "error")) {
            values[s, 1L] <- 1
            reason <- .sin_reason(fit)
            total$reasons[i, s, reason] <- total$reasons[i, s, reason] + 1L
          } else {
            successful[s] <- TRUE
            selected <- fit$selected
            if (!is.logical(selected) || length(selected) != N || anyNA(selected))
              stop("A successful fit returned invalid selection indicators.")
            R <- sum(selected); V <- sum(selected & alpha <= 0); S <- sum(selected & alpha > 0)
            values[s, ] <- c(0, V / max(R, 1), S / sum(alpha > 0))
            # The original baseline fit does not store its largest date weight.
            if (s == 1L) {
              residual <- sweep(Z - tcrossprod(projection$h, fit$theta), 2L, fit$d, "/")
              r2 <- rowSums(residual * residual)
              weight <- projection$h^2 / r2
              max_weight <- max(weight) / sum(weight)
            } else max_weight <- fit$max_date_weight
            if (length(max_weight) != 1L || !is.finite(max_weight)) stop("Invalid fitted date-weight diagnostic.")
            total$successful[i, s] <- total$successful[i, s] + 1L
            total$high_weight[i, s] <- total$high_weight[i, s] + (max_weight > 0.5)
            total$near_center[i, s] <- total$near_center[i, s] + (fit$min_radius_ratio < 0.01)
            total$weight_sum[i, s] <- total$weight_sum[i, s] + max_weight
            total$iterations_sum[i, s] <- total$iterations_sum[i, s] + fit$iterations
            total$iterations_max[i, s] <- max(total$iterations_max[i, s], fit$iterations)
            total$score_max[i, s] <- max(total$score_max[i, s], fit$location_residual)
            total$radius_ratio_min[i, s] <- min(total$radius_ratio_min[i, s], fit$min_radius_ratio)
          }
        }
        differences <- sweep(values[-1L, , drop = FALSE], 2L, values[1L, ], "-")
        both <- successful[1L] & successful[-1L]
        rescue <- !successful[1L] & successful[-1L]
        regress <- successful[1L] & !successful[-1L]
        total$n[i] <- total$n[i] + 1L
        total$sum[i, , ] <- total$sum[i, , ] + values
        total$sumsq[i, , ] <- total$sumsq[i, , ] + values^2
        total$pair_sum[i, , ] <- total$pair_sum[i, , ] + differences
        total$pair_sumsq[i, , ] <- total$pair_sumsq[i, , ] + differences^2
        total$both_sum[i, , ] <- total$both_sum[i, , ] + differences * both
        total$both_sumsq[i, , ] <- total$both_sumsq[i, , ] + differences^2 * both
        total$rescued_sum[i, , ] <- total$rescued_sum[i, , ] + values[-1L, , drop = FALSE] * rescue
        total$both_n[i, ] <- total$both_n[i, ] + both
        total$rescued[i, ] <- total$rescued[i, ] + rescue
        total$regressed[i, ] <- total$regressed[i, ] + regress
        total$lambda_sum[i, ] <- total$lambda_sum[i, ] + lambdas
        for (g in which(cfg$groups$membership[task$indices[i], ]))
          rep_groups[j, , g, ] <- rep_groups[j, , g, ] + values
      }
    }
  }
  list(indices = task$indices, replicates = task$replicates,
       total = total, rep_groups = rep_groups)
}

.sin_add <- function(total, part) {
  idx <- part$indices; add <- part$total
  for (name in names(total)) {
    dims <- dim(total[[name]])
    if (is.null(dims)) total[[name]][idx] <- total[[name]][idx] + add[[name]]
    else if (length(dims) == 3L) total[[name]][idx, , ] <- total[[name]][idx, , , drop = FALSE] + add[[name]]
    else if (name %in% c("iterations_max", "score_max"))
      total[[name]][idx, ] <- pmax(total[[name]][idx, , drop = FALSE], add[[name]])
    else if (name == "radius_ratio_min")
      total[[name]][idx, ] <- pmin(total[[name]][idx, , drop = FALSE], add[[name]])
    else total[[name]][idx, ] <- total[[name]][idx, , drop = FALSE] + add[[name]]
  }
  total
}

.sin_finalize <- function(total, cases, methods, c_grid, rep_groups, groups, reps, q) {
  summaries <- paired <- reasons <- list()
  metrics <- dimnames(total$sum)[[3L]]
  for (i in seq_len(nrow(cases))) {
    n <- total$n[i]
    for (m in seq_along(methods)) {
      row <- cases[i, , drop = FALSE]
      row$method <- methods[m]; row$c <- row$smoothing <- c(0, c_grid)[m]
      row$repetitions <- n; row$successful <- total$successful[i, m]
      row$failed <- n - row$successful
      for (j in seq_along(metrics)) {
        metric <- metrics[j]
        row[[metric]] <- total$sum[i, m, j] / n
        se <- .rfs_mcse(total$sum[i, m, j], total$sumsq[i, m, j], n)
        row[[paste0(metric, "_MCSE")]] <- se
        row[[paste0(metric, "_CI_low")]] <- max(0, row[[metric]] - 1.96 * se)
        row[[paste0(metric, "_CI_high")]] <- min(1, row[[metric]] + 1.96 * se)
      }
      row$FDR_above_nominal <- row$FDR > q
      row$FDR_lower_CI_above_nominal <- row$FDR_CI_low > q
      row$lambda_mean <- total$lambda_sum[i, m] / n
      row$high_date_weight_count <- total$high_weight[i, m]
      row$near_center_count <- total$near_center[i, m]
      row$high_date_weight_fraction <- if (row$successful) total$high_weight[i, m] / row$successful else NA_real_
      row$near_center_fraction <- if (row$successful) total$near_center[i, m] / row$successful else NA_real_
      row$max_date_weight_mean <- if (row$successful) total$weight_sum[i, m] / row$successful else NA_real_
      row$iterations_mean <- if (row$successful) total$iterations_sum[i, m] / row$successful else NA_real_
      row$iterations_max <- total$iterations_max[i, m]
      row$score_max <- if (row$successful) total$score_max[i, m] else NA_real_
      row$min_radius_ratio <- if (row$successful) total$radius_ratio_min[i, m] else NA_real_
      summaries[[length(summaries) + 1L]] <- row
      for (reason in dimnames(total$reasons)[[3L]]) {
        reasons[[length(reasons) + 1L]] <- data.frame(case_id = cases$case_id[i], method = methods[m],
          reason = reason, count = total$reasons[i, m, reason], stringsAsFactors = FALSE)
      }
    }
    for (m in seq_along(c_grid)) for (j in seq_along(metrics)) {
      row <- cases[i, , drop = FALSE]
      row$method <- methods[m + 1L]; row$c <- row$smoothing <- c_grid[m]; row$metric <- metrics[j]
      row$comparison <- paste0(methods[m + 1L], "-IN"); row$repetitions <- n
      row$difference <- total$pair_sum[i, m, j] / n
      row$MCSE <- .rfs_mcse(total$pair_sum[i, m, j], total$pair_sumsq[i, m, j], n)
      row$CI_low <- row$difference - 1.96 * row$MCSE
      row$CI_high <- row$difference + 1.96 * row$MCSE
      row$rescued <- total$rescued[i, m]; row$regressed <- total$regressed[i, m]
      row$both_successful <- total$both_n[i, m]
      row$both_success_difference <- if (row$both_successful)
        total$both_sum[i, m, j] / row$both_successful else NA_real_
      row$both_success_MCSE <- .rfs_mcse(total$both_sum[i, m, j], total$both_sumsq[i, m, j], row$both_successful)
      row$rescued_smooth_mean <- if (row$rescued) total$rescued_sum[i, m, j] / row$rescued else NA_real_
      paired[[length(paired) + 1L]] <- row
    }
  }
  group_summary <- group_paired <- list()
  for (g in seq_along(groups$names)) {
    for (m in seq_along(methods)) {
      values <- matrix(rep_groups[, m, g, ], reps, 3L) / groups$counts[g]
      row <- data.frame(group = groups$names[g], scenarios = groups$counts[g],
        method = methods[m], c = c(0, c_grid)[m], smoothing = c(0, c_grid)[m], repetitions = reps, independent_replicates = reps,
        attempted_panels = reps * groups$counts[g], stringsAsFactors = FALSE)
      for (j in seq_along(metrics)) {
        row[[metrics[j]]] <- mean(values[, j])
        row[[paste0(metrics[j], "_MCSE")]] <- if (reps > 1L) stats::sd(values[, j]) / sqrt(reps) else NA_real_
      }
      row$failed <- as.integer(round(sum(rep_groups[, m, g, 1L])))
      group_summary[[length(group_summary) + 1L]] <- row
    }
    for (m in seq_along(c_grid)) for (j in seq_along(metrics)) {
      difference <- (rep_groups[, m + 1L, g, j] - rep_groups[, 1L, g, j]) / groups$counts[g]
      se <- if (reps > 1L) stats::sd(difference) / sqrt(reps) else NA_real_
      group_paired[[length(group_paired) + 1L]] <- data.frame(group = groups$names[g],
        scenarios = groups$counts[g], method = methods[m + 1L], c = c_grid[m], smoothing = c_grid[m], metric = metrics[j],
        comparison = paste0(methods[m + 1L], "-IN"), independent_replicates = reps,
        difference = mean(difference), MCSE = se,
        CI_low = mean(difference) - 1.96 * se, CI_high = mean(difference) + 1.96 * se,
        stringsAsFactors = FALSE)
    }
  }
  answer <- lapply(list(summary = summaries, paired = paired, group_summary = group_summary,
                       group_paired = group_paired, failure_reasons = reasons), function(x) {
    out <- do.call(rbind, x); rownames(out) <- NULL; out
  })
  answer
}

.sin_export <- function(cl) {
  env <- environment(.sin_block); todo <- ".sin_block"; found <- character()
  while (length(todo)) {
    name <- todo[1L]; todo <- todo[-1L]
    if (name %in% found || !exists(name, envir = env, inherits = FALSE)) next
    found <- c(found, name); object <- get(name, envir = env, inherits = FALSE)
    if (is.function(object)) todo <- unique(c(todo, setdiff(codetools::findGlobals(object, merge = TRUE), found)))
  }
  parallel::clusterExport(cl, found, envir = env)
}

simulate_smooth_study <- function(reps = 2000L, workers = 1L, seed = 20261105L,
                                 c_grid = 0.5, q = .1, tol = 1e-8,
                                 max_iter = 1000L, block_size = 5L,
                                 scenarios = scenario_inventory(), progress = TRUE) {
  reps <- .rfs_integer(reps, "reps"); workers <- .rfs_integer(workers, "workers", 1L, 25L)
  seed <- .rfs_integer(seed, "seed", 0L); block_size <- .rfs_integer(block_size, "block_size")
  .fit_controls(q, tol, max_iter)
  if (!is.numeric(c_grid) || !length(c_grid) || any(!is.finite(c_grid)) ||
      any(c_grid <= 0) || anyDuplicated(c_grid)) stop("c_grid must contain distinct finite positive values.")
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) stop("progress must be TRUE or FALSE.")
  inventory <- scenario_inventory()
  if (!is.data.frame(scenarios) || !"case_id" %in% names(scenarios) || !nrow(scenarios) ||
      anyNA(scenarios$case_id) || anyDuplicated(scenarios$case_id) ||
      !all(scenarios$case_id %in% inventory$case_id)) stop("Invalid scenario subset.")
  scenarios <- inventory[match(scenarios$case_id, inventory$case_id), , drop = FALSE]
  rownames(scenarios) <- NULL
  methods <- .sin_methods(c_grid)
  groups <- .sin_groups(scenarios)
  rng_state <- .rfs_rng_snapshot(); on.exit(.rfs_rng_restore(rng_state), add = TRUE)
  seeds <- .rfs_seed_table(seed, reps)
  task_groups <- unique(scenarios[c("distribution", "T")]); tasks <- list()
  for (g in seq_len(nrow(task_groups))) {
    idx <- which(scenarios$distribution == task_groups$distribution[g] & scenarios$T == task_groups$T[g])
    offset <- if (task_groups$T[g] == 60L) 0L else 4L
    code <- match(task_groups$distribution[g], c("normal", "t3", "mixture"))
    for (first in seq.int(1L, reps, by = block_size)) {
      rr <- seq.int(first, min(first + block_size - 1L, reps))
      tasks[[length(tasks) + 1L]] <- list(indices = idx, replicates = rr,
        nuisance_seeds = seeds[[offset + 1L]][rr], error_seeds = seeds[[offset + code + 1L]][rr])
    }
  }
  rm(seeds)
  actual_workers <- min(workers, length(tasks))
  cfg <- list(scenarios = scenarios, methods = methods, groups = groups, c_grid = c_grid,
              q = q, tol = tol, max_iter = max_iter)
  total <- .sin_empty(nrow(scenarios), methods)
  rep_groups <- array(0, c(reps, length(methods), length(groups$names), 3L))
  env_names <- c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")
  previous_env <- Sys.getenv(env_names, unset = NA_character_)
  do.call(Sys.setenv, as.list(stats::setNames(rep("1", length(env_names)), env_names)))
  on.exit({for (j in seq_along(env_names)) {
    if (is.na(previous_env[j])) Sys.unsetenv(env_names[j])
    else do.call(Sys.setenv, stats::setNames(list(previous_env[j]), env_names[j]))
  }}, add = TRUE)
  start <- proc.time()[["elapsed"]]; cl <- NULL
  if (actual_workers > 1L) {
    cl <- parallel::makePSOCKcluster(actual_workers, outfile = "")
    on.exit(parallel::stopCluster(cl), add = TRUE)
    if (isNamespace(environment(simulate_smooth_study))) {
      parallel::clusterCall(cl, function(paths, package) {
        .libPaths(paths)
        loadNamespace(package)
        NULL
      }, .libPaths(), getNamespaceName(environment(simulate_smooth_study)))
    } else {
      .sin_export(cl)
    }
  }
  if (progress) message(sprintf("Smooth IN: %d scenarios x %d repetitions x %d methods; %d workers; memory-only accumulation.",
     nrow(scenarios), reps, length(methods), actual_workers))
  batch_size <- max(1L, 2L * actual_workers)
  for (first in seq.int(1L, length(tasks), by = batch_size)) {
    last <- min(first + batch_size - 1L, length(tasks))
    batch <- tasks[seq.int(first, last)]
    blocks <- if (is.null(cl)) lapply(batch, .sin_block, cfg = cfg) else
      parallel::clusterApplyLB(cl, batch, .sin_block, cfg = cfg)
    for (part in blocks) {
      total <- .sin_add(total, part)
      rr <- part$replicates
      rep_groups[rr, , , ] <- rep_groups[rr, , , , drop = FALSE] + part$rep_groups
    }
    rm(blocks, part)
    if (progress) message(sprintf("%d/%d blocks; %d/%d panels; %.1f seconds.", last, length(tasks),
      sum(total$n), nrow(scenarios) * reps, proc.time()[["elapsed"]] - start))
  }
  stopifnot(all(total$n == reps))
  answer <- .sin_finalize(total, scenarios, methods, c_grid, rep_groups, groups, reps, q)
  # No replication-level accumulator is returned or persisted.
  rm(rep_groups, total)
  answer$scenarios <- scenarios
  answer$config <- list(repetitions = reps, workers_requested = workers, workers_used = actual_workers,
    seed = seed, c_grid = c_grid, methods = methods, q = q, tol = tol, max_iter = max_iter,
    block_size = block_size, scenarios = nrow(scenarios), panels = nrow(scenarios) * reps,
    method_evaluations = nrow(scenarios) * reps * length(methods),
    rng = "Original R L'Ecuyer-CMRG streams, identical to the previous R reproduction.",
    smoothing = "lambda=c*median(||D_Q^-1/2 M_[1,F]Y||); lambda fixed throughout each alpha fit",
    score = "mean(h*e/(r^2+lambda^2)); original normalized raw-score stopping criterion",
    standard_error = "A_h=mean(h^2/(r^2+lambda^2)); pi=P_F(h/(r^2+lambda^2)); B_c=mean((r*(h/(r^2+lambda^2)-pi))^2)",
    df_calibration = "sqrt((T-p-1)/T), identical for all fitted statistics",
    baseline = "Unchanged original fit_in_projected, no smoothing or alternate initialization",
    failure_policy = "Failed fit selects no funds; FDP=TDP=0 and failure=1; every attempted panel retained",
    group_MCSE = "First average settings within each repetition, then calculate MCSE across independent repetitions; cross-setting CRN retained",
    diagnostics = "high_date_weight=max(h^2/(r^2+lambda^2))/sum(h^2/(r^2+lambda^2))>0.5; near_center=min(r)/median(r)<0.01; fractions conditional on successful fits",
    persistence = "No simulated panels, replication records, chunk files, checkpoints, or RNG states written to disk or returned",
    elapsed_seconds = proc.time()[["elapsed"]] - start, R_version = R.version.string)
  answer
}
