# Figure data and performance measures follow the Monte Carlo section.
.paper_grid_facts <- function(study) {
  summary <- study$summary
  differences <- study$paired_differences
  ranges <- comparisons <- list()
  for (grid in c("signal_grid", "rho_sensitivity")) {
    cells <- summary[summary[[grid]], , drop = FALSE]
    pairs <- differences[differences$case_id %in% cells$case_id &
      differences$metric %in% c("FDR", "TDP") &
      differences$comparison %in% c("IN-LS", "IN-SS"), , drop = FALSE]
    for (group in split(cells, list(cells$distribution, cells$method), drop = TRUE)) {
      ranges[[length(ranges) + 1L]] <- data.frame(
        experiment = grid, distribution = group$distribution[1L],
        method = group$method[1L], settings = nrow(group),
        FDR_min = min(group$FDR), FDR_max = max(group$FDR),
        TDP_min = min(group$TDP), TDP_max = max(group$TDP))
    }
    for (group in split(pairs, list(pairs$distribution, pairs$comparison, pairs$metric),
                        drop = TRUE)) {
      comparisons[[length(comparisons) + 1L]] <- data.frame(
        experiment = grid, distribution = group$distribution[1L],
        comparison = group$comparison[1L], metric = group$metric[1L],
        settings = nrow(group), min_difference = min(group$difference),
        max_difference = max(group$difference),
        intervals_above_zero = sum(group$CI_low > 0),
        intervals_below_zero = sum(group$CI_high < 0))
    }
  }
  inverse <- summary[summary$method == "IN", , drop = FALSE]
  list(ranges = do.call(rbind, ranges), paired_ranges = do.call(rbind, comparisons),
       maximum_IN_FDR = inverse[which.max(inverse$FDR), , drop = FALSE],
       IN_above_nominal = sum(inverse$FDR > study$config$q),
       IN_lower_interval_above_nominal = sum(inverse$FDR_CI_low > study$config$q))
}

# Put the two components of the smoothing contrast on the paper's
# denominator (all attempted panels); conditional means can be missing
# where no panel belongs to the corresponding group.
.paper_smoothing_facts <- function(study) {
  paired_tdp <- study$paired[study$paired$metric == "TDP" &
    study$paired$smoothing == 0.5, , drop = FALSE]
  component <- function(n, value) {
    mean(ifelse(n > 0, n * value / paired_tdp$repetitions, 0))
  }
  list(overall = study$group_summary[study$group_summary$group == "Overall", ],
       overall_paired = study$group_paired[study$group_paired$group == "Overall", ],
       rescued_TDP_contribution = component(paired_tdp$rescued, paired_tdp$rescued_smooth_mean),
       both_success_TDP_contribution = component(paired_tdp$both_successful,
                                                 paired_tdp$both_success_difference),
       largest_TDP_loss = paired_tdp[which.min(paired_tdp$difference), , drop = FALSE])
}

reproduce_simulations <- function(reps = 2000L, workers = 1L, output_dir = NULL,
                                  progress = TRUE) {
  study <- simulate_study(reps = reps, workers = workers, seed = 20261105L,
                          progress = progress, smoothing = 0.5)
  smooth <- simulate_smooth_study(reps = reps, workers = workers, seed = 20261105L,
                                  progress = progress, c_grid = 0.5)
  figures <- if (!is.null(output_dir)) plot_simulation_figures(study, output_dir)
  list(study = study, smoothing = smooth, simulation_facts = .paper_grid_facts(study),
       smoothing_facts = .paper_smoothing_facts(smooth), figures = figures)
}

.paper_draw_grid <- function(data, xvar, fdr_top, nominal, path) {
  grDevices::pdf(path, width = 7, height = 6.8, family = "Times", useDingbats = FALSE)
  on.exit(grDevices::dev.off())
  graphics::par(mfrow = c(3, 2), mar = c(3.2, 4.3, 1.2, 0.6),
                oma = c(0.1, 0, 2.2, 0), mgp = c(2.4, 0.55, 0), tcl = -0.25,
                xaxs = "i", yaxs = "i", cex.axis = 0.9, cex.lab = 1)
  methods <- c("LS", "SS", "IN")
  colours <- c("#000000", "#1F77B4", "#D62728")
  styles <- c(2L, 1L, 4L)
  symbols <- c(18L, 16L, 17L)
  row_titles <- c(normal = "Gaussian", t3 = "Student t3", mixture = "Mixture normal")
  for (distribution in names(row_titles)) {
    row <- data[data$distribution == distribution, , drop = FALSE]
    for (metric in c("FDR", "TDP")) {
      ytop <- if (metric == "FDR") fdr_top else 1.02
      yticks <- if (metric == "FDR") pretty(c(0, ytop), n = 4) else seq(0, 1, 0.25)
      yticks <- yticks[yticks >= 0 & yticks <= ytop + 1e-9]
      xlim <- range(row[[xvar]])
      graphics::plot(NA, xlim = xlim + c(-1, 1) * diff(xlim) * .025,
        ylim = c(0, ytop), axes = FALSE,
        xlab = if (xvar == "delta") "Signal increment" else "Residual correlation",
        ylab = if (metric == "FDR") row_titles[[distribution]] else "")
      if (distribution == "normal") graphics::title(main =
          if (metric == "FDR") "False discovery rate" else "True discovery proportion",
          font.main = 1, cex.main = 1.05, line = 0.3)
      graphics::abline(h = yticks, col = "#DBDEE3", lwd = 0.6)
      graphics::axis(1); graphics::axis(2, at = yticks,
        labels = paste0(formatC(100 * yticks, format = "fg", digits = 3), "%"), las = 1)
      graphics::box(bty = "l")
      for (i in seq_along(methods)) {
        line <- row[row$method == methods[i], , drop = FALSE]
        line <- line[order(line[[xvar]]), , drop = FALSE]
        x <- line[[xvar]]
        lo <- line[[paste0(metric, "_CI_low")]]
        hi <- line[[paste0(metric, "_CI_high")]]
        if (all(is.finite(c(lo, hi)))) graphics::polygon(c(x, rev(x)), c(lo, rev(hi)),
          col = grDevices::adjustcolor(colours[i], 0.12), border = NA)
      }
      if (metric == "FDR") graphics::abline(h = nominal, lty = 3)
      for (i in seq_along(methods)) {
        line <- row[row$method == methods[i], , drop = FALSE]
        line <- line[order(line[[xvar]]), , drop = FALSE]
        graphics::lines(line[[xvar]], line[[metric]], type = "o",
                        col = colours[i], lty = styles[i], pch = symbols[i], cex = .55)
      }
    }
  }
  graphics::par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), new = TRUE)
  graphics::plot.new()
  graphics::legend("top", legend = paste0(methods, "-BH"), col = colours,
                   lty = styles, pch = symbols, horiz = TRUE, bty = "n", cex = .9)
}

plot_simulation_figures <- function(study, output_dir) {
  if (!is.character(output_dir) || length(output_dir) != 1L ||
      is.na(output_dir) || !nzchar(output_dir)) stop("Supply a PDF output directory.")
  cells <- study$summary
  if (!is.data.frame(cells) || !all(c("case_id", "method", "FDR", "TDP", "signal_grid",
      "rho_sensitivity", "FDR_CI_high", "FDR_CI_low", "TDP_CI_high", "TDP_CI_low",
      "distribution", "T", "rho", "delta") %in% names(cells)))
    stop("Supply the result of simulate_study().")
  keys <- expand.grid(case_id = scenario_inventory()$case_id,
                      method = c("LS", "SS", "IN"), stringsAsFactors = FALSE)
  if (nrow(cells) != nrow(keys) || anyDuplicated(cells[c("case_id", "method")]) ||
      nrow(merge(keys, cells[c("case_id", "method")])) != nrow(keys))
    stop("The three paper figures require all 78 scenarios and all three methods.")
  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE))
    stop("Cannot create the PDF directory.")
  limit <- function(x) max(.16, ceiling((max(c(x$FDR, x$FDR_CI_high), na.rm = TRUE) + .01) / .04) * .04)
  signal <- cells[cells$signal_grid, , drop = FALSE]
  correlation <- cells[cells$rho_sensitivity, , drop = FALSE]
  paths <- file.path(output_dir, c("signal_T60.pdf", "signal_T120.pdf", "correlation.pdf"))
  for (i in 1:2) .paper_draw_grid(signal[signal$T == c(60L, 120L)[i], ],
                                  "delta", limit(signal), study$config$q, paths[i])
  .paper_draw_grid(correlation, "rho", limit(correlation), study$config$q, paths[3L])
  invisible(paths)
}
