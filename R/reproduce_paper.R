reproduce_paper <- function(data, rf_file = NULL, output_dir = "paper-results",
                            reps = 2000L, workers = 1L, progress = TRUE) {
  if (!is.character(output_dir) || length(output_dir) != 1L ||
      is.na(output_dir) || !nzchar(output_dir)) {
    stop("Supply a directory for the four paper figures.", call. = FALSE)
  }
  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) {
    stop("Cannot create the figure directory.", call. = FALSE)
  }
  empirical <- reproduce_empirical(data, rf_file)
  simulations <- reproduce_simulations(reps = reps, workers = workers,
                                       output_dir = output_dir, progress = progress)
  empirical_pdf <- file.path(output_dir, "empirical_T69.pdf")
  plot_empirical_wealth(empirical, empirical_pdf)
  list(simulations = simulations, empirical = empirical,
       figures = c(simulations$figures, empirical_pdf))
}
