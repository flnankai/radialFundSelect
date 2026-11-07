# The Gao-Peng-Yan study archive and French risk-free observations are
# external research inputs. Only computations and aggregate outputs are here.
.paper_rf <- function(file, months) {
  if (is.null(file) || length(file) != 1L || !file.exists(file))
    stop("Supply the monthly French risk-free file (CSV or ZIP).", call. = FALSE)
  if (grepl("\\.zip$", file, ignore.case = TRUE)) {
    entries <- utils::unzip(file, list = TRUE)$Name
    csv <- entries[grepl("(^|/)F-F_Research_Data_Factors\\.CSV$", entries,
                         ignore.case = TRUE)]
    if (length(csv) != 1L)
      stop("Expected the monthly F-F_Research_Data_Factors CSV in the ZIP.",
           call. = FALSE)
    con <- unz(file, csv, open = "rb")
    on.exit(close(con))
    lines <- readLines(con, warn = FALSE)
  } else {
    lines <- readLines(file, warn = FALSE)
  }
  lines <- trimws(lines[nzchar(trimws(lines))])
  header <- which(toupper(gsub("[[:space:]]", "", lines)) == ",MKT-RF,SMB,HML,RF")
  if (length(header)) {
    rows <- lines[grepl("^[0-9]{6}[[:space:]]*,", lines)]
    table <- utils::read.csv(text = paste(rows, collapse = "\n"), header = FALSE,
                             colClasses = "character", strip.white = TRUE)
    ym <- format(months, "%Y%m")
    if (ncol(table) != 5L || anyDuplicated(table[[1L]]) ||
        !all(ym %in% table[[1L]]))
      stop("French monthly rows must contain each requested YYYYMM exactly once.",
           call. = FALSE)
    rf <- suppressWarnings(as.numeric(table[[5L]][match(ym, table[[1L]])]))
  } else {
    rf <- suppressWarnings(as.numeric(lines))
  }
  if (length(rf) != length(months) || any(!is.finite(rf)) ||
      any(rf %in% c(-99.99, -999)) || any(rf <= -100))
    stop("Require 368 monthly RF observations in percent, February 1987--September 2017.",
         call. = FALSE)
  rf / 100
}

.paper_fund_inputs <- function(data, rf_file) {
  if (length(data) != 1L || !is.character(data) || !file.exists(data))
    stop("Supply CaseStudy/mutual_fund.mat or the directory with the aligned CSVs.",
         call. = FALSE)
  months <- seq(as.Date("1987-03-01"), by = "month", length.out = 368L) - 1L
  if (dir.exists(data)) {
    funds <- t(as.matrix(utils::read.csv(file.path(data, "fund.csv"),
                                         header = FALSE)))
    factors <- as.matrix(utils::read.csv(file.path(data, "factor.csv"),
                                         header = FALSE))
    if (is.null(rf_file)) rf_file <- file.path(data, "rf.csv")
  } else {
    if (!requireNamespace("R.matlab", quietly = TRUE))
      stop("Install the optional R.matlab package to read mutual_fund.mat.",
           call. = FALSE)
    panel <- R.matlab::readMat(data)
    if (!identical(dim(panel$R), c(668L, 726L)) ||
        !identical(dim(panel$F), c(668L, 4L)))
      stop("Expected original archive matrices R (668 x 726) and F (668 x 4).",
           call. = FALSE)
    funds <- panel$R[301:668, , drop = FALSE]
    factors <- panel$F[301:668, , drop = FALSE]
  }
  if (!identical(dim(funds), c(368L, 726L)) ||
      !identical(dim(factors), c(368L, 4L)) ||
      !is.numeric(funds) || !is.numeric(factors) ||
      any(!is.finite(funds)) || any(!is.finite(factors)) || any(funds < -1))
    stop("Expected 368 aligned months, 726 decimal fund series and four decimal factors.",
         call. = FALSE)
  colnames(factors) <- c("MKT", "SMB", "HML", "MOM")
  list(funds = funds, factors = factors, rf = .paper_rf(rf_file, months), dates = months)
}

.paper_holding_wealth <- function(returns, cash) {
  # Each column starts with one dollar. Equal initial allocations drift.
  if (!ncol(returns)) {
    value <- cumprod(1 + cash)
    return(list(wealth = value, stress = value, missing = 0L))
  }
  values <- stress <- returns
  missing <- 0L
  for (j in seq_len(ncol(returns))) {
    r <- returns[, j]
    first <- match(0, r)
    if (!is.na(first)) {
      missed <- first:length(r)
      r[missed] <- cash[missed]
      missing <- missing + 1L
    }
    values[, j] <- cumprod(1 + r)
    stress[, j] <- values[, j]
    if (!is.na(first)) stress[missed, j] <- 0
  }
  list(wealth = rowMeans(values), stress = rowMeans(stress), missing = missing)
}

.paper_drawdown <- function(wealth) max(1 - wealth / cummax(wealth))

.paper_invest <- function(panel, horizon, common_months) {
  training <- 69L
  methods <- c("LS", "SS", "IN")
  blocks <- (length(panel$dates) - training) %/% horizon
  months <- blocks * horizon
  path <- stress <- matrix(1, months + 1L, 3L, dimnames = list(NULL, methods))
  periods <- vector("list", blocks)
  for (b in seq_len(blocks)) {
    offset <- (b - 1L) * horizon
    train <- offset + seq_len(training)
    hold <- offset + training + seq_len(horizon)
    eligible <- colSums(panel$funds[train, , drop = FALSE] == 0) == 0
    hist <- panel$funds[train, eligible, drop = FALSE]
    fits <- fit_three(sweep(hist, 1L, panel$rf[train], "-"),
                      panel$factors[train, , drop = FALSE], q = 0.05,
                      tol = 1e-8, max_iter = 1000L, capture_errors = FALSE,
                      smoothing = 0.5)
    held <- panel$funds[hold, eligible, drop = FALSE]
    counts <- missing <- integer(length(methods))
    previous <- offset + 1L
    next_rows <- previous + seq_len(horizon)
    for (j in seq_along(methods)) {
      selected <- fits[[methods[j]]]$selected
      block <- .paper_holding_wealth(held[, selected, drop = FALSE], panel$rf[hold])
      path[next_rows, j] <- path[previous, j] * block$wealth
      stress[next_rows, j] <- stress[previous, j] * block$stress
      counts[j] <- sum(selected)
      missing[j] <- block$missing
    }
    periods[[b]] <- data.frame(holding = horizon, period = b, method = methods,
                               train_start = panel$dates[train[1L]],
                               train_end = panel$dates[train[training]],
                               hold_start = panel$dates[hold[1L]],
                               hold_end = panel$dates[hold[horizon]],
                               eligible = sum(eligible), selected = counts,
                               missing = missing)
  }
  periods <- do.call(rbind, periods)
  calendar <- panel$dates[training + (0:months)]
  summary <- do.call(rbind, lapply(seq_along(methods), function(j) {
    d <- periods[periods$method == methods[j], ]
    w <- path[, j]
    common <- w[seq_len(common_months + 1L)]
    data.frame(holding = horizon, method = methods[j], periods = blocks,
               months = months, first = calendar[2L], last = calendar[months + 1L],
               min_eligible = min(d$eligible), max_eligible = max(d$eligible),
               empty = sum(d$selected == 0L), missing = sum(d$missing),
               terminal_wealth = w[months + 1L], cagr = w[months + 1L]^(12 / months) - 1,
               max_drawdown = .paper_drawdown(w),
               stress_cagr = stress[months + 1L, j]^(12 / months) - 1,
               common_last = calendar[common_months + 1L],
               common_cagr = common[common_months + 1L]^(12 / common_months) - 1,
               common_max_drawdown = .paper_drawdown(common), row.names = NULL)
  }))
  to_frame <- function(x) data.frame(holding = horizon, date = calendar, x)
  list(summary = summary, periods = periods, wealth = to_frame(path),
       stress_wealth = to_frame(stress))
}

reproduce_empirical <- function(data, rf_file = NULL, figure = NULL) {
  panel <- .paper_fund_inputs(data, rf_file)
  horizons <- c(3L, 6L, 9L)
  common_months <- min(((length(panel$dates) - 69L) %/% horizons) * horizons)
  designs <- lapply(horizons, function(h) .paper_invest(panel, h, common_months))
  result <- lapply(names(designs[[1L]]), function(key)
    do.call(rbind, lapply(designs, "[[", key)))
  names(result) <- names(designs[[1L]])
  if (!is.null(figure)) plot_empirical_wealth(result, figure)
  result
}

plot_empirical_wealth <- function(result, file) {
  paths <- result$wealth
  horizons <- c(3L, 6L, 9L)
  if (!is.data.frame(paths) ||
      !all(c("holding", "date", "LS", "SS", "IN") %in% names(paths)) ||
      !all(horizons %in% paths$holding))
    stop("Supply the result of reproduce_empirical().", call. = FALSE)
  grDevices::pdf(file, width = 480/72, height = 505/72, family = "Times",
                 title = "Rolling mutual-fund selection", useDingbats = FALSE)
  on.exit(grDevices::dev.off())
  graphics::layout(matrix(1:4), heights = c(1, 1, 1, 0.14))
  colors <- c(LS = "#202020", SS = "#2562B2", IN = "#C63032")
  line_types <- c(LS = "42", SS = "13", IN = "solid")
  ticks <- as.Date(paste0(seq(1995, 2015, 5), "-01-01"))
  yrange <- range(c(0.8, 3.5, as.matrix(paths[, names(colors)])))
  for (i in seq_along(horizons)) {
    h <- horizons[i]
    p <- paths[paths$holding == h, ]
    graphics::par(mar = c(1.9, 3.8, 1.35, 0.65), mgp = c(2.4, 0.5, 0),
                  tcl = -0.2, cex.axis = 0.78, cex.lab = 0.85)
    graphics::plot(as.numeric(p$date), p$LS, type = "n", xlab = "",
                   ylab = "Portfolio wealth", axes = FALSE, xaxs = "i", yaxs = "i",
                   xlim = as.numeric(range(paths$date)), ylim = yrange)
    ytick <- seq(1, 3.5, by = 0.5)
    graphics::abline(h = ytick, col = "#DADADA", lwd = 0.5)
    graphics::axis(1L, at = as.numeric(ticks), labels = format(ticks, "%Y"))
    graphics::axis(2L, at = ytick, labels = sprintf("%.1f", ytick), las = 1)
    for (m in names(colors))
      graphics::lines(as.numeric(p$date), p[[m]], col = colors[m],
                      lty = line_types[m], lwd = if (m == "IN") 1.5 else 1.1)
    graphics::mtext(sprintf("(%s) %d-month holding period", letters[i], h),
                    side = 3, adj = 0, line = 0.2, font = 2, cex = 0.89)
    dates_label <- if (h == 6L) "Nov 1992 - Apr 2017" else "Nov 1992 - Jul 2017"
    graphics::mtext(dates_label, side = 3, adj = 1, line = 0.2, cex = 0.77)
  }
  graphics::par(mar = rep(0, 4))
  graphics::plot.new()
  graphics::legend("center", legend = c("IN-BH", "SS-BH", "LS-BH"),
                   col = rev(colors), lty = rev(line_types), lwd = c(1.5, 1.1, 1.1),
                   horiz = TRUE, bty = "n", cex = 0.9)
  invisible(file)
}
