select_in_bh <- function(Y, factors = NULL, q = 0.1, smoothing = 0.5,
                         tol = 1e-8, max_iter = 1000L) {
  select_funds(Y, factors, method = "IN", q = q, tol = tol,
               max_iter = max_iter, smoothing = smoothing)
}

select_ss_bh <- function(Y, factors = NULL, q = 0.1,
                         tol = 1e-8, max_iter = 1000L) {
  select_funds(Y, factors, method = "SS", q = q, tol = tol,
               max_iter = max_iter)
}

select_ls_bh <- function(Y, factors = NULL, q = 0.1) {
  select_funds(Y, factors, method = "LS", q = q)
}
