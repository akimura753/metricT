ext <- function(f) system.file("extdata", f, package = "metricT")
json_num <- function(txt, key) {
  m <- regmatches(txt, regexpr(sprintf("\"%s\": *\\[[^]]*\\]", key), txt))
  suppressWarnings(as.numeric(strsplit(gsub(sprintf("\"%s\": *\\[|\\]|\\s", key), "", m), ",")[[1]]))
}

test_that("connectivity matches the Python reference", {
  ep <- mt_read_npy(ext("tiny_epochs.npy"))
  expect_equal(dim(ep), c(6L, 4L, 128L))
  con <- mt_connectivity(ep, 64, list(theta = c(4, 8), alpha = c(8, 13), beta = c(13, 30)))
  txt <- paste(readLines(ext("tiny_conn_ref.json"), warn = FALSE), collapse = "")
  for (b in c("theta", "alpha", "beta")) {
    blk <- regmatches(txt, regexpr(sprintf("\"%s\": *\\{[^}]*\\}", b), txt))
    expect_equal(unname(con[[b]]$wpli), json_num(blk, "wpli"), tolerance = 1e-12)
    expect_equal(unname(con[[b]]$msc), json_num(blk, "msc"), tolerance = 1e-12)
  }
  expect_true(all(unlist(con) >= 0 & unlist(con) <= 1 + 1e-12))
})

test_that("Metric-T replays the Python permutations exactly", {
  z <- mt_read_npz(ext("B_features.npz"))
  expect_equal(as.vector(table(z$group)[c("AD", "FTD", "CN")]), c(9L, 8L, 7L))
  pidx <- mt_read_npy(ext("B_perm_index.npy")) + 1L
  res <- metric_t(z$features, z$group, comparisons = list(c("AD", "CN"), c("FTD", "CN"), c("AD", "FTD")),
                  control = "CN", perm_index = pidx)
  txt <- paste(readLines(ext("B_ref.json"), warn = FALSE), collapse = "")
  t <- res$table
  expect_equal(t$DC_wPLI, json_num(txt, "DC_wPLI"), tolerance = 1e-10)
  expect_equal(t$DC_MSC, json_num(txt, "DC_MSC"), tolerance = 1e-10)
  expect_equal(t$T, json_num(txt, "T_obs"), tolerance = 1e-10)
  expect_equal(t$p_raw, json_num(txt, "p_raw"), tolerance = 1e-12)
  expect_equal(t$p_holm, json_num(txt, "p_holm"), tolerance = 1e-12)
  expect_equal(t$q_BH, json_num(txt, "q_BH"), tolerance = 1e-12)
  expect_equal(t$p_WY, json_num(txt, "p_WY"), tolerance = 1e-12)
})

test_that("DC, reversal flag and CSV round trip behave", {
  g1 <- matrix(c(1, 1, 0, 0), 2); g2 <- matrix(c(0, 0, 1, 1), 2)
  expect_equal(mt_dc(g1, g2), 50)
  sim <- mt_simulate(n = c(A = 8, B = 8), n_ch = 4, n_epochs = 6, seed = 2)
  r1 <- metric_t(sim$features, sim$group, n_perm = 199, seed = 1)
  r2 <- metric_t(sim$features, sim$group, n_perm = 199, seed = 1)
  expect_identical(r1$table, r2$table)
  expect_true(all(r1$table$p_raw > 0 & r1$table$p_raw <= 1))
  expect_equal(r1$table$reversal, (r1$table$DC_wPLI - 50) * (r1$table$DC_MSC - 50) < 0)
  f <- tempfile(fileext = ".csv")
  mt_write_features_csv(sim$features, sim$group, f)
  back <- mt_read_features_csv(f)
  r3 <- metric_t(back$features, back$group, n_perm = 199, seed = 1)
  expect_equal(r3$table$T, r1$table$T)
  # swapping the groups mirrors DC around 50 when there are no exact ties in means
  r4 <- metric_t(sim$features, sim$group, comparisons = list(c("B", "A")), n_perm = 99, seed = 1)
  expect_equal(r4$table$DC_wPLI, 100 - r1$table$DC_wPLI)
  expect_equal(r4$table$T, -r1$table$T)
})

test_that("epoching follows length and overlap", {
  x <- matrix(seq_len(2000), ncol = 2)
  e <- mt_epoch(x, sfreq = 100, epoch_len = 2, overlap = 0.5)
  expect_equal(dim(e), c(9L, 2L, 200L))
  expect_equal(e[2, 1, 1], 101)
})
