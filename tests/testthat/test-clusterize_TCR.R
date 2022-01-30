Sys.setenv("R_TESTS" = "")

test_that("input check works", {
  expect_error(clusterize_TCR(example_TCR_df, chains="", tmp_folder=".", id_col="id", ncores=2))

  test_df = example_TCR_df %>% dplyr::rename(cdr3_beta = junction_beta)
  expect_error(clusterize_TCR(test_df, chains="AB", tmp_folder=".", id_col="id", ncores=2))

  test_df = example_TCR_df %>% dplyr::select(-contains("alpha"))
  expect_error(clusterize_TCR(test_df, chains="AB", tmp_folder=".", id_col="id", ncores=2))

  test_df = example_TCR_df %>% dplyr::select(-contains("beta"))
  expect_error(clusterize_TCR(test_df, chains="B", tmp_folder=".", id_col="id", ncores=2))
})

test_that("edge case works", {
  test_df = example_TCR_df[0,]
  expect_error(clusterize_TCR(test_df, chains="AB", tmp_folder=".", id_col="id", ncores=2))

  test_df = example_TCR_df[1,]
  expect_error(clusterize_TCR(test_df, chains="AB", tmp_folder=".", id_col="id", ncores=2))
})

test_that("filters work", {
  test_df = data.frame(junction_beta = c("C", "CASS*", "CASSS", "CASSS", "CASSS", "CASSS"),
                       v_beta = c("TRBV1", "TRBV1", "TRV", "TRBV1", "TRBV1", "TRBV1"),
                       j_beta = c("TRBJ2-4", "TRBJ2-4", "TRBJ2-4", "TRBJ", "TRBJ2-4", "TRBJ2-4"),
                       id = 1:6)
  out_len = nrow(clusterize_TCR(test_df, chains="B", tmp_folder=".", id_col="id", ncores=2))
  expect_equal(out_len, 2)
})

test_that("no error occurs", {
  expect_error(clusterize_TCR(example_TCR_df, chains="AB", tmp_folder=".", id_col="id", ncores=2), NA)
})
