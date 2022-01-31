test_that("input check works", {
  expect_error(add_cdr3_anchors(example_TCR_df, chains=T, species = "human"))
  expect_error(add_cdr3_anchors(example_TCR_df, chains="B", species = "human"))
  expect_error(add_cdr3_anchors(example_TCR_df, chains="B", species = "homo"))

  test_df <- example_TCR_df
  test_df$cdr3_beta <- substr(test_df$junction_beta, 2, nchar(test_df$junction_beta) - 2)
  expect_error(add_cdr3_anchors(test_df, chains="A", species = "human"))
  expect_error(add_cdr3_anchors(test_df, chains="B", species = "human"), NA)
})

test_that("gene filtering works", {
  test_df <- data.frame(cdr3_beta = c("SMTH", "SMTHS"),
                        j_beta = c("TRBJ2-4", "TRBJ"))
  res = add_cdr3_anchors(test_df, chains="B", species = "human")
  expect_equal(nrow(res), 2)
  expect_equal(sum(is.na(res$junction_beta)), 1)
})
