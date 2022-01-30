test_that("result is not NULL", {
  test_df = data.frame(junction_beta = c("CASSS", "CARSS", "CASRS", "CASSR"),
                       v_beta = c("TRBV1", "TRBV1", "TRBV1", "TRBV1"),
                       j_beta = c("TRBJ2-4", "TRBJ2-4", "TRBJ2-4", "TRBJ2-4"),
                       id = 1:4)

  sequence_dt <- data.table::as.data.table(test_df)
  sequence_dt[, receptor_id := id]
  data.table::setkey(sequence_dt, receptor_id)
  expect_equal(is.null(calculate_scores(sequence_dt, chains="B", tmp_folder=".",
                                        scores_filename=NA, ncores=2)),
               FALSE)
})

