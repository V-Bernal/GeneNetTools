test_that("z-score_shrunk works", {
  js <- "params/zscore.json"
  lp <- validate_json_file(js)
  validate_parameters(js,"zscore_shrunk_schema.json")
  c_zscore_shrunk(lp)
})
