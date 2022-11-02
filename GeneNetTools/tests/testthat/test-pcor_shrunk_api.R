test_that("shrunk works", {
  js <- "params/shrunk.json"
  lp <- validate_json_file(js)
  validate_parameters(js,"pcor_shrunk_schema.json")
  c_pcor_shrunk(lp)
})
