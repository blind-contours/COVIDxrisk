library(here)
create_superlearner <- function(){

  ## source the custom learners built for poisson outcomes
  mean_lrnr <- Lrnr_mean$new()

  lrnr_ranger100 <- make_learner(Lrnr_ranger, num.trees = 100)
  lrnr_gam <- make_learner(Lrnr_gam)
  lrnr_polspline <- make_learner(Lrnr_polspline)


  # choose base learners
  lrnr_glm <- make_learner(Lrnr_glm)

  lrnr_ranger10 <- make_learner(Lrnr_ranger, num.trees = 10)
  lrnr_ranger50 <- make_learner(Lrnr_ranger, num.trees = 50)
  # lrnr_hal_simple <- make_learner(Lrnr_hal9001, max_degree = 2, n_folds = 2)
  lrnr_lasso <- make_learner(Lrnr_glmnet) # alpha default is 1
  lrnr_ridge <- make_learner(Lrnr_glmnet, alpha = 0)
  lrnr_elasticnet <- make_learner(Lrnr_glmnet, alpha = .5)

  # choose base learners
  lrnr_ranger100 <- make_learner(Lrnr_ranger, num.trees = 100)
  lrnr_ranger200 <- make_learner(Lrnr_ranger, num.trees = 200)
  lrnr_ranger300 <- make_learner(Lrnr_ranger, num.trees = 300)
  lrnr_ranger400 <- make_learner(Lrnr_ranger, num.trees = 400)
  lrnr_ranger500 <- make_learner(Lrnr_ranger, num.trees = 500)


  grid_params <- list(max_depth = c(2, 4, 6, 8, 10, 12),
                      eta = c(0.001, 0.01, 0.1, 0.2, 0.3),
                      nrounds = c(20, 50, 100, 200,400))

  grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)
  # params_default <- list(nthread = getOption("sl.cores.learners", 1))
  xgb_learners <- apply(grid, MARGIN = 1, function(params_tune) {
    do.call(Lrnr_xgboost$new, c(as.list(params_tune)))})

  full_lrn_earth_1 <- Lrnr_earth$new(degree = 1)
  full_lrn_earth_2 <- Lrnr_earth$new(degree = 2)
  full_lrn_earth_3 <- Lrnr_earth$new(degree = 3)
  full_lrn_earth_4 <- Lrnr_earth$new(degree = 4)

  full_lrn_poly_1 <- Lrnr_polspline$new(knots = 1)
  full_lrn_poly_2 <- Lrnr_polspline$new(knots = 2)
  full_lrn_poly_3 <- Lrnr_polspline$new(knots = 3)
  full_lrn_poly_4 <- Lrnr_polspline$new(knots = 4)

  stack <- make_learner(
    Stack,
    mean_lrnr,
    lrnr_glm,
    lrnr_ranger10,
    lrnr_ranger50,
    lrnr_ranger100,
    lrnr_ranger200,
    lrnr_ranger300,
    lrnr_ranger400,
    lrnr_ranger500,
    lrnr_lasso,
    lrnr_ridge,
    lrnr_elasticnet,
    lrnr_ranger100,
    xgb_learners[[1]],
    xgb_learners[[5]],
    xgb_learners[[10]],
    xgb_learners[[15]],
    xgb_learners[[20]],
    xgb_learners[[25]],
    xgb_learners[[31]],
    xgb_learners[[32]],
    xgb_learners[[33]],
    xgb_learners[[34]],
    xgb_learners[[40]],
    lrnr_ranger10,
    xgb_learners[[50]],
    xgb_learners[[60]],
    xgb_learners[[61]],
    xgb_learners[[70]],
    xgb_learners[[80]],
    xgb_learners[[90]],
    xgb_learners[[100]],
    xgb_learners[[110]],
    full_lrn_earth_1,
    full_lrn_earth_2,
    full_lrn_earth_3,
    full_lrn_earth_4

  )

  discrete_sl_metalrn <- Lrnr_cv_selector$new()

  discrete_sl <- Lrnr_sl$new(
    learners = stack,
    metalearner = discrete_sl_metalrn
  )

  return(discrete_sl)
}
