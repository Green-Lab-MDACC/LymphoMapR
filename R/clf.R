# naive bayes
nb_clf <- function(mt, y) {
  stopifnot(all.equal(rownames(mt), names(y)))
  stopifnot(is.factor(y))
  nb <- naive_bayes(mt, y)
  return(nb)
}


# naive bayes classifer evaluation via cross validation
nb_clf_cv <- function(mt, y, cv_folds = 10, iseed = 10101) {
  stopifnot(all.equal(rownames(mt), names(y)))
  stopifnot(is.factor(y))
  
  # naive bayes: cross validation performance
  message(paste0("Evaluating naive Bayes model performance using ", cv_folds, "-fold cross validation"))
  set.seed(iseed)
  shuffled_indices <- sample(1:nrow(mt))
  stride <- ceiling(nrow(mt)/cv_folds)
  s <- 1
  tb <- tibble()
  while (s <= nrow(mt)) {
    e <- min(s + stride - 1, nrow(mt))
    ival <- shuffled_indices[c(s:e)]
    itrain <- shuffled_indices[-c(s:e)]
    
    fit = naive_bayes(mt[itrain, ], y[itrain])
    pred.res = predict(fit, mt[ival, ], type = "class")
    tb <- tb %>%
      bind_rows(tibble(sample_name = rownames(mt[ival, ]),
                       actual = y[ival],
                       predicted = as.vector(pred.res)))
    
    s <- e + 1
  }
  
  archetypes <- levels(y)
  tb <- tb %>%
    mutate(actual = factor(actual,
                           levels = ),
           predicted = factor(predicted,
                              levels = archetypes))
  
  # confusion mat
  confmat <- caret::confusionMatrix(tb$predicted, tb$actual)
  print(confmat)
  return(confmat)
}


# LymphoMAP prediction using a trained naive bayes model
nb_clf_predict <- function(nb, mt) {
  stopifnot(all.equal(colnames(mt), colnames(nb$data$x)))
  
  # pred labels
  pred.res <- predict(nb, mt, type = "class")
  tb <- tibble(sample_id = rownames(mt),
               LymphoMAP = as.vector(pred.res))
  tb_pred <- tb
  
  # pred probabilities
  pred.prob <- predict(nb, mt, type = "prob")
  rownames(pred.prob) <- rownames(mt)
  df <- pred.prob %>% as.data.frame()
  tb <- tibble(sample_id = rownames(df)) %>%
    bind_cols(as_tibble(df))
  tb_prob <- tb
  
  # get prob for assigned archetype
  tb <- tb_prob %>%
    gather("LymphoMAP", "LymphoMAP_prob", -sample_id)
  tb_pred <- tb_pred %>%
    left_join(tb, by = c("sample_id", "LymphoMAP")) %>%
    left_join(tb_prob, by = "sample_id") %>%
    mutate(LymphoMAP = factor(LymphoMAP, levels = levels(nb$data$y)))
  return(tb_pred)
}


