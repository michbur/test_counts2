library(dpcR)
library(binom)

sapply(0L:20*20, function(num_mol) {
  number_of_exps <- 1000
  adpcr1 <- sim_adpcr(m = num_mol, n = 1000, times = 100000, pos_sums = FALSE, 
                      n_panels = number_of_exps)
  dat <- binarize(adpcr1)
  
  sum <- summary(dat, print = FALSE)[["summary"]]
  
  
  positives <- colSums(dat > 0, na.rm = TRUE)
  total <- slot(dat, "n")
  
  sidak_ci <- dpcR:::fl(binom.confint(positives, 
                                      total, 
                                      conf.level = (1 - 0.05)^(1/ncol(dat)),
                                      "wilson")[, 4L:6])[, 2L:3]
  
  nadj_ci <- dpcR:::fl(binom.confint(positives, 
                                     total, 
                                     conf.level = 0.95,
                                     "wilson")[, 4L:6])[, 2L:3]
  
  dube_ci <- sum[sum[["method"]] == "dube", c("lambda.low", "lambda.up")]
  bhat_ci <- sum[sum[["method"]] == "bhat", c("lambda.low", "lambda.up")]
  
  c(sidak = sum(apply(cbind(rep(num_mol, number_of_exps)/total, sidak_ci), 1, 
                      function(i)
                        abs(i[1]) > abs(i[2]) && i[1] < i[3]))/number_of_exps,
    nadj = sum(apply(cbind(rep(num_mol, number_of_exps)/total, nadj_ci), 1, 
                     function(i)
                       abs(i[1]) > abs(i[2]) && i[1] < i[3]))/number_of_exps,
    dube = sum(apply(cbind(rep(num_mol, number_of_exps)/total, dube_ci), 1, 
                     function(i)
                       abs(i[1]) > abs(i[2]) && i[1] < i[3]))/number_of_exps,
    bhat = sum(apply(cbind(rep(num_mol, number_of_exps)/total, bhat_ci), 1, 
                     function(i)
                       abs(i[1]) > abs(i[2]) && i[1] < i[3]))/number_of_exps)
})
