library(dpcR)
library(binom)

sapply(0L:20*20, function(num_mol) {
  number_of_exps <- 100
  dat <- sim_ddpcr_bkm(num_mol, n_exp = number_of_exps, type = "tnp")
  
  sum <- summary(dat, print = FALSE)[["summary"]]
  
  
  positives <- slot(dat, ".Data")
  total <- slot(dat, "n")
  
  sidak_ci <- dpcR:::fl(binom.confint(as.vector(positives), 
                                      total, 
                                      conf.level = (1 - 0.05)^(1/length(dat)),
                                      "wilson")[, 4L:6])[, 2L:3]
  
  nadj_ci <- dpcR:::fl(binom.confint(as.vector(positives), 
                                     total, 
                                     conf.level = 0.95,
                                     "wilson")[, 4L:6])[, 2L:3]
  
  dube_ci <- sum[sum[["method"]] == "dube", c("lambda.low", "lambda.up")]
  bhat_ci <- sum[sum[["method"]] == "bhat", c("lambda.low", "lambda.up")]
  
  c(sidak = sum(apply(cbind(rep(num_mol, number_of_exps), sidak_ci), 1, 
                      function(i)
                        abs(i[1]) > abs(i[2]) && i[1] < i[3]))/number_of_exps,
    nadj = sum(apply(cbind(rep(num_mol, number_of_exps), nadj_ci), 1, 
                     function(i)
                       abs(i[1]) > abs(i[2]) && i[1] < i[3]))/number_of_exps,
    dube = sum(apply(cbind(rep(num_mol, number_of_exps), dube_ci), 1, 
                     function(i)
                       abs(i[1]) > abs(i[2]) && i[1] < i[3]))/number_of_exps,
    bhat = sum(apply(cbind(rep(num_mol, number_of_exps), bhat_ci), 1, 
                     function(i)
                       abs(i[1]) > abs(i[2]) && i[1] < i[3]))/number_of_exps)
})
