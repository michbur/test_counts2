library(dpcR)
library(pbapply)
sim_dat <- pblapply(1L:15*10, function(m1) {
  res <- lapply(1L:15*10, function(m2) {   
    adpcr1 <- sim_adpcr(m = m1, n = 765, times = 1e6, pos_sums = FALSE, n_panels = 3)
    adpcr2 <- sim_adpcr(m = m2, n = 765, times = 1e6, pos_sums = FALSE, n_panels = 3)
    
    adpcrs <- bind_dpcr(adpcr1, adpcr2)
    
    compb <- test_counts(adpcrs, "binomial")
    compr <- test_counts(adpcrs, "ratio")
    compp <- test_counts(adpcrs, "prop")
    
    lapply(list(glm = compb, 
                prop = compp, 
                ratio = compr), function(single_test) as.numeric(coef(single_test)[["group"]]))
  })
  names(res) <- paste0("m2.", 1L:15*10)
  res
})

names(sim_dat) <- paste0("m1.", 1L:15*10)
