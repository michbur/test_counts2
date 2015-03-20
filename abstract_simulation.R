library(dpcR)
library(pbapply)
library(reshape2)

#copmarision of adpcr seting
adpcr_comp <- pbsapply(1L:300, function(dummy) {
  sim_dat2 <- lapply(1L:10*10, function(m1) {
    m2_vector <- m1 + 1L:5*10
    res <- lapply(m2_vector, function(m2) {   
      adpcr1 <- sim_adpcr(m = m1, n = 765, times = 1e6, pos_sums = FALSE, n_panels = 3)
      adpcr2 <- sim_adpcr(m = m2, n = 765, times = 1e6, pos_sums = FALSE, n_panels = 3)
      
      adpcrs <- bind_dpcr(adpcr1, adpcr2)
      
      compb <- test_counts(adpcrs, "binomial")
      compr <- test_counts(adpcrs, "ratio")
      
      lapply(list(glm = compb, 
                  ratio = compr), function(single_test) as.numeric(coef(single_test)[["group"]]))
    })
    names(res) <- paste0("m2.", m2_vector)
    res
  })
  
  names(sim_dat2) <- paste0("m1.", 1L:10*10)
  c(mean((3 - sapply(sim_dat2, function(i)
    sapply(i, function(j) {
      sum(j[["glm"]][1L:3] != j[["glm"]][4L:6])
    })))^2), mean((3 - sapply(sim_dat2, function(i)
      sapply(i, function(j) {
        sum(j[["ratio"]][1L:3] != j[["ratio"]][4L:6])
      })))^2))
})


rownames(adpcr_comp) <- c("GLM", "MT")
madpcr_comp <- melt(sqrt(adpcr_comp)/6)
colnames(madpcr_comp) <- c("method", "repetition", "value")
save(madpcr_comp, file = "adpcr_comp.RData")
