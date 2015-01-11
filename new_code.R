library(dpcR)

adpcr1 <- sim_adpcr(m = 200, n = 765, times = 1e5, pos_sums = FALSE, n_panels = 2)
adpcr2 <- sim_adpcr(m = 400, n = 765, times = 1e5, pos_sums = FALSE, n_panels = 2)
adpcr3 <- sim_adpcr(m = 600, n = 765, times = 1e5, pos_sums = FALSE, n_panels = 2)

adpcrs <- bind_dpcr(adpcr1, adpcr2, adpcr3)

badpcrs <- binarize(adpcrs)



#p-value for confidence levels
p_val <- p.adjust(rep(0.05, ncol(badpcrs)), "BH") 
conf_ints <- dpcR:::fl(binom.confint(colSums(badpcrs), 
                                     slot(badpcrs, "n"), 
                                     conf.level = 1 - p_val[1],
                                     "wilson")[, 4L:6])