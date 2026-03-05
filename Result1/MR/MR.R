library(TwoSampleMR)
library(dplyr)

all_exp <- c("finn-b-F5_ALLANXIOUS","ebi-a-GCST90018840","ukb-e-2050_AFR","ebi-a-GCST90018869",
             "ukb-d-20548_2","ebi-a-GCST90018919","finn-b-F5_BIPO")
out <- c("ebi-a-GCST90027158",
         "ukb-d-HEARTFAIL",
         "finn-b-I9_HYPTENS",
         "ieu-a-1108",
         "ieu-b-4966",
         "ebi-a-GCST90018894",
         "ebi-a-GCST90014023"
)

result_all <- list()
hetero_all <- list()
loo_all <- list()

k <- 1
for (i in seq_along(all_exp)) {
  exp <- all_exp[i]

  exposure_dat3 <- extract_instruments(exp, p1 = 1e-4, clump = TRUE, r2 = 0.0001, kb = 100000)


  if (is.null(exposure_dat3) || nrow(exposure_dat3) == 0) next

  for (j in seq_along(out)) {
    outcomes <- out[j]

    outcome_dat <- extract_outcome_data(snps = exposure_dat3$SNP, outcomes = outcomes)
    if (is.null(outcome_dat) || nrow(outcome_dat) == 0) next

    dat <- harmonise_data(exposure_dat = exposure_dat3, outcome_dat = outcome_dat)
    if (is.null(dat) || nrow(dat) == 0) next


    res <- mr(dat)
    result_all[[k]] <- res


    het <- mr_heterogeneity(dat)
    hetero_all[[k]] <- het

    # Leave-one-out

    if (length(unique(dat$SNP)) >= 3) {
      loo <- mr_leaveoneout(dat)
      loo_all[[k]] <- loo
    }

    k <- k + 1
  }
}


result_all_df <- bind_rows(result_all)
hetero_all_df <- bind_rows(hetero_all)
loo_all_df <- bind_rows(loo_all)


result_ivw <- result_all_df %>% filter(method == "Inverse variance weighted")

write.csv(result_ivw, "mr_ivw_results.csv", row.names = FALSE)
write.csv(hetero_all_df, "mr_heterogeneity.csv", row.names = FALSE)
write.csv(loo_all_df, "mr_leaveoneout.csv", row.names = FALSE)
