ieugwasr::get_opengwas_jwt()
Sys.getenv("R_ENVIRON_USER")

library(TwoSampleMR)

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
result_all <- as.data.frame(matrix(ncol=9,nrow=0))
colnames(result_all) <- colnames(res)
for (i in 1:length(all_exp)){
  exp <- all_exp[i]
  exposure_dat3<- extract_instruments(exp,p1=1e-4,clump =T,r2=0.0001,kb=100000)
  for (j in 1:length(out)){
    outcomes <- out[j]
    outcome_dat <-extract_outcome_data(snps=exposure_dat3$SNP, outcomes= outcomes)
    dat <- harmonise_data(
      exposure_dat = exposure_dat3,
      outcome_dat = outcome_dat)
    res <- mr(dat)
    
    result_all <- rbind(result_all,res)
  }
}

data <- all_results[[1]]
for ( i in 2:length(all_results)){
  data <- rbind(data,all_results[[i]])
}

data2 <-  result_all %>%
  dplyr::filter(method == "Inverse variance weighted")
write.csv(data2,"data3.csv")
############
data <- read.csv("data(1).csv")
library(ggplot2)
data2 <- data[,c(1,3,6,8)]
colnames(data2) <- c("x","y","b","p")
data <- data2
##########
library(ggplot2)
library(reshape2)
heatmap_data <- dcast(data2, x ~ y, value.var = "p")
rownames(heatmap_data) <- heatmap_data$x
heatmap_data$x <- NULL
heatmap_matrix <- as.matrix(heatmap_data)
ggplot(data, aes(x = x, y = y)) +
  geom_point(aes(size = 10^abs(b), color = -log10(p)), alpha = 0.7) +
  scale_color_gradientn(values = scales::rescale(c(0,0.01,0.03,0.05,0.1,0.15,0.2,0.25,0.4,1)),
                        colors = c("#39489f", "#39bbec", "#f9ed36", "#f38466", "#b81f25")) +
  scale_size_continuous(range = c(5, 10)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none") +  
  guides(color = "none", size = "none")  



