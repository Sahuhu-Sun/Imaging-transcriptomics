library(readxl)
library(EWCE) 
library(ggplot2)
library(gridExtra)
library(dplyr)
library(purrr)

#### cell specificity  #### 
data0 <- read_excel("PLS_res.xlsx", sheet = "Sheet1")

percen<-0.3 # top30% genes
Level<-1 # level 1 Cell annotation
rep<-10000 # 10,000 repetitions

pos <- data0 %>% filter(Weights > 0)
neg <- data0 %>% filter(Weights < 0)
plots <- list()
datannames <- list()
for (i in 1:2) {
  
  if (i == 1) {
    re <- 'Up'
    data <- pos
  } else {
    re <- 'Down'
    data <- neg
  }
  
  
  for (j in 1:3) {
    # j is cell type dateset 
    # 'AHBA' 1 # DRONC 2 # Karolinska Institute
    print(i,j)
    
    if (j == 1) {
      dataname <- 'AHBA'
      ctd_path<-"ctd_AHBA_human.rds"
      ctd <- readRDS(ctd_path)
      sct<-'human'
    } else if (j == 2) {
      dataname <- 'DroNC-seq'
      ctd_path<-"ctd_DRONC_human.rds"
      ctd <- readRDS(ctd_path)
      sct<-'human'
    } else {
      dataname <- 'Karolinska'
      ctd <- ewceData::ctd()
      sct<-'mouse'
    }
    
    # Select the top percent of genes
    n <- nrow(data)
    top10_n <- ceiling(n * percen)
    top10_data <- data %>%
      mutate(abs_w = abs(Weights)) %>%  
      arrange(desc(abs_w)) %>%     
      slice(1:top10_n) %>% # 
      select(Genes)                                      
    
    hits <- top10_data$Genes
    full_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd, hits = hits,
                                                    sctSpecies = sct, genelistSpecies = "human",
                                                    reps = rep, annotLevel = Level )
    p <- EWCE::ewce_plot(total_res = full_results$results,mtc_method = "BH")
    plots <- c(plots, list(p$plain))
    datannames <- c(datannames, paste0(dataname,'-',re))
  }
}

# get plots values
library(openxlsx)
wb <- createWorkbook()

for (i in seq_along(plots)) {
  plot_data <- plots[[i]]$data
  sheet_name <- datannames[[i]]
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, x = plot_data)
}
saveWorkbook(wb, file = "plots_data_celltype.xlsx", overwrite = TRUE)
arranged_plot <- grid.arrange(grobs = plots, ncol = 3)




#### tissue specificity  #### 
ctd_path<-"ctd_Region.rda"
ctd <- EWCE::load_rdata(ctd_path)
data0 <- read_excel("PLS_res.xlsx", sheet = "Sheet1")

percen<-0.3 # top30% genes
Level<-1 # level 1 Cell annotation
rep<-10000 # 10,000 repetitions

pos <- data0 %>% filter(Weights > 0)
neg <- data0 %>% filter(Weights < 0)
datannames <- list('Up','Down')
plots <- list()
for (i in 1:2) {
  if (i == 1) {
    re <- 'Up'
    data <- pos
  } else {
    re <- 'Down'
    data <- neg
  }
  
  sct<-'human'
  n <- nrow(data)
  top10_n <- ceiling(n * percen)
  top10_data <- data %>%
    mutate(abs_w = abs(Weights)) %>%  
    arrange(desc(abs_w)) %>%     
    slice(1:top10_n) %>% 
    select(Genes)                                   
  
  hits <- top10_data$Genes
  full_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd, hits = hits,
                                                  sctSpecies = sct, genelistSpecies = "human",
                                                  reps = rep, annotLevel = Level )
  p <- EWCE::ewce_plot(total_res = full_results$results,mtc_method = "BH")
  plots <- c(plots, list(p$plain))
  
}
library(openxlsx)
wb <- createWorkbook()
for (i in seq_along(plots)) {
  plot_data <- plots[[i]]$data
  sheet_name <- datannames[[i]]
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, x = plot_data)
}
saveWorkbook(wb, file = "plots_data_region.xlsx", overwrite = TRUE)



#### developmental stage specificity  #### 
ctd_path<-"ctd_Lifespan.rda"
ctd <- EWCE::load_rdata(ctd_path)
data0 <- read_excel("PLS_res.xlsx", sheet = "Sheet1")

percen<-0.3 # top30% genes
Level<-1 # level 1 Cell annotation
rep<-10000 # 10,000 repetitions

pos <- data0 %>% filter(Weights > 0)
neg <- data0 %>% filter(Weights < 0)
datannames <- list('Up','Down')
plots <- list()
for (i in 1:2) {
  if (i == 1) {
    re <- 'Up'
    data <- pos
  } else {
    re <- 'Down'
    data <- neg
  }
  
  sct<-'human'
  n <- nrow(data)
  top10_n <- ceiling(n * percen)
  top10_data <- data %>%
    mutate(abs_w = abs(Weights)) %>%  
    arrange(desc(abs_w)) %>%     
    slice(1:top10_n) %>% 
    select(Genes)                                   
  
  hits <- top10_data$Genes
  full_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd, hits = hits,
                                                  sctSpecies = sct, genelistSpecies = "human",
                                                  reps = rep, annotLevel = Level )
  p <- EWCE::ewce_plot(total_res = full_results$results,mtc_method = "BH")
  plots <- c(plots, list(p$plain))
  
}
library(openxlsx)
wb <- createWorkbook()
for (i in seq_along(plots)) {
  plot_data <- plots[[i]]$data
  sheet_name <- datannames[[i]]
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, x = plot_data)
}
saveWorkbook(wb, file = "plots_data_stage.xlsx", overwrite = TRUE)


















