#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%% 
#%%%%%%%%%%%%%%% FastGxC Manuscript Figures
#%%%%%%%%%%%%%%% Lena Krockenberger, Andrew Lu and Brunilda Balliu
#%%%%%%%%%%%%%%% Thursday, January 14 2021
#%%%%%%%%%%%%%%% 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rm(list=ls())

# imports
library(dplyr)
library(ggplot2)
library(ggpattern)
library(tidyverse)
library(cowplot)
library(RColorBrewer)
library(scales)
library(ggrepel)
library(mvtnorm)
library(reshape2)
library(pheatmap)
library(ggcorrplot)
library(colorspace)
library(formattable)
library(htmltools)
library(webshot) #webshot::install_phantomjs()   
library(magrittr)
library(ggpubr)
library(eulerr)
library(ggbreak)
library(patchwork)
library(stats)
library(data.table)
library(openxlsx)
library(scales)
library(mppa)
library(WriteXLS)

setwd("./FastGxC/manuscript/scripts/")
setwd("../Figures/")
source(file = '../scripts/00_functions.R')


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%% GTEx and Methods Colors and Abbreviations
# NOTE: ALWAYS RUN THIS
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(1){
  ################################ set up ################################
  gtex_colors.tib <- read_csv("../Input_Files/FigureS1_abbreviations_and_colors/GTEx_v8.AllTissues.Colors_Abbreviations.csv", show_col_types = FALSE) %>% dplyr::select(tissue, color_hex)
  gtex_colors.vec <- gtex_colors.tib$color_hex
  names(gtex_colors.vec) <- gtex_colors.tib$tissue
  
  sc_colors.tib <- read_csv("../Input_Files/FigureS1_abbreviations_and_colors/scMeta.AllContexts.Colors_Abbreviations.csv", show_col_types = FALSE) %>% dplyr::select(tissue, color_hex)
  sc_colors.vec <- sc_colors.tib$color_hex
  names(sc_colors.vec) <- sc_colors.tib$tissue
  
  gtex_abbrev.tib <- read_csv("../Input_Files/FigureS1_abbreviations_and_colors/GTEx_v8.AllTissues.Colors_Abbreviations.csv", show_col_types = FALSE) %>% dplyr::select(tissue, abbreviation)
  gtex_abbrev.vec <- gtex_abbrev.tib$abbreviation
  names(gtex_abbrev.vec) <- gtex_abbrev.tib$tissue
  
  sc_abbrev.tib <- read_csv("../Input_Files/FigureS1_abbreviations_and_colors/scMeta.AllContexts.Colors_Abbreviations.csv", show_col_types = FALSE) %>% dplyr::select(tissue, abbreviation)
  sc_abbrev.vec <- sc_abbrev.tib$abbreviation
  names(sc_abbrev.vec) <- sc_abbrev.tib$tissue
  
  ############################# use these to add color or abbreviations ################################
  gtex_colors_scale <- scale_colour_manual(values = gtex_colors.vec, guide = "none")
  gtex_fill_scale <- scale_fill_manual(values = gtex_colors.vec, guide = "none")
  gtex_abb_scale <- scale_x_discrete(labels = gtex_abbrev.vec)
  
  sc_colors_scale <- scale_colour_manual(values = sc_colors.vec, guide = "none")
  sc_fill_scale <- scale_fill_manual(values = sc_colors.vec, guide = "none")
  sc_abb_scale <- scale_x_discrete(labels = sc_abbrev.vec)
  
  both_colors_scale <- scale_colour_manual(values = c(gtex_colors.vec, sc_colors.vec), guide = "none")
  both_fill_scale <- scale_fill_manual(values = c(gtex_colors.vec, sc_colors.vec), guide = "none")
  both_abb_scale <- scale_x_discrete(labels = c(gtex_abbrev.vec, sc_abbrev.vec))
  
  ########### standardize colors for methods: CxC, FastGxC, LM-GxC, LMM-GxC #############################
  manuscript_colors_vec = c(
    "CxC" = "#c87e7e",
    "FastGxC"="#56A3E9",
    "LM-GxC" = "#000000",
    "LMM-GxC" = "#009E73",
    "Metasoft_CxC_FE"= "#c87e7e",
    "Metasoft_CxC_RE"= "#c87e7e",
    "Metasoft_CxC_Mvalues_sigRE2" = "#c87e7e",
    "Metasoft_CxC_HetQ" = "#c87e7e",
    "Metasoft_CxC_Mvalues"= "#c87e7e",
    "Metasoft_CxC_HetQ_sigRE2" = "#c87e7e",
    "Metasoft_CxC_RE2" = "#c87e7e",
    "Metasoft_FastGxC_FE"= "#56A3E9",
    "Metasoft_FastGxC_RE"= "#56A3E9",
    "Metasoft_FastGxC_Mvalues_sigRE2" = "#56A3E9",
    "Metasoft_FastGxC_HetQ" = "#56A3E9",
    "Metasoft_FastGxC_Mvalues"= "#56A3E9",
    "Metasoft_FastGxC_HetQ_sigRE2" = "#56A3E9",
    "Metasoft_FastGxC_RE2" = "#56A3E9",
    "MTFERE2_FE"= "#B6A3FF",
    "MTFERE2_RE"= "#B6A3FF",
    "MTFERE2_Mvalues_sigRE2" = "#B6A3FF",
    "MTFERE2_HetQ" = "#B6A3FF",
    "MTFERE2_Mvalues"= "#B6A3FF",
    "MTFERE2_HetQ_sigRE2" = "#B6A3FF",
    "MTFERE2_RE2" = "#B6A3FF")
  manuscript_colors_vec_noLMGxC = c(
    "CxC" = "#c87e7e",
    "FastGxC"="#56A3E9",
    "LMM-GxC" = "#009E73",
    "Metasoft_CxC_FE"= "#c87e7e",
    "Metasoft_CxC_RE"= "#c87e7e",
    "Metasoft_CxC_Mvalues_sigRE2" = "#c87e7e",
    "Metasoft_CxC_HetQ" = "#c87e7e",
    "Metasoft_CxC_Mvalues"= "#c87e7e",
    "Metasoft_CxC_HetQ_sigRE2" = "#c87e7e",
    "Metasoft_CxC_RE2" = "#c87e7e",
    "Metasoft_FastGxC_FE"= "#56A3E9",
    "Metasoft_FastGxC_RE"= "#56A3E9",
    "Metasoft_FastGxC_Mvalues_sigRE2" = "#56A3E9",
    "Metasoft_FastGxC_HetQ" = "#56A3E9",
    "Metasoft_FastGxC_Mvalues"= "#56A3E9",
    "Metasoft_FastGxC_HetQ_sigRE2" = "#56A3E9",
    "Metasoft_FastGxC_RE2" = "#56A3E9",
    "MTFERE2_FE"= "#B6A3FF",
    "MTFERE2_RE"= "#B6A3FF",
    "MTFERE2_Mvalues_sigRE2" = "#B6A3FF",
    "MTFERE2_HetQ" = "#B6A3FF",
    "MTFERE2_Mvalues"= "#B6A3FF",
    "MTFERE2_HetQ_sigRE2" = "#B6A3FF",
    "MTFERE2_RE2" = "#B6A3FF")
  
  
  manuscript_shape_vec = c(
    "CxC" = 16,
    "FastGxC"= 16,
    "LM-GxC" = 16,
    "LMM-GxC" = 16,
    "Metasoft_CxC_FE"= 17,
    "Metasoft_CxC_RE"= 3,
    "Metasoft_CxC_Mvalues_sigRE2" = 4,
    "Metasoft_CxC_HetQ" = 5,
    "Metasoft_CxC_Mvalues"= 6,
    "Metasoft_CxC_HetQ_sigRE2" = 7,
    "Metasoft_CxC_RE2" = 4,
    "Metasoft_FastGxC_FE"= 17,
    "Metasoft_FastGxC_RE"= 3,
    "Metasoft_FastGxC_Mvalues_sigRE2" = 4,
    "Metasoft_FastGxC_HetQ" = 5,
    "Metasoft_FastGxC_Mvalues"= 6,
    "Metasoft_FastGxC_HetQ_sigRE2" = 7,
    "Metasoft_FastGxC_RE2" = 4,
    "MTFERE2_FE"= 17,
    "MTFERE2_RE"= 3,
    "MTFERE2_Mvalues_sigRE2" = 4,
    "MTFERE2_HetQ" = 5,
    "MTFERE2_Mvalues"= 6,
    "MTFERE2_HetQ_sigRE2" = 7,
    "MTFERE2_RE2" = 4)
  manuscript_shape_vec_noLMGxC = c(
    "CxC" = 16,
    "FastGxC"=16,
    "LMM-GxC" = 16,
    "Metasoft_CxC_FE"= 17,
    "Metasoft_CxC_RE"= 3,
    "Metasoft_CxC_Mvalues_sigRE2" = 4,
    "Metasoft_CxC_HetQ" = 5,
    "Metasoft_CxC_Mvalues"= 6,
    "Metasoft_CxC_HetQ_sigRE2" = 7,
    "Metasoft_CxC_RE2" = 4,
    "Metasoft_FastGxC_FE"= 17,
    "Metasoft_FastGxC_RE"= 3,
    "Metasoft_FastGxC_Mvalues_sigRE2" = 4,
    "Metasoft_FastGxC_HetQ" = 5,
    "Metasoft_FastGxC_Mvalues"= 6,
    "Metasoft_FastGxC_HetQ_sigRE2" = 7,
    "Metasoft_FastGxC_RE2" = 4,
    "MTFERE2_FE"= 17,
    "MTFERE2_RE"= 3,
    "MTFERE2_Mvalues_sigRE2" = 4,
    "MTFERE2_HetQ" = 5,
    "MTFERE2_Mvalues"= 6,
    "MTFERE2_HetQ_sigRE2" = 7,
    "MTFERE2_RE2" = 4)
  
  methods_color_scale <- scale_colour_manual(values = manuscript_colors_vec)
  methods_fill_scale <- scale_fill_manual(values = manuscript_colors_vec)
  
  methods_color_scale_noLMGxC = scale_colour_manual(values = manuscript_colors_vec_noLMGxC)
  
  methods_shape_scale <- scale_shape_manual(values = manuscript_shape_vec)
  methods_shape_scale_noLMGxC = scale_shape_manual(values = manuscript_shape_vec)
  
  blank_theme = theme_bw() +
    theme(
      # panel.grid.major = element_blank(),
      # panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black")
    )
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%% Figure 1 : Illustration of FastGxC method
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(1){
  
  N=500
  maf=.3
  set.seed(09301987)
  genos = rbinom(N, 2, maf)
  zval = 2.5 # 1.96
  
  # assume all contexts have correlation of w_corr within people
  # w_corr=0
  w_corr=.3
  nT=5
  sigma = matrix(w_corr,nrow=nT,ncol=nT)
  diag(sigma) = 2
  
  Betas=data.frame(rbind(rep(0, nT),
                         rep(0.5, nT),
                         c(rep(0,nT-1),.5),
                         c(rep(.5, nT-1), 1),
                         c(-2,-1,0,1,2),
                         seq(0,1,length.out=nT)), scenario= c("No heterogeneity\n(No eQTL)", 
                                                              "No heterogeneity", 
                                                              "Single-context\nheterogeneity\nnoShared",
                                                              "Single-context\nheterogeneity",
                                                              "Extensive\nheterogeneity\nnoShared",
                                                              "Extensive\nheterogeneity"))
  
  my_plot_data_all_scenarios= NULL
  for(i in 1:nrow(Betas)){
    betas=unlist(Betas[i,1:nT])
    mus=rep(0, nT)
    
    Y = matrix(0,nrow=N,ncol=nT, dimnames = list(1:N, paste0("T",1:nT))) 
    for (j in 1:nT)  Y[,j] = mus[j] + genos*betas[j]
    
    # get noise per individual
    for(j in 1:N) Y[j,] = Y[j,] + rmvnorm(1, rep(0,nT), sigma)
    
    # data frame with expression, genotypes, and context information
    data_mat=melt(Y)
    colnames(data_mat) =  c("id", "Tissue", "E")
    data_mat=merge(x = data_mat, y = data.frame(id=1:N, G=genos), by = "id")
    
    # data frame with within-individual variability in expression, genotypes, and context information
    wY = Y - apply(Y,1,mean)
    wdata_mat=melt(wY)
    colnames(wdata_mat) =  c("id", "Tissue", "wE")
    data_mat=merge(x = wdata_mat, y = data_mat, by = c("id","Tissue"))
    
    CbC_lm_beta_and_se=data.frame(t(sapply(X = 1:nT, FUN = function(j) coef(summary(lm(Y[,j]~genos)))["genos",c("Estimate", "Std. Error")])), check.names = F)
    CbC_lm_beta_and_se[,"ci_l"]= CbC_lm_beta_and_se[,"Estimate"] - zval * CbC_lm_beta_and_se[,"Std. Error"]
    CbC_lm_beta_and_se[,"ci_u"]= CbC_lm_beta_and_se[,"Estimate"] + zval * CbC_lm_beta_and_se[,"Std. Error"]
    CbC_lm_beta_and_se=CbC_lm_beta_and_se[1:(nT+1),] 
    
    CbC_wlm_beta_and_se=data.frame(t(sapply(X = 1:nT, FUN = function(j) coef(summary(lm(wY[,j]~genos)))["genos",c("Estimate", "Std. Error")])), check.names = F)
    CbC_wlm_beta_and_se[,"ci_l"]=CbC_wlm_beta_and_se[,"Estimate"] - zval * CbC_wlm_beta_and_se[,"Std. Error"]
    CbC_wlm_beta_and_se[,"ci_u"]=CbC_wlm_beta_and_se[,"Estimate"] + zval * CbC_wlm_beta_and_se[,"Std. Error"]
    
    avg = matrix(apply(Y,1,mean), ncol = 1) # make matrix for average expression across tissues
    dimnames(avg) = list(1:N, c("T_AVG"))
    CbC_wlm_avg=data.frame(t(sapply(X = 1:1, FUN = function(j) coef(summary(lm(avg[,j]~genos)))["genos",c("Estimate", "Std. Error")])),check.names = F)
    CbC_wlm_avg[,"ci_l"]=CbC_wlm_avg[,"Estimate"] - zval * CbC_wlm_avg[,"Std. Error"]
    CbC_wlm_avg[,"ci_u"]=CbC_wlm_avg[,"Estimate"] + zval * CbC_wlm_avg[,"Std. Error"]
    
    
    #CbC_lm_beta_and_se[, "Context"] = c("Lung(-specific)", "Liver(-specific)", "Heart(-specific)", "Brain(-specific)", "Blood(-specific)", "Shared")
    #CbC_wlm_beta_and_se[, "Context"] = c("Lung(-specific)", "Liver(-specific)", "Heart(-specific)", "Brain(-specific)", "Blood(-specific)")
    CbC_lm_beta_and_se[, "Context"] = c("Lung", "Liver", "Heart", "Brain", "Blood", "Shared")
    CbC_wlm_beta_and_se[, "Context"] = c("Lung", "Liver", "Heart", "Brain", "Blood")
    CbC_wlm_avg[, "Context"]="Shared"
    
    my_plot_data=rbind(CbC_lm_beta_and_se %>% mutate(type="CxC"), 
                       CbC_wlm_beta_and_se%>% mutate(type="FastGxC"),
                       CbC_wlm_avg%>% mutate(type="FastGxC"))
    
    my_plot_data %<>% mutate(type=factor(x = type, levels = c("CxC","FastGxC")), scenario=Betas$scenario[i])
    
    my_plot_data_all_scenarios = rbind(my_plot_data_all_scenarios,my_plot_data)
  }
  
  my_plot_data_all_scenarios=my_plot_data_all_scenarios %>% 
    mutate(scenario = factor(x = scenario,
                             levels = c("No heterogeneity\n(No eQTL)", 
                                        "No heterogeneity", 
                                        "Single-context\nheterogeneity\nnoShared",
                                        "Single-context\nheterogeneity",
                                        "Extensive\nheterogeneity\nnoShared",
                                        "Extensive\nheterogeneity"),
                             labels = c("No heterogeneity\n(No eQTL)", 
                                        "No heterogeneity", 
                                        "Single-context\nheterogeneity\nnoShared",
                                        "Single-context\nheterogeneity",
                                        "Extensive\nheterogeneity\nnoShared",
                                        "Extensive\nheterogeneity")),
           
           Context = factor(x = Context, levels = c("Lung", "Liver", "Heart", "Brain", "Blood", "Shared"))
    )
  
  Fig1B=ggplot(my_plot_data_all_scenarios %>% 
                 filter(scenario!="Extensive\nheterogeneity\nnoShared") %>% 
                 droplevels(), 
               aes(y=Context, x=Estimate, xmin=ci_l, xmax=ci_u, color=Context))+ 
    facet_wrap(scenario~type, nrow = 1) + #, scales = "free_x"
    geom_point(size=3) + 
    geom_errorbarh(size = .5, height = .2, color = "gray50") +
    xlab("eQTL effect") + 
    geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5) +
    theme_bw()+
    theme(
      # panel.grid.major = element_blank(),
      # panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.title.x = element_text(size=20),
      axis.title.y = element_blank(),
      axis.text = element_text(size=13, color = 'black'),
      #strip.text.x = element_text(size = 10),
      strip.text.x = element_blank(),
      legend.position = "none",
      legend.text = element_text(size=15))+
    xlab("eQTL effect size") +
    scale_color_manual(name="", values=c("Shared"="black",
                                         "Liver" = "#598dc9", 
                                         "Lung" = "#5aaa46", 
                                         "Heart" = "#f2b342", 
                                         "Brain" = "#EF4869", 
                                         "Blood" = "#EA91BD"))+
    scale_x_continuous(labels = scales::number_format(accuracy = 0.1), 
                       breaks = scales::pretty_breaks(n = 2))
  #scale_x_continuous(breaks = c(0.0,1), limits = c(-0.7,1.2), labels = scales::number_format(accuracy = 0.1))
  
  ggsave(filename = 'Fig01_Method_Illustration/Fig1B.png', plot = Fig1B, width = 14, height = 3)  
  
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%% Figure 2 & Sup Figures: Simulation Study
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(1){
  param_setup = NULL
  
  #### Data preprocess 
  if(1){
    sim_res = fread(file = paste0("../Input_Files/Figure2_Simulation/simulation_results_typeIErrorPower_reviews.txt"), sep = "\t", data.table = F) %>%
      mutate(nC = paste0("#Contexts:", nC),
             N=paste0("#Individuals:", N),
             nC=factor(x = nC, levels=paste0("#Contexts:", c(8,49))),
             N=factor(x = N, levels=paste0("#Individuals:", c(100,698))),
             missing=factor(x = missing, levels=c("Complete data", "OneK1K missing", "GTEx missing")),
             method=factor(x = method, levels = c("CxC", "CxC-unadjasted", "LM-GxC", "LMM-GxC", "FastGxC",
                                                  "Metasoft_CxC_FE",
                                                  "Metasoft_CxC_RE",
                                                  "Metasoft_CxC_Mvalues_sigRE2",
                                                  "Metasoft_CxC_HetQ",
                                                  "Metasoft_CxC_Mvalues",
                                                  "Metasoft_CxC_HetQ_sigRE2",
                                                  "Metasoft_CxC_RE2",
                                                  "Metasoft_FastGxC_FE",
                                                  "Metasoft_FastGxC_RE",
                                                  "Metasoft_FastGxC_Mvalues_sigRE2",
                                                  "Metasoft_FastGxC_HetQ",
                                                  "Metasoft_FastGxC_Mvalues",
                                                  "Metasoft_FastGxC_HetQ_sigRE2",
                                                  "Metasoft_FastGxC_RE2",
                                                  "MTFERE2_FE",
                                                  "MTFERE2_RE",
                                                  "MTFERE2_Mvalues_sigRE2",
                                                  "MTFERE2_HetQ",
                                                  "MTFERE2_Mvalues",
                                                  "MTFERE2_HetQ_sigRE2",
                                                  "MTFERE2_RE2"))) #, "METASOFT", "Meta-Tissue"
    
    sim_res_bias = fread(file = "../Input_Files/Figure2_Simulation/simulation_results_effects_summarized_2025.txt", sep = "\t", data.table = F) %>% 
      mutate(nC = paste0("#Contexts:", nC),
             N=paste0("#Individuals:", N),
             nC=factor(x = nC, levels=paste0("#Contexts:", c(8,49))),
             N=factor(x = N, levels=paste0("#Individuals:",c(100,698))),
             missing=factor(x = missing, levels=c("Complete data", "OneK1K missing", "GTEx missing")),
             parspace=factor(x=parspace,levels=c("No Heterogeneity","Single-context Heterogeneity","Two context Heterogeneity","Weaker extensive heterogeneity", "Extensive heterogeneity")),
             shared=factor(x=shared,levels=c("No shared","Shared"))
      )
    
    run_time=read_csv("../Input_Files/Figure2_Simulation/run_time.csv") %>%
      mutate(model=factor(x = model, levels = c("TbT", "LM-GxE", "LMM-GxE", "FastGxE"), labels=c("CxC", "LM-GxC", "LMM-GxC", "FastGxC")))
    
  }
  
  main_methods= c("CxC", "LM-GxC", "LMM-GxC", "FastGxC", 
                  "Metasoft_CxC_FE",
                  "Metasoft_CxC_RE",
                  "Metasoft_CxC_Mvalues_sigRE2",
                  "Metasoft_CxC_HetQ",
                  "Metasoft_CxC_Mvalues",
                  "Metasoft_CxC_HetQ_sigRE2",
                  "Metasoft_CxC_RE2",
                  "Metasoft_FastGxC_FE",
                  "Metasoft_FastGxC_RE",
                  "Metasoft_FastGxC_Mvalues_sigRE2",
                  "Metasoft_FastGxC_HetQ",
                  "Metasoft_FastGxC_Mvalues",
                  "Metasoft_FastGxC_HetQ_sigRE2",
                  "Metasoft_FastGxC_RE2",
                  "MTFERE2_FE",
                  "MTFERE2_RE",
                  "MTFERE2_Mvalues_sigRE2",
                  "MTFERE2_HetQ",
                  "MTFERE2_Mvalues",
                  "MTFERE2_HetQ_sigRE2",
                  "MTFERE2_RE2")
  main_methods_noLM_GxC = c("CxC", "LMM-GxC", "FastGxC",
                            "Metasoft_CxC_FE",
                            "Metasoft_CxC_RE",
                            "Metasoft_CxC_Mvalues_sigRE2",
                            "Metasoft_CxC_HetQ",
                            "Metasoft_CxC_Mvalues",
                            "Metasoft_CxC_HetQ_sigRE2",
                            "Metasoft_CxC_RE2",
                            "Metasoft_FastGxC_FE",
                            "Metasoft_FastGxC_RE",
                            "Metasoft_FastGxC_Mvalues_sigRE2",
                            "Metasoft_FastGxC_HetQ",
                            "Metasoft_FastGxC_Mvalues",
                            "Metasoft_FastGxC_HetQ_sigRE2",
                            "Metasoft_FastGxC_RE2",
                            "MTFERE2_FE",
                            "MTFERE2_RE",
                            "MTFERE2_Mvalues_sigRE2",
                            "MTFERE2_HetQ",
                            "MTFERE2_Mvalues",
                            "MTFERE2_HetQ_sigRE2",
                            "MTFERE2_RE2")
  d=0.05
  pcs=4
  questions=c("Is there an eQTL?", "Is the eQTL context-specific?", "Does context c drive heterogeneity?")
  
  custom_theme <- theme(
    panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    legend.text = element_text(size=20), 
    legend.background = element_blank(),
    strip.text = element_text(size = 15),
    axis.title = element_text(size=20),
    axis.text = element_text(size=20, color = 'black')
  )
  
  #### Main Figure Under Single Context Heterogeneity  
  if(1){
    
    # A: Type I error rate to identify an eQTL and test if the eQTL is context-specific 
    if(TRUE){
      p_null = rbind(sim_res %>% filter(Tissues %in% c("global_eQTL") & 
                                          scenario_nr %in% c(1:10) &
                                          method %in% main_methods & 
                                          missing=="GTEx missing" & 
                                          N=="#Individuals:698" & 
                                          nC== "#Contexts:49" & 
                                          shared == "No shared") %>% 
                       mutate(parspace = questions[1]),
                     sim_res %>% filter(Tissues %in% c("global") & 
                                          scenario_nr %in% c(1:10) &
                                          method %in% main_methods & 
                                          missing=="GTEx missing" & 
                                          N=="#Individuals:698" & 
                                          nC== "#Contexts:49") %>% 
                       mutate(parspace = questions[2])) %>%
        mutate(parspace = factor(parspace, levels= c(questions[1:2]))) %>%
        droplevels() %>% 
        ggplot(aes(x=w_corr, y = PD, color = method, linetype = shared)) + 
        facet_grid(~parspace) +
        geom_point(size=pcs,position=position_dodge(d)) +
        geom_errorbar(aes(ymin=LCI, ymax=UCI), width=.1, linewidth = .5, 
                      position=position_dodge(d)) +
        geom_hline(yintercept = 0.05, linetype="dashed") + 
        geom_line(position=position_dodge(d)) +
        theme_bw() + 
        custom_theme + 
        theme(
          legend.position = "top",
          legend.title = element_blank())  +   
        scale_colour_manual(values = manuscript_colors_vec) +
        xlab(expression("Intra-individual correlation")) +
        ylab("Global Type I Error Rate")  + 
        scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
                           breaks = seq(0,0.2,.05), limits = c(0,0.22))  +
        scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8)) +
        guides(
          color = guide_legend(order = 1),    # Color legend appears first
          linetype = guide_legend(order = 2), # Linetype legend appears second
          shape = guide_legend(order = 3)     # Shape legend appears last
        ) 
      
      p_null
      
      Figure2_legend = get_legend(p_null)
      
      
    }
    
    ## B: Power to identify an eQTL and test if the eQTL is context-specific
    if(TRUE){
      p_het = rbind(
        sim_res %>% filter(Tissues == "global_eQTL" & scenario_nr %in% c(11:15) &
                             method %in% main_methods & missing=="GTEx missing" & 
                             N=="#Individuals:698" & nC== "#Contexts:49") %>% 
          mutate(parspace = questions[1]), 
        sim_res %>% filter(Tissues %in% c("global") & scenario_nr %in% c(11:20) &
                             method %in% main_methods & missing=="GTEx missing" & 
                             N=="#Individuals:698" & nC== "#Contexts:49") %>% 
          mutate(parspace = questions[2])) %>%
        mutate(parspace = factor(parspace, levels= c(questions[1:2]))) %>%
        droplevels() %>% 
        ggplot(aes(x=w_corr, y = PD, color = method, linetype = shared)) + 
        facet_grid(~parspace) +
        geom_point(size=pcs,position=position_dodge(d)) +
        geom_errorbar(aes(ymin=LCI, ymax=UCI), width=.1, linewidth = .5, position=position_dodge(d)) +
        geom_hline(yintercept = 0.80, linetype="dashed") + 
        geom_line(position=position_dodge(d)) +
        theme_bw() + custom_theme + 
        theme(
          legend.position = c(0.5, 0.5),
          legend.title = element_blank())  +   
        scale_colour_manual(values = manuscript_colors_vec, guide="none") +
        guides(linetype = "none") +  
        xlab(expression("Intra-individual correlation")) +
        ylab("Global Power")  + 
        scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(0,1,0.2), limits = c(0, 1)) +
        scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8))
      
      
      p_het
      
      
      
      
      
    }
    
    ## C: Marginal power to identify context driving the heterogeneity
    if(TRUE){
      d=0
      
      p_het_tbt=sim_res %>% filter(Tissues != "global" & scenario_nr %in% c(11:20) &
                                     method %in% main_methods & missing=="GTEx missing" & 
                                     N=="#Individuals:698" & nC== "#Contexts:49" & 
                                     Tissues %in% paste0("p.value_T",c(1,2,3,4,49))) %>%
        mutate(parspace = questions[3], shared=factor(shared)) %>%
        droplevels() %>%
        mutate(effect=ifelse(test = Tissues=="p.value_T49",
                             yes = "Most heterogeneous\ncontext(s)", 
                             no = "Least heterogeneous\ncontext(s)")) %>%
        
        ggplot(aes(x=w_corr, y = PD, color = method)) + 
        facet_grid(~parspace) +
        geom_point(aes(shape=effect),size=pcs, position=position_dodge(d)) +
        geom_line(aes(group = interaction(Tissues, shared), 
                      linetype = shared), position=position_dodge(d)) +
        theme_bw() + custom_theme + 
        theme(
          legend.position = c(0.25, 0.6),
          legend.title = element_blank())  +   
        scale_colour_manual(values = manuscript_colors_vec, guide="none") +
        xlab(expression("Intra-individual correlation")) +
        ylab("Marginal Power")  + 
        scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
                           breaks = seq(0,1,0.2), limits = c(0, 1)) +
        scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8)) + 
        geom_hline(yintercept = c(0.80), linetype="dashed")  +
        guides(linetype = "none") 
      
      p_het_tbt  
    }
    
    ## D: Projected run time for 250 individuals  
    if(TRUE){
      n_indi = 250
      
      ### get projected results for n_tests = 200e6
      n_tested_to_project = 200e6
      lm_results = run_time %>% 
        filter(n_individuals == n_indi) %>% 
        dplyr::select(n_tests,model,avg_days) %>%  
        group_by(model) %>% 
        summarise(avg_days=10^predict(object=lm(log10(avg_days) ~ log10(n_tests)), 
                                      newdata = data.frame(n_tests = c(n_tested_to_project)))) %>% 
        mutate(n_tests=n_tested_to_project)  %>% dplyr::select(model,n_tests,avg_days)
      
      
      pd=position_dodge(0.1)
      
      Figure2D = run_time %>% 
        filter(n_individuals == n_indi) %>% dplyr::select(n_tests,model,avg_days) %>% 
        ggplot(aes(x = log10(n_tests),y = avg_days, color = model)) +
        geom_point(position = pd,size=pcs) +
        geom_line(aes(group = model), position = pd) +
        #### add projected point
        geom_point(data=lm_results,aes(x = log10(n_tests),y = avg_days, color = model), size=pcs, position = pd) +
        geom_line(data=bind_rows(lm_results,run_time %>% filter(n_individuals == n_indi) %>% 
                                   dplyr::select(n_tests,model,avg_days) %>% 
                                   filter(n_tests==1e7)),
                  aes(group=model), linetype="dashed",position = pd) +
        ### add label for time @ projected points
        geom_label_repel(data=lm_results %>% 
                           mutate(time_label = case_when(model == "FastGxC" ~ paste0(round(avg_days*24*60*60)," secs"),
                                                         model == "CxC" ~ paste0(round(avg_days*24*60*60)," secs"),
                                                         model == "LM-GxC" ~ paste0(round(avg_days/30)," months"),
                                                         model == "LMM-GxC" ~ paste0(round(avg_days/365)," years"))),
                         aes(x = log10(n_tests),y = avg_days,label=time_label,color=model),
                         point.padding = 3, size=6, show.legend = FALSE, segment.alpha = 0, 
                         arrow = arrow(length = unit(0.015, "npc")))+
        scale_y_continuous("Run time in days (log)",trans = log10_trans(),
                           breaks = trans_breaks("log10", function(x) 10^x),
                           labels = trans_format("log10", math_format(10^.x))) +
        scale_x_continuous("Number of tests (log)", breaks = seq(2,9), labels = math_format(10^.x), limits = c(1.8,9)) +
        theme_bw() +
        custom_theme +
        theme(legend.position="none"
              # axis.title = element_text(size=20,color="black"),
              # axis.text = element_text(size=20,color="black"),
              # axis.line = element_line(colour = "black"),
              # plot.title = element_text(hjust = 0.5,size=15),
              # legend.title=element_text(size=15),
              # legend.text=element_text(size=23)
        ) +
        methods_color_scale +
        xlab("Number of Tests") +
        labs(color='') 
    }
    
    ## Combined Figure 2
    if(TRUE){
      Figure2=plot_grid(
        Figure2_legend, 
        ggdraw() +
          draw_plot(plot = p_null + theme(legend.position = "none"), 
                    x = 0.0, y = 0.50, width = 0.50, height = 0.50) +
          draw_plot(plot = p_het, 
                    x = 0.5, y = 0.50, width = 0.50, height = 0.50) +
          draw_plot(plot = p_het_tbt, 
                    x = 0.0, y = 0.00, width = 0.50, height = 0.50) +
          draw_plot(plot = Figure2D, 
                    x = 0.5, y = 0.00, width = 0.50, height = 0.50) +
          draw_plot_label(label = c('A','B','C','D'),  
                          x = c(0,0.5,0,.5), 
                          y = c(1,1,0.5,0.5), 
                          size = 30),
        rel_heights = c(0.05, 1), 
        ncol = 1)
      
      ggsave(plot = Figure2, filename = paste0('Fig02_Simulations',param_setup,'.pdf'),height = 15, width = 15)
    }
    
    ## manuscript claims 
    if(TRUE){
      # FastGxC is this much more powerful than existing methods
      
      
      power_increase = sim_res %>% filter(method %in% c("FastGxC", "CxC") & Tissues == "global_eQTL") %>% dplyr::select(method, w_corr, scenario_nr, N, nC, missing, PD) %>% pivot_wider(names_from = method, values_from = PD) %>%
        mutate(ratio = FastGxC/CxC)
      max_power_increase = max(power_increase$ratio)
      
      speed_increase = run_time %>% filter(model %in% c("FastGxC", "LMM-GxC")) %>% dplyr::select(n_individuals, n_tissues, n_tests, model, avg_years) %>% pivot_wider(names_from = model, values_from = avg_years) %>%
        mutate(ratio = `LMM-GxC`/FastGxC)
      speed_increase = max(speed_increase$ratio)
      
      ## manuscript claims - 
      # FastGxC is this much faster than existing methods:
      FastGxC_speed = run_time %>% 
        filter(n_individuals == n_indi & n_tests == 10000000 & model == "FastGxC") %>% dplyr::select(n_tests,model,avg_days) 
      LMM_GxC_speed = run_time %>% 
        filter(n_individuals == n_indi & n_tests == 10000000 & model == "LMM-GxC") %>% dplyr::select(n_tests,model,avg_days)
      speed_increase = LMM_GxC_speed[,"avg_days"]/FastGxC_speed[,"avg_days"]
    }    
    
    
  }
  
  # Sup Figure S1: Global Type I Error Rate 
  setwd("reviews/")
  if(1){
    
    ggsave(filename="FigureS01_global_T1ErRate.png",
           plot =  ggplot(sim_res %>% 
                            filter(Tissues == "global_eQTL" & scenario_nr %in% c(1:5) &
                                     (method %in% main_methods) ) %>% #
                            droplevels(), 
                          aes(x=w_corr, y = PD, color = method)) + 
             geom_point(aes(shape = method), size=3, position=position_dodge(d)) + 
             geom_errorbar(aes(ymin=LCI, ymax=UCI), width=.1, position=position_dodge(d)) +
             geom_line(position=position_dodge(d)) +
             theme_bw() + 
             custom_theme  +   
             xlab("Within-individual correlation") +
             ylab("Global Type I Error Rate to Identify eQTL")  + 
             scale_x_continuous(breaks = seq(0,.8,.2)) + 
             scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
             geom_hline(yintercept = 0.05, linetype="dashed") +
             methods_color_scale  + 
             methods_shape_scale + 
             facet_grid(N+nC~missing) +
             theme(
               legend.position = "top",
               legend.title = element_blank()), 
           width = 15, height = 17)
  }
  
  # Sup Figure S2: Global Type I Error Rate GxC
  if(1){
    alpha = 0.05
    
    ggsave(filename="FigureS02_global_T1ErRate_gxc.png",
           plot =  ggplot(sim_res %>% 
                            filter(Tissues == "global" & scenario_nr %in% c(1:10) &
                                     (method %in% main_methods) ) %>% #
                            droplevels(), 
                          aes(x=w_corr, y = PD, color = method, linetype=shared)) + 
             geom_point(aes(shape = method), size=3, position=position_dodge(d)) + 
             geom_errorbar(aes(ymin=LCI, ymax=UCI), width=.1, position=position_dodge(d)) +
             geom_line(position=position_dodge(d)) +
             theme_bw() + 
             custom_theme  +   
             xlab("Within-individual correlation") +
             ylab("Global Type I Error Rate to Identify if eQTL is Context-Specific")  + 
             scale_x_continuous(breaks = seq(0,.8,.2)) + 
             scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
             geom_hline(yintercept = alpha, linetype="dashed") +
             methods_color_scale  + 
             methods_shape_scale + 
             facet_grid(N+nC~missing) +
             theme(
               legend.position = "top",
               legend.title = element_blank(),
               legend.text = element_text(size=8)), 
           width = 12, height = 15)
  }
  
  # Sup Figure S3: Global Power - shared only
  if(1){
    alpha = 0.80
    ggsave(filename="FigureS03_global_Power_OnlyShared.png",
           plot =  ggplot(sim_res %>% 
                            filter(Tissues == "global_eQTL" & 
                                     scenario_nr %in% c(6:10) & 
                                     N=="#Individuals:100" &
                                     (method %in% main_methods) ) %>%
                            droplevels() %>%
                            mutate(shared=factor(shared, levels = "Shared", labels="Shared-only")), 
                          aes(x=w_corr, y = PD, color = method)) + 
             geom_point(aes(shape = method), size=3, position=position_dodge(d)) + 
             geom_errorbar(aes(ymin=LCI, ymax=UCI), width=.1, position=position_dodge(d)) +
             geom_line(position=position_dodge(d)) +
             theme_bw() + 
             custom_theme  +   
             xlab("Within-individual correlation") +
             ylab("Global Power to Identify eQTL")  + 
             scale_x_continuous(breaks = seq(0,.8,.2)) + 
             scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(0,1,0.2), limits = c(0, 1)) +
             geom_hline(yintercept = .8, linetype="dashed") + 
             methods_color_scale  + 
             methods_shape_scale + 
             facet_grid(shared+nC~missing) +
             theme(
               legend.position = "top",
               legend.title = element_blank()), 
           width = 12, height = 10)  
  }  
  
  # Sup Figure S4: Global Power   ## Single context heterogeneity 
  if(1){
    alpha = 0.80
    ggsave(filename="FigureS04_global_Power_singleContext.png",
           plot =  ggplot(sim_res %>% 
                            #mutate(par_space2=gsub(pattern="heterogeneity", replacement = "spec", x = par_space2)) %>% 
                            filter(Tissues == "global_eQTL" & scenario_nr %in% c(11:20) & N!="#Individuals:698" &
                                     (method %in% main_methods_noLM_GxC) ) %>%
                            droplevels(), 
                          aes(x=w_corr, y = PD, color = method, linetype = shared)) + 
             geom_point(aes(shape = method),size=3, position=position_dodge(d)) + 
             geom_errorbar(aes(ymin=LCI, ymax=UCI), width=.1, position=position_dodge(d)) +
             geom_line(position=position_dodge(d)) +
             theme_bw() + custom_theme  +   
             xlab("Within-individual correlation") +
             ylab("Global Power to Identify eQTL")  + 
             scale_x_continuous(breaks = seq(0,.8,.2)) + 
             scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(0,1,0.2), limits = c(0, 1)) +
             geom_hline(yintercept = .8, linetype="dashed") + 
             methods_color_scale_noLMGxC  + 
             methods_shape_scale_noLMGxC + 
             facet_grid(N+nC~missing) +
             theme(
               legend.position = "top",
               legend.title = element_blank()), 
           width = 12, height = 10)  
  }  
  
  # Sup Figure S5: Global Power   ## Extensive heterogeneity 
  if(1){
    alpha = 0.80
    ggsave(filename="FigureS05_global_Power_ExtensiveHet.png",
           plot =  ggplot(sim_res %>% 
                            filter(Tissues == "global_eQTL" & scenario_nr %in% c(36:40) & N!="#Individuals:698" &
                                     method %in% main_methods_noLM_GxC) %>%
                            droplevels(), 
                          aes(x=w_corr, y = PD, color = method)) + 
             geom_point(aes(shape = method),size=3, position=position_dodge(d)) + 
             geom_errorbar(aes(ymin=LCI, ymax=UCI), width=.1, position=position_dodge(d)) +
             geom_line(position=position_dodge(d)) +
             theme_bw() + custom_theme + 
             xlab("Within-individual correlation") +
             ylab("Global Power to Identify eQTL")  + 
             scale_x_continuous(breaks = seq(0,.8,.2)) + 
             geom_hline(yintercept = .8, linetype="dashed") + 
             scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(0,1,0.2), limits = c(0, 1)) +
             methods_color_scale_noLMGxC + 
             methods_shape_scale_noLMGxC + 
             facet_grid(N+nC~missing) +
             theme(
               legend.position = "top",
               legend.title = element_blank()), 
           width = 12, height = 10)  
  }
  
  # Sup Figure S6: Global Power GxC  ## Single context heterogeneity 
  if(1){
    alpha = 0.80
    ggsave(filename="FigureS06_global_Power_singleContext_gxc.png",
           plot =  ggplot(sim_res %>% 
                            #mutate(par_space2=gsub(pattern="heterogeneity", replacement = "spec", x = par_space2)) %>% 
                            filter(Tissues == "global" & scenario_nr %in% c(11:20) &
                                     (method %in% main_methods) ) %>%
                            droplevels(), 
                          aes(x=w_corr, y = PD, color = method)) + 
             geom_point(aes(shape = method), size=3, position=position_dodge(d)) + 
             geom_errorbar(aes(ymin=LCI, ymax=UCI), width=.1, position=position_dodge(d)) +
             geom_line(aes(linetype = shared),position=position_dodge(d)) +
             theme_bw() + custom_theme  +   
             xlab("Within-individual correlation") +
             ylab("Global Power to Test if eQTL is Context-Specific")  + 
             scale_x_continuous(breaks = seq(0,.8,.2)) + 
             scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(0,1,0.2), limits = c(0, 1)) +
             geom_hline(yintercept = .8, linetype="dashed") + 
             methods_color_scale  + 
             methods_shape_scale +
             facet_grid(N+nC~missing) +
             theme(
               legend.position = "top",
               legend.title = element_blank(),
               legend.text = element_text(size=8)), 
           width = 12, height = 15)  
  }  
  
  # Sup Figure S7: Global Power GxC## Two context heterogeneity
  if(1){
    alpha = 0.80
    ggsave(filename="FigureS07_global_Power_TwoContext_gxc.png",
           plot =  ggplot(sim_res %>% 
                            #mutate(par_space2=gsub(pattern="heterogeneity", replacement = "spec", x = par_space2)) %>% 
                            filter(Tissues == "global" & scenario_nr %in% c(21:30) &
                                     (method %in% main_methods) ) %>%
                            droplevels(), 
                          aes(x=w_corr, y = PD, color = method)) + 
             geom_point(aes(shape = method), size=3, position=position_dodge(d)) + 
             geom_errorbar(aes(ymin=LCI, ymax=UCI), width=.1, position=position_dodge(d)) +
             geom_line(aes(linetype=shared),position=position_dodge(d)) +
             theme_bw() + custom_theme  +   
             xlab("Within-individual correlation") +
             ylab("Global Power to Test if eQTL is Context-Specific")  + 
             scale_x_continuous(breaks = seq(0,.8,.2)) + 
             scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(0,1,0.2), limits = c(0, 1)) +
             geom_hline(yintercept = .8, linetype="dashed") + 
             methods_color_scale  + 
             methods_shape_scale + 
             facet_grid(N+nC~missing) +
             theme(
               legend.position = "top",
               legend.title = element_blank(),
               legend.text = element_text(size=8)), 
           width = 12, height = 15)  
  }  
  
  # Sup Figure S8: Global Power GxC ## Extensive heterogeneity 
  if(1){
    alpha = 0.80
    d=0.05
    ggsave(filename="FigureS08_global_Power_ExtensiveHet_gxc.png",
           plot =  ggplot(sim_res %>% 
                            filter(Tissues == "global" & scenario_nr %in% c(36:40) &
                                     method %in% main_methods) %>%
                            droplevels(), 
                          aes(x=w_corr, y = PD, color = method, shape = N)) + 
             geom_point(size=3, position=position_dodge(d)) + 
             geom_errorbar(aes(ymin=LCI, ymax=UCI), width=.1, position=position_dodge(d)) +
             geom_line(position=position_dodge(d)) +
             theme_bw() + custom_theme + 
             xlab("Within-individual correlation") +
             ylab("Global Power to Test if eQTL is Context-Specific")  + 
             scale_x_continuous(breaks = seq(0,.8,.2)) + 
             geom_hline(yintercept = .8, linetype="dashed") + 
             scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(0,1,0.2), limits = c(0, 1)) +
             methods_color_scale + 
             #methods_shape_scale + 
             facet_grid(nC~missing) +
             theme(
               legend.position = "top",
               legend.title = element_blank()), 
           width = 15, height = 15)  
  }  
  
  # Sup Figure S9: Marginal Type I Error Rate 
  if(1){
    alpha = 0.05
    ggsave(filename="FigureS09_marginal_T1ErRate.png",
           plot =  ggplot(sim_res %>% 
                            filter(Tissues != "global" & scenario_nr %in% c(1:10) &
                                     (method %in% "FastGxC") & Tissues %in% paste0("p.value_T",seq(1,49))) %>% #
                            droplevels(), 
                          aes(x=as.character(w_corr), y = PD, color = method, linetype = shared)) + 
             facet_grid(N+nC~missing) +
             geom_point(aes(group = Tissues), size=3, position=position_dodge(d)) + 
             geom_errorbar(aes(ymin=LCI, ymax=UCI, group = Tissues), width=.1, position=position_dodge(d)) +
             geom_line(aes(group = interaction(Tissues,shared)), position=position_dodge(d)) +
             theme_bw() + 
             custom_theme  + 
             xlab("Within-individual correlation") +
             ylab("Marginal Type I Error Rate to Test if Context c Drives Heterogeneity")  + 
             scale_x_discrete(breaks = seq(0,.8,.2)) + 
             scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(0,.5,.01), limits = c(0,0.05))+
             geom_hline(yintercept = alpha, linetype="dashed") +
             methods_color_scale  + 
             theme(legend.position = "top",
                   legend.title = element_blank()), width = 15, height = 15)  
  }
  
  # Sup Figure S10: Marginal Power ### Single-Context Heterogeneity
  if(1){
    alpha = 0.8
    d=0
    ggsave(filename="FigureS10_marginal_Power_singleContext.png",
           plot =  sim_res %>% filter(Tissues != "global" & scenario_nr %in% c(11:20) &
                                        Tissues %in% paste0("p.value_T",1:49)) %>%
             mutate(parspace = questions[3], shared=factor(shared)) %>%
             droplevels() %>%
             mutate(effect=ifelse(
               test = (nC == "#Contexts:49" & Tissues=="p.value_T49")|(nC == "#Contexts:8" & Tissues=="p.value_T8"),
               
               yes = "Most heterogeneous context(s)", 
               no = "Least heterogeneous context(s)")) %>%
             
             ggplot(aes(x=w_corr, y = PD, color = method)) + 
             facet_grid(N+nC~missing) +
             geom_point(aes(shape=effect),size=pcs, position=position_dodge(d)) +
             geom_line(aes(group = interaction(Tissues, shared), 
                           linetype = shared), position=position_dodge(d)) +
             theme_bw() + custom_theme + 
             theme(
               legend.position = "top",
               legend.title = element_blank())  +   
             scale_colour_manual(values = manuscript_colors_vec, guide="none") +
             xlab(expression("Intra-individual correlation")) +
             ylab("Marginal Power to Test if Context c Drives Heterogeneity")  + 
             scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
                                breaks = seq(0,1,0.2), limits = c(0, 1)) +
             scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8)) + 
             geom_hline(yintercept = alpha, linetype="dashed")  , 
           width = 15, height = 15)  
  }
  
  # Sup Figure S11: Marginal Power ### Two-Context Heterogeneity
  if(1){
    alpha = 0.8
    d=.01
    ggsave(filename="FigureS11_marginal_Power_TwoContexts.png",
           plot =  sim_res %>% filter(Tissues != "global" & scenario_nr %in% c(21:30) &
                                        Tissues %in% paste0("p.value_T",1:49)) %>%
             mutate(parspace = questions[3], shared=factor(shared)) %>%
             droplevels() %>%
             mutate(effect=ifelse(
               test = 
                 (nC == "#Contexts:49" & Tissues %in% c("p.value_T48","p.value_T49"))| 
                 (nC == "#Contexts:8" & Tissues%in%c("p.value_T7","p.value_T8")),
               
               
               yes = "Most heterogeneous context(s)", 
               no = "Least heterogeneous context(s)")) %>%
             
             ggplot(aes(x=w_corr, y = PD, color = method)) + 
             facet_grid(N+nC~missing) +
             geom_point(aes(shape=effect),size=pcs, position=position_dodge(d)) +
             geom_line(aes(group = interaction(Tissues, shared), 
                           linetype = shared), position=position_dodge(d)) +
             theme_bw() + custom_theme + 
             theme(
               legend.position = "top",
               legend.title = element_blank())  +   
             scale_colour_manual(values = manuscript_colors_vec, guide="none") +
             xlab(expression("Intra-individual correlation")) +
             ylab("Marginal Power to Test if Context c Drives Heterogeneity")  + 
             scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
                                breaks = seq(0,1,0.2), limits = c(0, 1)) +
             scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8)) + 
             geom_hline(yintercept = alpha, linetype="dashed"),
           width = 15, height = 15)  
  }
  
  #Sup Figure S12: Marginal power - Extensive Heterogeneity
  if(1){
    d=0
    ggsave(filename="FigureS12_marginal_Power_ExtensiveHet.png",
           plot =  sim_res %>% filter(Tissues != "global" & scenario_nr %in% c(36:40) &
                                        Tissues %in% paste0("p.value_T",1:49)) %>%
             mutate(parspace = questions[3], shared=factor(shared)) %>%
             droplevels() %>%
             ggplot(aes(x=w_corr, y = PD, color = method)) + 
             facet_grid(N+nC~missing) +
             geom_point(size=pcs, position=position_dodge(d)) +
             geom_line(aes(group = interaction(Tissues, shared), 
                           linetype = shared), position=position_dodge(d)) +
             theme_bw() + custom_theme + 
             theme(
               legend.position = "top",
               legend.title = element_blank())  +   
             scale_colour_manual(values = manuscript_colors_vec, guide="none") +
             xlab(expression("Intra-individual correlation")) +
             ylab("Marginal Power to Test if Context c Drives Heterogeneity")  + 
             scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
                                breaks = seq(0,1,0.2), limits = c(0, 1)) +
             scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8)) + 
             geom_hline(yintercept = alpha, linetype="dashed")  +
             guides(linetype = "none"), 
           width = 15, height = 15)  
  }
  
  # Sup Figure S13-14: Effect size comparison
  if(1){
    
    alpha = 0.05
    pcs = 1
    d = 1
    data_s13 = sim_res_bias %>% 
      filter(method == "FastGxC" & nC=="#Contexts:8" & 
               N=="#Individuals:698" & w_corr %in% c(0.8) & scenario_nr %in% c(11:20,36:40)) %>% 
      mutate(context = gsub("T", "", context)) %>% mutate(context = factor(context, levels=as.character(seq(1,49)))) %>% 
      mutate(w_corr = paste0("Within individual correlation:", w_corr)) %>% 
      pivot_longer(cols = c(mean_effect, betas), 
                   names_to = "effect_type", 
                   values_to = "effect_size") %>%
      
      # Rename effect_type for legend labels
      mutate(effect_type = recode(effect_type, 
                                  "mean_effect" = "Estimated Effect",
                                  "betas" = "True Effect"))
    
    ggsave(filename="FigureS13_bias.png",
           plot = ggplot(data_s13, aes(x=context, y = effect_size, color = effect_type)) + 
             facet_grid(parspace+shared~missing) +
             geom_point(data = data_s13 %>% filter(effect_type == "Estimated Effect"), size=pcs,position=position_dodge(d)) +
             geom_errorbar(data = data_s13 %>% filter(effect_type == "Estimated Effect"), aes(ymin=LCI, ymax=UCI), width=0.2, size = 1, position=position_dodge(d))+
             geom_segment(data = data_s13 %>% filter(effect_type == "True Effect"), aes(x = as.numeric(context) + 0.1, xend = as.numeric(context) + 0.3, 
                                                                                        y = effect_size, yend = effect_size), color = "red", size = 1)+
             
             #geom_hline(yintercept = 0, col="red")  + theme(legend.position = "none")  + 
             ylab("eQTL effect size with 95% confidence interval") + 
             #scale_y_continuous(limits = c(-0.05,0.05))  + 
             scale_color_manual(values = c("True Effect" = "red", "Estimated Effect" = "#636363"))+
             theme_bw()+
             theme(legend.position = "top",
                   legend.title=element_blank(),
                   panel.spacing = unit(2, "lines"),
                   axis.text.x = element_text(size=20)) +
             custom_theme + 
             xlab("Context"),
           width = 15, height = 15)
  }
  
  # Sup Figure S12: Run time for 1000 individuals
  if(1){
    n_indi = 1000
    
    ### get projected results for n_tests = 200e6
    n_tested_to_project = 200e6
    lm_results = run_time %>% 
      filter(n_individuals == n_indi) %>% 
      dplyr::select(n_tests,model,avg_days) %>%  
      group_by(model) %>% 
      summarise(avg_days=10^predict(object=lm(log10(avg_days) ~ log10(n_tests)), 
                                    newdata = data.frame(n_tests = c(n_tested_to_project)))) %>% 
      mutate(n_tests=n_tested_to_project)  %>% dplyr::select(model,n_tests,avg_days)
    
    ggsave(filename="FigureS12_runtime_1000.png",
           plot = run_time %>% 
             filter(n_individuals == n_indi) %>% dplyr::select(n_tests,model,avg_days) %>% 
             ggplot(aes(x = log10(n_tests),y = avg_days, color = model)) +
             geom_point(position = pd,size=pcs) +
             geom_line(aes(group = model), position = pd) +
             #### add projected point
             geom_point(data=lm_results,aes(x = log10(n_tests),y = avg_days, color = model), size=pcs, position = pd) +
             geom_line(data=bind_rows(lm_results,run_time %>% filter(n_individuals == n_indi) %>% 
                                        dplyr::select(n_tests,model,avg_days) %>% 
                                        filter(n_tests==1e7)),
                       aes(group=model), linetype="dashed",position = pd) +
             ### add label for time @ projected points
             geom_label_repel(data=lm_results %>% 
                                mutate(time_label = case_when(model == "FastGxC" ~ paste0(round(avg_days*24*60*60)," secs"),
                                                              model == "CxC" ~ paste0(round(avg_days*24*60*60)," secs"),
                                                              model == "LM-GxC" ~ paste0(round(avg_days/30)," months"),
                                                              model == "LMM-GxC" ~ paste0(round(avg_days/365)," years"))),
                              aes(x = log10(n_tests),y = avg_days,label=time_label,color=model),
                              point.padding = 3, size=6, show.legend = FALSE, segment.alpha = 0, 
                              arrow = arrow(length = unit(0.015, "npc")))+
             scale_y_continuous("Run time in days (log)",trans = log10_trans(),
                                breaks = trans_breaks("log10", function(x) 10^x),
                                labels = trans_format("log10", math_format(10^.x))) +
             scale_x_continuous("Number of tests (log)", breaks = seq(2,9), labels = math_format(10^.x), limits = c(1.8,9)) +
             theme_bw() +
             custom_theme +
             theme(legend.position="none"
                   # axis.title = element_text(size=20,color="black"),
                   # axis.text = element_text(size=20,color="black"),
                   # axis.line = element_line(colour = "black"),
                   # plot.title = element_text(hjust = 0.5,size=15),
                   # legend.title=element_text(size=15),
                   # legend.text=element_text(size=23)
             ) +
             methods_color_scale +
             xlab("Number of Tests") +
             labs(color=''), width = 5, height = 5)
  }
  
}

################################################################
################## Figure 3 with eGenes and supp Fig 16 ########
################################################################

if(1){
  # 3A barplot of sharing vs specificity
  # make input data Bulk: DONE
  if(1){
    
    ## create eGene input file
    if(1){
      egenes_gtex = read_csv("../Input_Files/Figure3_Performance/eGenes.v8.EUR.all_tissues.residualized_exp_types.txt") %>%
        mutate(exp_type = case_when(exp_type == "normalized_and_residualized_expression" ~ "CxC",
                                    exp_type == "normalized_and_residualized_expression_heterogeneous" ~ "Specific",
                                    exp_type == "normalized_and_residualized_expression_homogeneous" ~ "Shared",
                                    TRUE ~ exp_type)) %>% mutate(cohort = "GTEx")
      egenes_sc = read_csv("../Input_Files/Figure3_Performance/eGenes.scMeta.all_contexts.residualized_exp_types.txt") %>%
        mutate(exp_type = case_when(exp_type == "mean_norm_res_exp" ~ "CxC",
                                    exp_type == "mean_norm_res_exp.specific" ~ "Specific",
                                    exp_type == "mean_norm_res_exp.shared" ~ "Shared",
                                    TRUE ~ exp_type)) %>% mutate(cohort = "SC")
      all_egenes = rbind(egenes_gtex, egenes_sc)
      
      within_tis = all_egenes %>% 
        filter(exp_type != "CxC" ) ## only fastgxc
      
      counted = within_tis %>% group_by(gene, cohort) %>% summarize(n_exp_types = n_distinct(exp_type), 
                                                                    group = case_when(n_distinct(exp_type) == 2 ~ "ts&sh",
                                                                                      n_distinct(exp_type) == 1 & unique(exp_type) == "Shared" ~ "sh",
                                                                                      n_distinct(exp_type) == 1 & unique(exp_type) == "Specific" ~ "sp",
                                                                                      TRUE ~ NA_character_)) %>% distinct %>% 
        group_by(cohort, group) %>% mutate(value = n()) %>% dplyr::select(c("group", "value", "cohort")) %>% distinct %>%
        group_by(cohort) %>% mutate(percentage = value/sum(value)) %>% mutate(group = fct_relevel(group,c("sp", "sh", "ts&sh")))
      
      fwrite(counted, file = "../Input_Files/Figure3_Performance/eGenes_counted.csv", sep = ",")
    }
    
    counted = fread("../Input_Files/Figure3_Performance/eGenes_counted.csv", sep = ",", data.table = F)
    fig3A = ggplot(counted, aes(fill=group, y=percentage, x=cohort)) + 
      geom_bar_pattern(aes(pattern = group), position="stack", stat="identity", color = "black", width = 0.7) +
      scale_pattern_manual(values = c("none", "stripe", "weave"), labels = c( "sh" = "Shared",
                                                                              "sp" = "Specific",
                                                                              "ts&sh" = "Both")) +
      scale_fill_manual(values = c("white", "white", "white"), labels = c("sh" = "Shared", 
                                                                          "sp" = "Specific",
                                                                          "ts&sh" = "Both")) + 
      scale_x_discrete(labels = c("GTEx" = "Tissues", 
                                  "SC" = "PBMC"))+
      guides( pattern = guide_legend(),
              fill = guide_legend())+
      theme_bw() + 
      theme(legend.position="top",
            axis.title = element_text(face = "bold", size=25,color="black"),
            axis.text = element_text(size=20,color="black"),
            axis.text.x = element_text(face = "bold"),
            axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            plot.title = element_text(hjust = 0.5,size=15),
            legend.title=element_blank(),
            legend.text=element_text(size=18),
            strip.text.y = element_text(size = 0),
            strip.background = element_blank(),
            panel.spacing = unit(2, "lines")) + 
      xlab("") + ylab("Proportion of eGenes")
    
    
  }
  
  ######### manuscript claims ######
  
  
  # 3B 
  # make input data Bulk: DONE
  if(1){
    if(1){
      eGenes = read_csv("../Input_Files/Figure3_Performance/eGenes.v8.EUR.all_tissues.residualized_exp_types.txt")
      eGenes_counted = eGenes %>% group_by(exp_type, tissue) %>% summarize(n_egenes = n(), .groups = "drop")
      
      fastgxe_tissue_nsamples = read_csv("../Input_Files/Figure3_Performance/GTEx_v8_NSamples_by_Tissue_and_Method.csv") %>% 
        filter(method != ".v8.EUR.normalized_and_residualized_expression") %>% 
        dplyr::select(tissue, n_samples)
      
      fastgxe_eGenes_counted = eGenes_counted %>% 
        mutate(exp_type = case_when(exp_type == "normalized_and_residualized_expression" ~ "tbt",
                                    exp_type == "normalized_and_residualized_expression_homogeneous" ~ "fastgxe",
                                    exp_type == "normalized_and_residualized_expression_heterogeneous" ~ "fastgxe",
                                    T ~ "unknown?")) %>% 
        filter(exp_type == "fastgxe") %>% 
        dplyr::select(tissue, n_egenes) %>% 
        left_join(fastgxe_tissue_nsamples, by=c("tissue")) %>% 
        mutate(is_avg = case_when(tissue == "AverageTissue" ~ "Shared", T ~ "Tissue-Specific"))
      
      fastgxe_eGenes_counted %>% write_csv("../Input_Files/Figure3_Performance/figure3B_GTEx_eGenes_data.csv")
    }
    
    ## make input file single cell with number of cells: DONE
    if(0){
      sc_eGenes = read_csv("../Input_Files/Figure3_Performance/eGenes.scMeta.all_contexts.residualized_exp_types.txt", col_types = cols(.default = col_character()))
      sc_eGenes_counted = sc_eGenes %>% 
        group_by(exp_type, tissue) %>% 
        summarize(n_egenes = n(), .groups="drop")
      
      fastgxe_sc_tissue_ncells = read_csv("../Input_Files/Figure3_Performance/CLUES_OneK1K_NCells_by_Context.csv")
      
      fastgxe_sceGenes_counted = sc_eGenes_counted %>% 
        mutate(exp_type = case_when(exp_type == "mean_norm_res_exp" ~ "tbt",
                                    exp_type == "mean_norm_res_exp.shared" ~ "fastgxe",
                                    exp_type == "mean_norm_res_exp.specific" ~ "fastgxe",
                                    T ~ "unknown?")) %>% 
        filter(exp_type == "fastgxe") %>% 
        dplyr::select(tissue, n_egenes) %>% 
        merge(fastgxe_sc_tissue_ncells, by.x=c("tissue"), by.y = "final_celltypes") %>% 
        mutate(is_avg = case_when(tissue == "AverageContext" ~ "Shared", T ~ "Tissue-Specific"))
      
      fastgxe_sceGenes_counted %>% write_csv("../Input_Files/Figure3_Performance/figure3B_sc_eGenes_data.csv")
    }
  }
  
  fastgxe_eGenes_counted = read_csv("../Input_Files/Figure3_Performance/figure3B_GTEx_eGenes_data.csv")
  fastgxe_sceGenes_counted = read_csv("../../manuscript/Input_Files/Figure3_Performance/figure3B_sc_eGenes_data.csv")
  fastgxe_eGenes_counted$group = "Tissues"
  fastgxe_sceGenes_counted$group = "PBMC"
  
  
  # 3C sharing and specificity of GTEx eGenes across tissues ##### made in powerpoint
  if(1){
    
    # plot: number of eqtls that are found in x number of tissues plot
    if(1){
      
      # make input file Bulk: DONE
      if(0){
        eGenes = read_csv("/u/project/bballiu/bballiu/FastGxC/FastGxC_Manuscript/results/eGenes.v8.EUR.all_tissues.residualized_exp_types.txt")
        
        eGenes_tissues_counted = eGenes %>% filter(exp_type == "normalized_and_residualized_expression_heterogeneous") %>%
          mutate(exp_type = case_when(exp_type == "normalized_and_residualized_expression" ~ "tbt",
                                      exp_type == "normalized_and_residualized_expression_homogeneous" ~ "fastgxe",
                                      exp_type == "normalized_and_residualized_expression_heterogeneous" ~ "fastgxe",
                                      T ~ "unknown?"))  %>% 
          filter(exp_type=="fastgxe") %>% 
          dplyr::select(gene, tissue) %>% distinct %>% 
          group_by(gene) %>% 
          summarize(n_tissues = n(), .groups="drop") %>% 
          dplyr::select(gene, n_tissues)
        
        counted2 = eGenes_tissues_counted %>% 
          dplyr::select(n_tissues, gene) %>% 
          group_by(n_tissues) %>% 
          summarize(n_egenes = n(), .groups="drop") %>% 
          dplyr::select(n_tissues, n_egenes)
        
        counted2 %>% write_csv("manuscript/input_files/Figure3_Performance/figure3C_n_tissues_eGenes_data.csv")
        
        bind_rows(
          counted2 %>% filter(n_tissues<=19) %>% mutate(n_tissues=as.character(n_tissues)),
          tribble(~n_tissues, ~n_egenes, ">=20", counted2 %>% filter(n_tissues>=20) %>% pull(n_egenes) %>% sum())
        ) %>% write_csv("manuscript/input_files/Figure3_Performance/figure3C_n_tissues_eGenes_data.csv")
        
      }
      
      # make input file single cell: DONE
      if(1){
        sc_eGenes = read_csv("../Input_Files/Figure3_Performance/eGenes.scMeta.all_contexts.residualized_exp_types.txt", col_types = cols(.default = col_character()))
        
        sc_eGenes_tissues_counted = sc_eGenes %>% filter(exp_type == "mean_norm_res_exp.specific") %>%
          mutate(exp_type = case_when(exp_type == "mean_norm_res_exp" ~ "tbt",
                                      exp_type == "mean_norm_res_exp.shared" ~ "fastgxe",
                                      exp_type == "mean_norm_res_exp.specific" ~ "fastgxe",
                                      T ~ "unknown?"))  %>% 
          filter(exp_type=="fastgxe") %>% 
          dplyr::select(gene, tissue) %>% distinct %>% 
          group_by(gene) %>% 
          summarize(n_tissues = n(), .groups="drop") %>% 
          dplyr::select(gene, n_tissues)
        
        counted2 = sc_eGenes_tissues_counted %>% 
          dplyr::select(n_tissues, gene) %>% 
          group_by(n_tissues) %>% 
          summarize(n_egenes = n(), .groups="drop") %>% 
          dplyr::select(n_tissues, n_egenes)
        
        counted2 %>% write_csv("../Input_Files/Figure3_Performance/figure3C_sc_n_tissues_eGenes_data.csv")
        
        #bind_rows(
        #  counted2 %>% filter(n_tissues<=4) %>% mutate(n_tissues=as.character(n_tissues)),
        #  tribble(~n_tissues, ~n_egenes, ">=5", counted2 %>% filter(n_tissues>=5) %>% pull(n_egenes) %>% sum())
        #) %>% write_csv("manuscript/input_files/Figure3_sc_Performance/figure3C_sc_n_tissues_eGenes_data.csv")
        
      }
      
      
      fig3C_hist_bulk_df = bind_rows(
        read_csv("../Input_Files/Figure3_Performance/figure3C_n_tissues_eGenes_data.csv", col_types = cols()) %>% 
          filter(n_tissues %in% seq(1,9)),
        tribble(~n_tissues, ~n_egenes,
                ">=10", read_csv("../Input_Files/Figure3_Performance/figure3C_n_tissues_eGenes_data.csv", col_types = cols()) %>% 
                  filter(!n_tissues %in% seq(1,9)) %>% pull(n_egenes) %>% sum())
      ) %>% 
        mutate(n_tissues = fct_relevel(n_tissues,c(as.character(seq(1,9)),">=10"))) %>% 
        #read_csv("manuscript/input_files/Figure3_GTEx_Performance/figure3C_n_tissues_data.csv", col_types = cols()) %>% 
        #mutate(n_tissues = fct_relevel(n_tissues,c(as.character(seq(1,19)),">=20"))) %>% 
        mutate(group = "Tissues") 
      
      
      fig3C_hist_sc_df = read_csv("../../manuscript/Input_Files/Figure3_Performance/figure3C_sc_n_tissues_eGenes_data.csv", col_types = cols()) %>%
        mutate(n_tissues = as.character(n_tissues)) %>% mutate(n_tissues = fct_relevel(n_tissues,c(as.character(seq(1,8))))) %>%
        # read_csv("manuscript/input_files/Figure3_GTEx_Performance/figure3C_n_tissues_data.csv", col_types = cols()) %>% 
        # mutate(n_tissues = fct_relevel(n_tissues,c(as.character(seq(1,19)),">=20"))) 
        mutate(group = "PBMC")
      
      
      ## facet wrapped figure 3C
      long_3C_df = rbind(fig3C_hist_bulk_df, fig3C_hist_sc_df)
      wrapped_3C_hist = long_3C_df %>% mutate(group=factor(group,levels=c("Tissues", "PBMC"))) %>%
        ggplot(aes(x=n_tissues, y=n_egenes, alpha = 0.5)) + 
        geom_bar(stat="identity") + 
        theme_bw() +
        theme(legend.position = "none",
              axis.text = element_text(size=23),
              axis.line = element_line(colour = "black"),
              axis.title = element_text(face = "bold", size=23,color="black"),
              panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.text = element_text(face = "bold", size = 25),
              panel.spacing = unit(5, "lines")) +
        #scale_y_continuous(labels = scientific_10)+
        xlab("Number of contexts") + ylab("Number of specific eGenes") + 
        facet_grid(~group, scales = "free_x", space = "free_x", labeller = label_value)
    }
    
    # plot: for single-tissue ts-eQTLs: which tissues are they found in?
    if(1){
      
      # make input file Bulk: 
      if(0){
        eGenes = read_csv("../Input_Files/Figure3_Performance/eGenes.v8.EUR.all_tissues.residualized_exp_types.txt")
        
        eGenes_single_tis_counted = eGenes %>% 
          filter(exp_type == "normalized_and_residualized_expression_heterogeneous" ) %>% 
          mutate(exp_type = case_when(exp_type == "normalized_and_residualized_expression" ~ "tbt",
                                      exp_type == "normalized_and_residualized_expression_homogeneous" ~ "fastgxe",
                                      exp_type == "normalized_and_residualized_expression_heterogeneous" ~ "fastgxe",
                                      T ~ "unknown?"))  %>% 
          dplyr::select(gene, tissue) %>% distinct %>% 
          group_by(gene) %>% 
          filter(n() == 1) %>% 
          ungroup %>% 
          dplyr::select(tissue, gene) %>% 
          group_by(tissue) %>% 
          summarize(n_egenes = n(), .groups="drop") %>% 
          dplyr::select(tissue, n_egenes)
        
        eGenes_single_tis_counted %>% 
          write_csv("../Input_Files/Figure3_Performance/figure3C_single_tissue_eGenes_data.csv") 
        
      }
      
      # make input file single cell:
      if(0){
        sc_eGenes = read_csv("/u/project/bballiu/bballiu/FastGxC/FastGxC_Manuscript/results/eGenes.scMeta.all_contexts.residualized_exp_types.txt", col_types = cols(.default = col_character()))
        
        eGenes_single_tis_counted = sc_eGenes %>% 
          filter(exp_type == "mean_norm_res_exp.specific" ) %>% 
          mutate(exp_type = case_when(exp_type == "mean_norm_res_exp" ~ "tbt",
                                      exp_type == "mean_norm_res_exp.shared" ~ "fastgxe",
                                      exp_type == "mean_norm_res_exp.specific" ~ "fastgxe",
                                      T ~ "unknown?"))  %>% 
          dplyr::select(gene, tissue) %>% distinct %>% 
          group_by(gene) %>% 
          filter(n() == 1) %>% 
          ungroup %>% 
          dplyr::select(tissue, gene) %>% 
          group_by(tissue) %>% 
          summarize(n_egenes = n(), .groups="drop") %>% 
          dplyr::select(tissue, n_egenes)
        
        eGenes_single_tis_counted %>% 
          write_csv("../Input_Files/Figure3_Performance/figure3C_sc_single_tissue_eGenes_data.csv") 
        
      }
      
      fig3C_bulk_df = read_csv("../Input_Files/Figure3_Performance/figure3C_single_tissue_eGenes_data.csv", col_types=cols()) %>% 
        filter(tissue!="AverageTissue") %>% 
        dplyr::select(tissue, n_egenes) %>% mutate(group = "Tissues") 
      
      fig3C_sc_df = read_csv("../Input_Files/Figure3_Performance/figure3C_sc_single_tissue_eGenes_data.csv", col_types=cols()) %>% 
        filter(tissue!="AverageContext") %>% 
        dplyr::select(tissue, n_egenes) %>% mutate(group = "PBMC") 
      
      ##### fig3C combined with fig3B
      ### BULK:
      ## get number of egenes that are not single tissue 
      non_single_tiss_bulk = fastgxe_eGenes_counted[-50,"n_egenes"] - fig3C_bulk_df$n_egenes
      combined_3BC_bulk_df = bind_rows(fastgxe_eGenes_counted[-50,] %>% dplyr::select(c("tissue", "n_egenes", "group")) %>% 
                                         mutate(component = "all_specific", n_egenes = non_single_tiss_bulk$n_egenes),
                                       fig3C_bulk_df %>% mutate(component = "tiss_specific"), 
                                       fastgxe_eGenes_counted[50,] %>% dplyr::select(c("tissue", "n_egenes", "group")) %>% 
                                         mutate(component = "tiss_specific"), 
                                       data.frame(tissue = "AverageTissue", n_egenes = 0, group = "Tissues", component = "all_specific"))
      
      combined_3BC_bulk = combined_3BC_bulk_df %>% ggplot(aes(x=fct_reorder(tissue,n_egenes,.desc = T), y=n_egenes, fill = tissue, alpha = component)) + 
        geom_bar(aes(fill = tissue), position="stack", stat="identity", width = 0.7) +
        gtex_abb_scale + gtex_fill_scale + 
        scale_alpha_manual(values = c(0.5, 1), labels = c("Context specific", "Single context specific"))+
        theme_bw() +
        theme(legend.position = c(0.6,0.7),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
              axis.title = element_text(face = "bold", size=25,color="black"),
              axis.text = element_text(size=21,color="black"),
              axis.line = element_line(colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              plot.title = element_text(hjust = 0.5,size=15),
              legend.title=element_blank(),
              legend.text=element_text(size=25),
              strip.text.y = element_text(size = 0),
              #strip.background = element_blank(),
              strip.text = element_text(face = "bold", size = 25),
              panel.spacing = unit(2, "lines")) + 
        facet_wrap(~ group, labeller = label_value)+
        scale_y_cut(breaks = c(7501,19999), which = c(1,2,3), scales = c(0.2,0,0.8), space = c(0.5))+
        scale_y_continuous( #labels = c(0, 2500, 5000, 7500, 20000, 22500), 
          breaks = c(0, 2500, 5000, 7500, 20000, 22500), limits = c(0,22500))+
        xlab("") + ylab("Number of eGenes")
      
      ### Single Cell:
      non_single_tiss_sc = fastgxe_sceGenes_counted[-1,"n_egenes"] - fig3C_sc_df$n_egenes
      combined_3BC_sc_df = bind_rows(fastgxe_sceGenes_counted[-1,] %>% dplyr::select(c("tissue", "n_egenes", "group")) %>% 
                                       mutate(component = "all_specific", n_egenes = non_single_tiss_sc$n_egenes),
                                     fig3C_sc_df %>% mutate(component = "tiss_specific"),
                                     fastgxe_sceGenes_counted[1,] %>% dplyr::select(c("tissue", "n_egenes", "group")) %>% 
                                       mutate(component = "tiss_specific"),
                                     data.frame(tissue = "AverageContext", n_egenes = 0, group = "PBMC", component = "all_specific"))
      
      combined_3BC_sc = combined_3BC_sc_df %>% ggplot(aes(x=fct_reorder(tissue,n_egenes,.desc = T), y=n_egenes, fill = tissue, alpha = component)) + 
        geom_bar(aes(fill = tissue), position="stack", stat="identity", width = 0.7) +
        sc_abb_scale + sc_fill_scale + 
        scale_alpha_manual(values = c(0.5, 1), labels = c("Context specific", "Single context specific"))+
        theme_bw() +
        theme(legend.position= "right",
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
              axis.text.y = element_blank(),
              axis.title = element_text(size=25,color="black"),
              axis.text = element_text(size=21,color="black"),
              axis.line = element_line(colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              axis.line.y = element_blank(),
              axis.ticks.y = element_blank(),
              plot.title = element_text(hjust = 0.5,size=15),
              legend.title=element_blank(),
              legend.text=element_text(size=15),
              strip.text.y = element_text(size = 0),
              #strip.background = element_blank(),
              strip.text = element_text(face = "bold", size = 25),
              panel.spacing = unit(2, "lines")) + 
        #facet_grid_sc(rows = vars(is_avg), scales = list(y = scales_y)) + 
        facet_wrap(~group,labeller = label_value)+
        #scale_y_continuous(breaks = c(0, 2500, 5000, 7500), limits = c(0,7500)) +
        scale_y_cut(breaks = c(7501,19999), which = c(1,2,3), scales = c(0.2,0,0.8), space = c(0.5))+
        scale_y_continuous(
          breaks = c(0, 2500, 5000, 7500, 20000, 22500), limits = c(0,22500))+
        xlab("") + ylab("")
      
    } #1200x700
    
    joint_3BC = aplot::plot_list(combined_3BC_bulk + theme(legend.position = "none"), combined_3BC_sc + theme(legend.position = "none"), widths = c(0.7, 0.3)) 
    combined_BC_legend = get_legend(combined_3BC_bulk + theme(legend.position = "top"))
    
    Figure3=plot_grid(
      ggdraw() +
        draw_plot(plot = fig3A, x = 0.01, y = 0.53, width = 0.215, height = 0.45) +
        draw_plot(plot = wrapped_3C_hist, x = 0.23, y = 0.53, width = 0.75, height = 0.45) +
        draw_plot(plot = combined_BC_legend, x = 0, y = 0)+
        draw_plot(plot = joint_3BC, x = 0.01, y = 0.0, width = 0.98, height = 0.49) +
        draw_plot_label(label = c('A','B','C'),  x = c(0.0,0.22,0.0), y = c(0.99,0.99,0.5), size = 35),
      rel_heights = c(0.05, 1), 
      ncol = 1)
    
    #print(Figure3)
    ggsave(plot = Figure3,
           filename = 'Fig03_sc_bulk_eGenes_results.pdf',
           dpi = 300,
           height = 20, width = 22, units="in")
    
  }
  
  #### Figure S15: sample size with number of eGenes
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(1){
    shared_only_gtex = 6749
    shared_only_sc = 1644
    fastgxe_sc_tissue_ncells = read_csv("../Input_Files/Figure3_Performance/CLUES_OneK1K_NCells_by_Context.csv")
    eGenes_samp_size_bulk = all_egenes %>% mutate(cohort = case_when(cohort == "GTEx" ~ "Tissues",
                                                                     cohort == "SC" ~ "PBMC",
                                                                     TRUE ~ cohort)) %>% filter(exp_type %in% c("Specific", "CxC")) %>% group_by(cohort, exp_type, tissue) %>% 
      mutate(exp_type = case_when(exp_type %in% c("Specific", "Shared") ~ "FastGxC",
                                  TRUE ~ exp_type)) %>%
      mutate(n = n()) %>% dplyr::select(c("exp_type", "tissue", "cohort", "n")) %>% distinct %>% filter(cohort == "Tissues") %>%
      merge(fastgxe_tissue_nsamples, by = "tissue") #%>% mutate(n = case_when(exp_type == "FastGxC" ~ n+shared_only_gtex,
    #                         TRUE ~ n))
    eGenes_samp_size_bulk = bind_rows(eGenes_samp_size_bulk %>% mutate(facet_lab = 1), fig3C_bulk_df %>% mutate(exp_type = "FastGxC Single Context") %>% 
                                        rename(cohort = "group") %>% rename(n = n_egenes) %>%
                                        merge(fastgxe_tissue_nsamples, by = "tissue") %>% mutate(facet_lab = 2))
    
    
    
    eGenes_samp_size_sc = all_egenes %>% mutate(cohort = case_when(cohort == "GTEx" ~ "Tissues",
                                                                   cohort == "SC" ~ "PBMC",
                                                                   TRUE ~ cohort)) %>% filter(exp_type %in% c("Specific", "CxC")) %>% group_by(cohort, exp_type, tissue) %>% 
      mutate(exp_type = case_when(exp_type %in% c("Specific", "Shared") ~ "FastGxC",
                                  TRUE ~ exp_type)) %>%
      mutate(n = n()) %>% dplyr::select(c("exp_type", "tissue", "cohort", "n")) %>% distinct %>% filter(cohort == "PBMC") %>% 
      merge(fastgxe_sc_tissue_ncells, by.x = "tissue", by.y = "final_celltypes") #%>% mutate(n = case_when(exp_type == "FastGxC" ~ n+shared_only_sc,
    #                          TRUE ~ n))
    eGenes_samp_size_sc = bind_rows(eGenes_samp_size_sc %>% mutate(facet_lab = 1), fig3C_sc_df %>% mutate(exp_type = "FastGxC Single Context") %>% 
                                      rename(cohort = "group") %>% rename(n = n_egenes) %>%
                                      merge(fastgxe_sc_tissue_ncells, by.x = "tissue", by.y = "final_celltypes") %>% mutate(facet_lab = 2))
    
    
    samp_size_bulk = ggplot(eGenes_samp_size_bulk, aes(x = n_samples, y = n, color = exp_type, alpha = exp_type))+
      geom_point(size = 5) + scale_color_manual(values = c("CxC" = "#c87e7e","FastGxC" = "#56A3E9", 'FastGxC Single Context' = "#56A3E9")) +
      scale_alpha_manual(values = c(1,0.3, 1))+
      geom_smooth(method='lm', se=FALSE, linewidth = 3) + 
      guides(color = guide_legend("exp_type"),
             alpha = guide_legend("exp_type"))+
      theme_bw() +
      #scale_y_continuous(labels = scientific_10)+
      #scale_x_continuous(labels = scientific_10)+
      facet_grid(facet_lab~cohort, scales = "free")+
      theme(
        strip.text.x = element_text(size = 30, face = "bold"),
        strip.text.y = element_blank(),
        #legend.position=c(0.7,0.95),
        legend.box = 'horizontal',
        # legend.position=c(0.75,0.99),
        legend.position="top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.background = element_blank(),
        # panel.border = element_blank(),
        #axis.title.x = element_blank(),
        #axis.title.y=element_blank(),
        #axis.text.y=element_blank(),
        #axis.ticks.y=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(face = "bold", size=25,color="black"),
        axis.text = element_text(size=20,color="black"),
        plot.title = element_text(hjust = 0.5,size=15),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title=element_blank(),
        legend.text=element_text(size=25))+
      xlab("Number of samples") + ylab("Number of eGenes")+
      geom_abline(slope=1, intercept = 0)
    
    cell_num_sc = ggplot(eGenes_samp_size_sc, aes(x = log10(num_cells), y = n, color = exp_type, alpha = exp_type))+
      geom_point(size = 5) + scale_color_manual(values = c("CxC" = "#c87e7e","FastGxC" = "#56A3E9", 'FastGxC Single Context' = "#56A3E9")) +
      scale_alpha_manual(values = c(1,0.3, 1))+
      geom_smooth(method='lm', se=FALSE, linewidth = 3) + 
      #scale_x_continuous(breaks = seq(0,960000, by = 160000), labels = (seq(0,960000, by = 160000)))+
      #scale_y_continuous(labels = scientific_10)+
      theme_bw() +
      facet_grid(facet_lab~cohort, scales = "free")+
      theme(
        strip.text.x = element_text(size = 30, face = "bold"),
        strip.text.y = element_blank(),
        legend.box = 'horizontal',
        # legend.position=c(0.75,0.99),
        legend.position="top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.background = element_blank(),
        # panel.border = element_blank(),
        #axis.title.x = element_blank(),
        axis.title.y=element_blank(),
        #axis.text.y=element_blank(),
        #axis.ticks.y=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(face = "bold", size=25,color="black"),
        axis.text = element_text(size=20,color="black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5,size=15),
        legend.title=element_blank(),
        legend.text=element_text(size=25))+
      xlab("log10 Number of cells") + ylab("Number of eGenes")+
      geom_abline(slope=1, intercept = 0)
    
    
    legend = get_legend(samp_size_bulk + theme(legend.position = "top"))
    supp_fig5 = plot_grid(legend, plot_grid(
      samp_size_bulk + theme(legend.position = "none"), cell_num_sc + theme(legend.position = "none"),
      rel_widths = c(1.2, 1)
    ), ncol = 1, rel_heights = c(1,12))
    
    ggsave(plot = supp_fig5,
           filename = 'FigureS15_samp_size.jpg', 
           width = 12,
           height = 6)
    
    ################ get metrics for paper (correlations)
    ## fastgxc correlations 
    tiss_samp_size = eGenes_samp_size_bulk %>% filter(exp_type == "FastGxC") %>% dplyr::select(n_samples)
    fastgxc_egenes = eGenes_samp_size_bulk %>% filter(exp_type == "FastGxC") %>% dplyr::select(n)
    fastgxc_singTiss_egenes = eGenes_samp_size_bulk %>% filter(exp_type == "FastGxC Single Context") %>% dplyr::select(n)
    cxc_egenes = eGenes_samp_size_bulk %>% filter(exp_type == "CxC") %>% dplyr::select(n)
    fastgxc_bulk_cor = cor(tiss_samp_size$n_samples, fastgxc_egenes$n, method = "spearman")
    cxc_bulk_cor = cor(tiss_samp_size$n_samples, cxc_egenes$n, method = "spearman")
    singTiss_bulk_cor = cor(tiss_samp_size$n_samples, fastgxc_singTiss_egenes$n, method = "spearman")
    
    num_cells = eGenes_samp_size_sc %>% filter(exp_type == "FastGxC") %>% dplyr::select(num_cells)
    fastgxc_egenes = eGenes_samp_size_sc %>% filter(exp_type == "FastGxC") %>% dplyr::select(n)
    fastgxc_singTiss_egenes = eGenes_samp_size_sc %>% filter(exp_type == "FastGxC Single Context") %>% dplyr::select(n)
    cxc_egenes = eGenes_samp_size_sc %>% filter(exp_type == "CxC") %>% dplyr::select(n)
    fastgxc_sc_cor = cor(num_cells$num_cells, fastgxc_egenes$n, method = "spearman")
    cxc_sc_cor = cor(num_cells$num_cells, cxc_egenes$n, method = "spearman")
    singTiss_sc_cor = cor(num_cells$num_cells, fastgxc_singTiss_egenes$n, method = "spearman")
    
    ### test for difference in correlation for paper ####
  }
}

#### Figure S14: FastGxC and CxC comparison
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(1){
  ### single cell
  egenes_sc = read_csv("../Input_Files/Figure3_Performance/eGenes.scMeta.all_contexts.residualized_exp_types.txt") %>%
    mutate(exp_type = case_when(exp_type == "mean_norm_res_exp" ~ "CxC",
                                exp_type == "mean_norm_res_exp.specific" ~ "FastGxC",
                                exp_type == "mean_norm_res_exp.shared" ~ "FastGxC",
                                TRUE ~ exp_type)) %>% mutate(cohort = "SC")
  egenes_counted = egenes_sc %>% group_by(gene) %>% summarize(n_exp_types = n_distinct(exp_type), 
                                                              group = case_when(n_distinct(exp_type) == 2 ~ "both",
                                                                                n_distinct(exp_type) == 1 & unique(exp_type) == "FastGxC" ~ "FastGxC only",
                                                                                n_distinct(exp_type) == 1 & unique(exp_type) == "CxC" ~ "CxC only",
                                                                                TRUE ~ NA_character_)) %>% distinct %>% 
    filter(group %in% c("FastGxC only", "CxC only")) %>% group_by(group) %>% mutate(value = n()) %>% dplyr::select(c("group", "value")) %>% distinct %>% mutate(cohort = "eGenes")
  
  ## eAssociations file too large to be saved on Github - only run once this takes a while
  if(0){
    eqtls_sc = read_csv("../../data/eAssociations.scMeta.all_contexts.residualized_exp_types.txt") %>%
      mutate(exp_type = case_when(exp_type == "mean_norm_res_exp" ~ "CxC",
                                  exp_type == "mean_norm_res_exp_homogeneous" ~ "FastGxC",
                                  exp_type == "mean_norm_res_exp.heterogeneous" ~ "FastGxC",
                                  TRUE ~ exp_type)) %>% mutate(cohort = "SC")
    
    eqtls_counted = eqtls_sc %>% group_by(snp, gene) %>% summarize(n_exp_types = n_distinct(exp_type), 
                                                                   group = case_when(n_distinct(exp_type) == 2 ~ "both",
                                                                                     n_distinct(exp_type) == 1 & unique(exp_type) == "FastGxC" ~ "FastGxC only",
                                                                                     n_distinct(exp_type) == 1 & unique(exp_type) == "CxC" ~ "CxC only",
                                                                                     TRUE ~ NA_character_)) %>% distinct %>% 
      filter(group %in% c("FastGxC only", "CxC only")) %>% group_by(group) %>% mutate(value = n()) %>% dplyr::select(c("group", "value")) %>% distinct %>% mutate(cohort = "eQTLs")
    
    esnps_counted = eqtls_sc %>% group_by(snp) %>% summarize(n_exp_types = n_distinct(exp_type), 
                                                             group = case_when(n_distinct(exp_type) == 2 ~ "both",
                                                                               n_distinct(exp_type) == 1 & unique(exp_type) == "FastGxC" ~ "FastGxC only",
                                                                               n_distinct(exp_type) == 1 & unique(exp_type) == "CxC" ~ "CxC only",
                                                                               TRUE ~ NA_character_)) %>% distinct %>% 
      filter(group %in% c("FastGxC only", "CxC only")) %>% group_by(group) %>% mutate(value = n()) %>% dplyr::select(c("group", "value")) %>% distinct %>% mutate(cohort = "eSNPs") %>%
      mutate(group = factor(group, levels = c("FastGxC only", "CxC only")))
    fwrite(eqtls_counted, file = "../Input_Files/Figure3_Performance/eqtls_sc_counted.csv", sep = ",")
    fwrite(esnps_counted, file = "../Input_Files/Figure3_Performance/esnps_sc_counted.csv", sep = ",")
  }
  eqtls_counted = fread("../Input_Files/Figure3_Performance/eqtls_sc_counted.csv", sep = ",", data.table = F)
  esnps_counted = fread("../Input_Files/Figure3_Performance/esnps_sc_counted.csv", sep = ",", data.table = F)
  
  
  #### bulk
  
  egenes_bulk = read_csv("../Input_Files/Figure3_Performance/eGenes.v8.EUR.all_tissues.residualized_exp_types.txt") %>%
    mutate(exp_type = case_when(exp_type == "normalized_and_residualized_expression" ~ "CxC",
                                exp_type == "normalized_and_residualized_expression_heterogeneous" ~ "FastGxC",
                                exp_type == "normalized_and_residualized_expression_homogeneous" ~ "FastGxC",
                                TRUE ~ exp_type)) %>% mutate(cohort = "SC")
  egenes_counted_bulk = egenes_bulk %>% group_by(gene) %>% summarize(n_exp_types = n_distinct(exp_type), 
                                                                     group = case_when(n_distinct(exp_type) == 2 ~ "both",
                                                                                       n_distinct(exp_type) == 1 & unique(exp_type) == "FastGxC" ~ "FastGxC only",
                                                                                       n_distinct(exp_type) == 1 & unique(exp_type) == "CxC" ~ "CxC only",
                                                                                       TRUE ~ NA_character_)) %>% distinct %>% 
    filter(group %in% c("FastGxC only", "CxC only")) %>% group_by(group) %>% mutate(value = n()) %>% dplyr::select(c("group", "value")) %>% distinct %>% mutate(cohort = "eGenes")
  
  ## eAssociations file too large to be saved on Github - only run once this takes a while
  if(0){
    eqtls_bulk = read_csv("../../data/eAssociations.v8.EUR.all_tissues.residualized_exp_types.all_stats.txt") %>%
      mutate(exp_type = case_when(exp_type == "normalized_and_residualized_expression" ~ "CxC",
                                  exp_type == "normalized_and_residualized_expression_homogeneous" ~ "FastGxC",
                                  exp_type == "normalized_and_residualized_expression_heterogeneous" ~ "FastGxC",
                                  TRUE ~ exp_type)) %>% mutate(cohort = "GTEx")
    
    eqtls_counted_bulk = eqtls_bulk %>% group_by(snp, gene_ensembl) %>% summarize(n_exp_types = n_distinct(exp_type), 
                                                                                  group = case_when(n_distinct(exp_type) == 2 ~ "both",
                                                                                                    n_distinct(exp_type) == 1 & unique(exp_type) == "FastGxC" ~ "FastGxC only",
                                                                                                    n_distinct(exp_type) == 1 & unique(exp_type) == "CxC" ~ "CxC only",
                                                                                                    TRUE ~ NA_character_)) %>% distinct %>% 
      filter(group %in% c("FastGxC only", "CxC only")) %>% group_by(group) %>% mutate(value = n()) %>% dplyr::select(c("group", "value")) %>% distinct %>% mutate(cohort = "eQTLs")
    
    esnps_counted_bulk = eqtls_bulk %>% group_by(snp) %>% summarize(n_exp_types = n_distinct(exp_type), 
                                                                    group = case_when(n_distinct(exp_type) == 2 ~ "both",
                                                                                      n_distinct(exp_type) == 1 & unique(exp_type) == "FastGxC" ~ "FastGxC only",
                                                                                      n_distinct(exp_type) == 1 & unique(exp_type) == "CxC" ~ "CxC only",
                                                                                      TRUE ~ NA_character_)) %>% distinct %>% 
      filter(group %in% c("FastGxC only", "CxC only")) %>% group_by(group) %>% mutate(value = n()) %>% dplyr::select(c("group", "value")) %>% distinct %>% mutate(cohort = "eSNPs") %>%
      mutate(group = factor(group, levels = c("FastGxC only", "CxC only")))
    fwrite(eqtls_counted_bulk, file = "../Input_Files/Figure3_Performance/eqtls_bulk_counted.csv", sep = ",")
    fwrite(esnps_counted_bulk, file = "../Input_Files/Figure3_Performance/esnps_bulk_counted.csv", sep = ",")
  }
  eqtls_counted_bulk = fread("../Input_Files/Figure3_Performance/eqtls_sc_counted.csv", sep = ",", data.table = F)
  esnps_counted_bulk = fread("../Input_Files/Figure3_Performance/esnps_sc_counted.csv", sep = ",", data.table = F)
  
  
  
  supp4_egenes = bind_rows(bind_rows(egenes_counted, eqtls_counted, esnps_counted) %>% mutate(study = "PBMC"), 
                           bind_rows(egenes_counted_bulk, eqtls_counted_bulk, esnps_counted_bulk) %>% mutate(study = "Tissues")) %>%
    filter(cohort == "eGenes") %>% mutate(study = factor(study, levels = c("Tissues", "PBMC"))) %>% ggplot(aes(fill=group, y=value, x=cohort))+ 
    geom_bar(position="dodge", stat="identity") +
    scale_fill_manual(name="",values = c("CxC only"="#c87e7e","FastGxC only"="#56A3E9")) +
    facet_grid(rows = vars(study), scale = "free", labeller = label_value)+
    theme_bw() + 
    theme(legend.position="none",
          axis.title = element_text(face = "bold", size=25,color="black"),
          axis.text = element_text(size=20,color="black"),
          axis.text.x = element_text(face = "bold"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(),
          plot.title = element_text(hjust = 0.5,size=15),
          legend.title=element_blank(),
          legend.text=element_text(size=18),
          strip.text = element_blank(),
          #strip.background = element_blank(),
          panel.spacing = unit(2, "lines")) + 
    xlab("") + ylab("Number of discoveries")
  
  supp4_eqtls = bind_rows(bind_rows(egenes_counted, eqtls_counted, esnps_counted) %>% mutate(study = "PBMC"), 
                          bind_rows(egenes_counted_bulk, eqtls_counted_bulk, esnps_counted_bulk) %>% mutate(study = "Tissues")) %>%
    filter(cohort == "eQTLs") %>% mutate(study = factor(study, levels = c("Tissues", "PBMC"))) %>% ggplot(aes(fill=group, y=value, x=cohort))+ 
    geom_bar(position="dodge", stat="identity") +
    scale_fill_manual(name="",values = c("CxC only"="#c87e7e","FastGxC only"="#56A3E9")) +
    facet_grid(rows = vars(study), scale = "free", labeller = label_value)+
    theme_bw() + 
    theme(legend.position="top",
          axis.title = element_text(face = "bold", size=25,color="black"),
          axis.text = element_text(size=20,color="black"),
          axis.text.x = element_text(face = "bold"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(),
          plot.title = element_text(hjust = 0.5,size=15),
          legend.title=element_blank(),
          legend.text=element_text(size=18),
          strip.text = element_blank(),
          #strip.background = element_blank(),
          panel.spacing = unit(2, "lines")) + 
    xlab("") + ylab("")
  
  supp4_esnps = bind_rows(bind_rows(egenes_counted, eqtls_counted, esnps_counted) %>% mutate(study = "PBMC"), 
                          bind_rows(egenes_counted_bulk, eqtls_counted_bulk, esnps_counted_bulk) %>% mutate(study = "Tissues")) %>%
    filter(cohort == "eSNPs") %>% mutate(study = factor(study, levels = c("Tissues", "PBMC"))) %>% ggplot(aes(fill=group, y=value, x=cohort))+ 
    geom_bar(position="dodge", stat="identity") +
    scale_fill_manual(name="",values = c("CxC only"="#c87e7e","FastGxC only"="#56A3E9")) +
    facet_grid(rows = vars(study), scale = "free", labeller = label_value)+
    theme_bw() + 
    theme(legend.position="none",
          axis.title = element_text(face = "bold", size=25,color="black"),
          axis.text = element_text(size=20,color="black"),
          axis.text.x = element_text(face = "bold"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(),
          plot.title = element_text(hjust = 0.5,size=15),
          legend.title=element_blank(),
          legend.text=element_text(size=18),
          strip.text = element_text(face = "bold", size = 25),
          #strip.background = element_blank(),
          panel.spacing = unit(2, "lines")) + 
    xlab("") + ylab("")
  
  legend = get_legend(supp4_eqtls)
  final_supp4 = plot_grid(plot_grid(
    legend
  ),
  ggdraw()+
    draw_plot(plot = supp4_egenes, x = 0.01, y = 0, width = 0.3) +
    draw_plot(plot = supp4_eqtls + theme(legend.position = "none"), x = 0.32, y = 0, width = 0.3) +
    draw_plot(plot = supp4_esnps, x = 0.64, y = 0,  width = 0.3),
  ncol = 1, rel_heights = c(0.5,10)
  )
  
  ggsave(plot = final_supp4,
         filename = 'FigureS14_FastGxC_vs_CxC.pdf', 
         width = 12,
         height = 6)
  
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%% Figure 4 : Correlation Heatmaps
#%%%%%%%%%%%%%%% Main: effect Sizes
#%%%%%%%%%%%%%%% Supps: intra-individual + genetic
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(1){
  library(pheatmap)
  
  # Main figure 4: effect size correlation heatmap
  if(1){
    
    # Fig4: correlation heatmap of effect sizes across tissues
    
    # fastgxe bulk and single cell
    if(1){
      
      # f.matr <- as.matrix(read.table("manuscript/input_files/Figure4_Correlation_Heatmap/fastgxe_hetonly_cor_spearman.txt",header=TRUE,row.names=1,check.names = F))
      f.matr <- as.matrix(read.table("../Input_Files/Figure4_Correlation_Heatmap/fastgxe_hetonly_cor_pearson.txt",header=TRUE,row.names=1,check.names = F))
      
      # switch name to abbreviation
      colnames(f.matr) <- gtex_abbrev.tib$abbreviation[match(colnames(f.matr),gtex_abbrev.tib$tissue)]
      rownames(f.matr) <- gtex_abbrev.tib$abbreviation[match(rownames(f.matr),gtex_abbrev.tib$tissue)]
      
      annot = data.frame(tissue = gtex_abbrev.tib$abbreviation)
      rownames(annot) = gtex_abbrev.tib$abbreviation
      
      ann_colors_vec <- unname(gtex_colors.vec)
      names(ann_colors_vec) <- as.character(gtex_abbrev.tib$abbreviation) 
      ann_colors = list(tissue = ann_colors_vec)
      
      pheat_breaks = seq(-1,1,0.2)
      # pheat_breaks = seq(-1,1,0.01)
      
      n_gaps = 12 # pearson = 12 gaps is perfect
      # n_gaps = 14 # spearman =  gaps is perfect
      
      fig4_heat = pheatmap(f.matr,
                           breaks = pheat_breaks,
                           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(pheat_breaks)-1),
                           cutree_rows = n_gaps, cutree_cols = n_gaps, # adds gaps
                           annotation_col = annot,
                           annotation_colors = ann_colors, 
                           annotation_names_col = F,       # removes label for "tissue"
                           annotation_legend = FALSE,      # removes legend for tissue color
                           treeheight_row = 0              # removes dendrogram on one side
      )
      ggsave(plot = fig4_heat,
             filename = 'Figure4_Correlation_Heatmap/Figure4_corheatmap_fastgxe.png',
             dpi = 300,
             height = 10, width = 10, units="in")
      
      #### SINGLE CELL FASTGXC
      # f.matr <- as.matrix(read.table("manuscript/input_files/Figure4_Correlation_Heatmap/fastgxe_hetonly_cor_spearman_sc.txt",header=TRUE,row.names=1,check.names = F))
      f.matr <- as.matrix(read.table("../Input_Files/Figure4_Correlation_Heatmap/scMeta_fastgxe_hetonly_cor_pearson_sc.txt",header=TRUE,row.names=1,check.names = F))
      
      # switch name to abbreviation
      colnames(f.matr) <- sc_abbrev.tib$abbreviation[match(colnames(f.matr),sc_abbrev.tib$tissue)]
      rownames(f.matr) <- sc_abbrev.tib$abbreviation[match(rownames(f.matr),sc_abbrev.tib$tissue)]
      
      annot = data.frame(tissue = sc_abbrev.tib$abbreviation)
      rownames(annot) = sc_abbrev.tib$abbreviation
      
      ann_colors_vec <- unname(sc_colors.vec)
      names(ann_colors_vec) <- as.character(sc_abbrev.tib$abbreviation) 
      ann_colors = list(tissue = ann_colors_vec)
      
      pheat_breaks = seq(-1,1,0.2)
      # pheat_breaks = seq(-1,1,0.01)
      
      n_gaps = 4 # pearson = 12 gaps is perfect
      # n_gaps = 14 # spearman =  gaps is perfect
      
      #callback = function(hc, mat){
      #  sv = svd((mat))$u[,2]
      #  dend = reorder(as.dendrogram(hc), wts = sv^2)
      #  as.hclust(dend)
      #}
      
      fig4_heat = pheatmap(f.matr,
                           breaks = pheat_breaks,
                           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(pheat_breaks)-1),
                           cutree_rows = n_gaps, cutree_cols = n_gaps, # adds gaps
                           annotation_col = annot,
                           annotation_colors = ann_colors, 
                           annotation_names_col = F,       # removes label for "tissue"
                           annotation_legend = FALSE,      # removes legend for tissue color
                           treeheight_row = 0,              # removes dendrogram on one side
                           #display_numbers = round(f.matr, 2),
                           #clustering_callback = callback,
                           fontsize = 20
      )
      
      col_dend <- fig4_heat[[2]]
      col_dend <- dendextend::rotate(col_dend, order = c("B", "NK", "CD4", "CD8", "PDC", "CDC", "NCMONO", "CMONO") )
      fig4_heat = pheatmap(f.matr,
                           breaks = pheat_breaks,
                           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(pheat_breaks)-1),
                           cutree_rows = n_gaps, cutree_cols = n_gaps, # adds gaps
                           annotation_col = annot,
                           annotation_colors = ann_colors, 
                           annotation_names_col = F,       # removes label for "tissue"
                           annotation_legend = FALSE,      # removes legend for tissue color
                           treeheight_row = 0,              # removes dendrogram on one side
                           #display_numbers = round(f.matr, 2),
                           clustering_callback = callback,
                           fontsize = 20,
                           cluster_cols=as.hclust(col_dend),
                           cluster_rows = as.hclust(col_dend)
      )
      
      ggsave(plot = fig4_heat,
             filename = 'Figure4_Correlation_Heatmap/Figure4_corheatmap_fastgxe.sc.png',
             dpi = 300,
             height = 10, width = 10, units="in")
    }
    
    # tbt
    if(1){
      
      library(pheatmap)
      
      # f.matr <- as.matrix(read.table("manuscript/input_files/Figure4_Correlation_Heatmap/fastgxe_hetonly_cor_spearman.txt",header=TRUE,row.names=1,check.names = F))
      f.matr <- as.matrix(read.table("..//Input_Files/Figure4_Correlation_Heatmap/tbt_cor_pearson.txt",header=TRUE,row.names=1,check.names = F))
      
      # switch name to abbreviation
      colnames(f.matr) <- gtex_abbrev.tib$abbreviation[match(colnames(f.matr),gtex_abbrev.tib$tissue)]
      rownames(f.matr) <- gtex_abbrev.tib$abbreviation[match(rownames(f.matr),gtex_abbrev.tib$tissue)]
      
      annot = data.frame(tissue = gtex_abbrev.tib$abbreviation)
      rownames(annot) = gtex_abbrev.tib$abbreviation
      
      ann_colors_vec <- unname(gtex_colors.vec)
      names(ann_colors_vec) <- as.character(gtex_abbrev.tib$abbreviation) 
      ann_colors = list(tissue = ann_colors_vec)
      
      pheat_breaks = seq(-1,1,0.2)
      # pheat_breaks = seq(-1,1,0.01)
      
      n_gaps = 12 # pearson = 12 gaps is perfect
      # n_gaps = 14 # spearman =  gaps is perfect
      
      fig4_heat = pheatmap(f.matr,
                           breaks = pheat_breaks,
                           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(pheat_breaks)-1),
                           cutree_rows = n_gaps, cutree_cols = n_gaps, # adds gaps
                           annotation_col = annot,
                           annotation_colors = ann_colors, 
                           annotation_names_col = F,       # removes label for "tissue"
                           annotation_legend = FALSE,      # removes legend for tissue color
                           treeheight_row = 0              # removes dendrogram on one side
      )
      
      
      ggsave(plot = fig4_heat,
             filename = 'Figure4_Correlation_Heatmap/Figure4_corheatmap_tbt.png',
             dpi = 300,
             height = 10, width = 10, units="in")
      
      
      # f.matr <- as.matrix(read.table("manuscript/input_files/Figure4_Correlation_Heatmap/fastgxe_hetonly_cor_spearman.txt",header=TRUE,row.names=1,check.names = F))
      f.matr <- as.matrix(read.table("../Input_Files/Figure4_Correlation_Heatmap/scMeta_tbt_cor_pearson_sc.txt",header=TRUE,row.names=1,check.names = F))
      
      # switch name to abbreviation
      colnames(f.matr) <- sc_abbrev.tib$abbreviation[match(colnames(f.matr),sc_abbrev.tib$tissue)]
      rownames(f.matr) <- sc_abbrev.tib$abbreviation[match(rownames(f.matr),sc_abbrev.tib$tissue)]
      
      annot = data.frame(tissue = sc_abbrev.tib$abbreviation)
      rownames(annot) = sc_abbrev.tib$abbreviation
      
      ann_colors_vec <- unname(sc_colors.vec)
      names(ann_colors_vec) <- as.character(sc_abbrev.tib$abbreviation) 
      ann_colors = list(tissue = ann_colors_vec)
      
      pheat_breaks = seq(-1,1,0.2)
      # pheat_breaks = seq(-1,1,0.01)
      
      n_gaps = 4 # pearson = 12 gaps is perfect
      # n_gaps = 14 # spearman =  gaps is perfect
      
      
      callback = function(hc, mat){
        sv = svd(t(mat))$v[,7]
        dend = reorder(as.dendrogram(hc), wts = sv)
        as.hclust(dend)
      }
      
      fig4_heat = pheatmap(f.matr,
                           breaks = pheat_breaks,
                           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(pheat_breaks)-1),
                           cutree_rows = n_gaps, cutree_cols = n_gaps, # adds gaps
                           annotation_col = annot,
                           annotation_colors = ann_colors, 
                           annotation_names_col = F,       # removes label for "tissue"
                           annotation_legend = FALSE,      # removes legend for tissue color
                           treeheight_row = 0,              # removes dendrogram on one side
                           #display_numbers = round(f.matr, 2),
                           clustering_callback = callback,
                           fontsize = 20
      )
      
      col_dend <- fig4_heat[[2]]
      col_dend <- dendextend::rotate(col_dend, order = c("B", "NK", "CD4", "CD8", "NCMONO", "CMONO", "PDC", "CDC") )
      fig4_heat = pheatmap(f.matr,
                           breaks = pheat_breaks,
                           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(pheat_breaks)-1),
                           cutree_rows = n_gaps, cutree_cols = n_gaps, # adds gaps
                           annotation_col = annot,
                           annotation_colors = ann_colors, 
                           annotation_names_col = F,       # removes label for "tissue"
                           annotation_legend = FALSE,      # removes legend for tissue color
                           treeheight_row = 0,              # removes dendrogram on one side
                           #display_numbers = round(f.matr, 2),
                           #clustering_callback = callback,
                           fontsize = 20,
                           cluster_cols=as.hclust(col_dend),
                           cluster_rows = as.hclust(col_dend)
      )
      
      ggsave(plot = fig4_heat,
             filename = 'Figure4_Correlation_Heatmap/Figure4_corheatmap_tbt.sc.png',
             dpi = 300,
             height = 10, width = 10, units="in")
    }
    
  }
  
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%% Figure 5 : Functional Characterization 
#%%%%%%%%%%%%%%% of GTEx eQTLs (Enrichment Analysis)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(1){
  
  # fig6A: show both sh-/ts-only + fastgxe/cbc only in one figure, as facets
  if(1){
    
    features_to_save = c("Promoter","Enhancer", "CTCF binding site",
                         "TF binding site","5' UTR","3' UTR","Upstream gene",
                         "Downstream gene")
    
    # fig_6A_left GTEx
    if(1){
      
      source("../scripts/geom_stripes.R")
      
      major_features_bulk = c("promoter","enhancer","promoter_flanking_region","CTCF_binding_site",
                              "TF_binding_site","5_prime_UTR_variant","3_prime_UTR_variant","upstream_gene_variant",
                              "downstream_gene_variant","missense_variant","synonymous_variant",
                              "intron_variant","non_coding_transcript_exon_variant")
      
      major_features_label_bulk = c("Promoter","Enhancer","Promoter-flanking","CTCF binding site",
                                    "TF binding site","5' UTR","3' UTR","Upstream gene",
                                    "Downstream gene","Missense","Synonymous",
                                    "Intron","NC transcript Exon")
      
      ag_data = read_csv("../Input_Files/Figure6_Enrichment/Enrichment.Tissue_Agnostic.SNPs_Matched_by_MAF.additional_SNP_sets.VEP_Annotations.FinalFisherResults.csv") %>% 
        mutate(set_type = case_when(set %in% c("FastGxE_all","TBT_all") ~ "compare all method results",
                                    set %in% c("FastGxE_only","FastGxE_TBT_intersect","TBT_only") ~ "compare by intersecting methods",
                                    set %in% c("HET_only","HOM_HET_intersect","HOM_only") ~ "compare within FastGxE",
                                    set %in% c("HET_all","HOM_all") ~ "FastGxE components"
        )) %>% mutate(odds_ratio = as.double(odds_ratio), conf_int.upper = as.double(conf_int.upper)) %>%
        filter(feature %in% major_features_bulk) %>% 
        mutate(feature = case_when(feature == "promoter" ~ "Promoter",
                                   feature == "enhancer" ~ "Enhancer",
                                   feature == "promoter_flanking_region" ~ "Promoter-flanking",
                                   feature == "CTCF_binding_site" ~ "CTCF binding site",
                                   feature == "TF_binding_site" ~ "TF binding site",
                                   feature == "5_prime_UTR_variant" ~ "5' UTR",
                                   feature == "3_prime_UTR_variant" ~ "3' UTR",
                                   feature == "upstream_gene_variant" ~ "Upstream gene",
                                   feature == "downstream_gene_variant" ~ "Downstream gene",
                                   feature == "missense_variant" ~ "Missense",
                                   feature == "synonymous_variant" ~ "Synonymous",
                                   feature == "intron_variant" ~ "Intron",
                                   feature == "non_coding_transcript_exon_variant" ~ "NC transcript Exon")) %>% 
        filter(set_type=="compare within FastGxE") %>% 
        filter(set %in% c("HET_only", "HOM_only"))  %>% 
        mutate(p.bh = p.adjust(p.value, method="BH")) %>%
        mutate(is_sign = case_when(p.bh <= 0.05 ~ "yes", T ~ "no")) 
      
    }
    
    ## Figure 6A single cell
    if(1){
      
      source("../scripts/geom_stripes.R")
      
      major_features_sc = c("promoter","enhancer","promoter_flanking_region","CTCF_binding_site",
                            "TF_binding_site","5_prime_UTR","3_prime_UTR","upstream_gene",
                            "downstream_gene","missense","synonymous",
                            "intron","non_coding_transcript_exon")
      
      major_features_label_sc = c("Promoter","Enhancer", "Promoter-flanking", "CTCF binding site",
                                  "TF binding site","5' UTR","3' UTR","Upstream gene",
                                  "Downstream gene","Missense","Synonymous",
                                  "Intron","NC transcript Exon")
      
      ag_data_sc_left = read_csv("../Input_Files/Figure6_Enrichment/Enrichment.Tissue_Agnostic.VEP.fisher_results_fdr_hom_het_only.sc.csv") %>% 
        filter(desc %in% major_features_sc) %>% 
        mutate(desc = case_when(desc == "promoter" ~ "Promoter",
                                desc == "enhancer" ~ "Enhancer",
                                desc == "promoter_flanking_region" ~ "Promoter-flanking",
                                desc == "CTCF_binding_site" ~ "CTCF binding site",
                                desc == "TF_binding_site" ~ "TF binding site",
                                desc == "5_prime_UTR" ~ "5' UTR",
                                desc == "3_prime_UTR" ~ "3' UTR",
                                desc == "upstream_gene" ~ "Upstream gene",
                                desc == "downstream_gene" ~ "Downstream gene",
                                desc == "missense" ~ "Missense",
                                desc == "synonymous" ~ "Synonymous",
                                desc == "intron" ~ "Intron",
                                desc == "non_coding_transcript_exon" ~ "NC transcript Exon")) %>% 
        filter(set %in% c("HET", "HOM"))  %>% 
        mutate(p.bh = p.adjust(p.value, method="BH")) %>%
        mutate(is_sign = case_when(p.bh <= 0.05 ~ "yes", T ~ "no")) 
      
    }
    
    fig_6A_left_combined = bind_rows(ag_data %>% dplyr::select(c("set", "feature", "odds_ratio", "conf_int.lower", "conf_int.upper", "p.bh", "is_sign")) %>%
                                       mutate(group = "Tissues") %>% mutate(feature = factor(feature, levels = rev(major_features_label_bulk))) %>% 
                                       mutate(set = case_when(set == "HET_only" ~ "HET",
                                                              set == "HOM_only" ~ "HOM",
                                                              TRUE ~ set)), 
                                     ag_data_sc_left %>% rename(feature = "desc") %>% dplyr::select(c("set", "feature", "odds_ratio", "conf_int.lower", "conf_int.upper", "p.bh", "is_sign")) %>%
                                       mutate(group = "PBMC") %>% mutate(feature = factor(feature, levels = rev(major_features_label_sc)))) %>%
      mutate(facet_lab = "FastGxC shared vs. specific") %>% filter(feature %in% features_to_save) %>% mutate(group = factor(group, levels = c("Tissues", "PBMC")))
    
    fig_6A_left_comb = fig_6A_left_combined %>% ggplot(aes(y=feature, x=odds_ratio, shape=set, fill = is_sign)) +
      facet_grid(group~facet_lab) + 
      geom_errorbarh(aes(xmax = conf_int.lower, xmin = conf_int.upper), size = .5) +
      geom_point(size=5, color = "#579ABE") +
      geom_vline(xintercept=1, linetype="dashed", color = "black") +
      scale_fill_manual(values=c("yes" = "#579ABE", "no" = "white"), 
                        labels = c("yes" = "FDR <=5%", "no" = "FDR >5%"), name = "") + 
      scale_shape_manual(values=c("HET" = 21, "HOM" = 24),
                         labels=c("HET" = "sp-eQTL", "HOM" = "sh-eQTL"), name = "") +
      guides(fill = guide_legend(nrow = 2, byrow = T, override.aes=list(shape=21))) +
      guides(shape = guide_legend(nrow = 2, byrow = T)) +
      theme_bw() +
      theme(
        legend.box = 'horizontal',
        
        # legend.position=c(0.75,0.75),
        legend.position="top",
        
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.background = element_blank(),
        # panel.border = element_blank(),
        strip.text.x = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size=20,color="black"),
        axis.text = element_text(size=17,color="black"),
        plot.title = element_text(hjust = 0.5,size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold", size = 20),
        strip.text.y = element_blank()) + 
      #scale_x_continuous(limits=c(0.7,1.4)) +
      # scale_y_continuous(limits=c(0.75,2)) +
      xlab("") + ylab("")
    
    
    # fig_6A_right Bulk
    if(1){
      
      major_features_bulk = c("promoter","enhancer","promoter_flanking_region","CTCF_binding_site",
                              "TF_binding_site","5_prime_UTR_variant","3_prime_UTR_variant","upstream_gene_variant",
                              "downstream_gene_variant","missense_variant","synonymous_variant",
                              "intron_variant","non_coding_transcript_exon_variant")
      
      major_features_label_bulk = c("Promoter","Enhancer","Promoter-flanking","CTCF binding site",
                                    "TF binding site","5' UTR","3' UTR","Upstream gene",
                                    "Downstream gene","Missense","Synonymous",
                                    "Intron","NC transcript Exon")
      
      # removed features where both methods are not significant
      # major_features = c("promoter","enhancer","promoter_flanking_region","CTCF_binding_site",
      #                    "3_prime_UTR_variant",
      #                    "downstream_gene_variant",
      #                    "intron_variant","non_coding_transcript_exon_variant")
      # 
      # major_features_label = c("Promoter","Enhancer","Promoter-flanking","CTCF binding site",
      #                          "3' UTR",
      #                          "Downstream gene",
      #                          "Intron","NC transcript")
      
      ag_data = read_csv("../Input_Files/Figure6_Enrichment/Enrichment.Tissue_Agnostic.SNPs_Matched_by_MAF.additional_SNP_sets.VEP_Annotations.FinalFisherResults.csv") %>% 
        mutate(set_type = case_when(set %in% c("FastGxE_all","TBT_all") ~ "compare all method results",
                                    set %in% c("FastGxE_only","FastGxE_TBT_intersect","TBT_only") ~ "compare by intersecting methods",
                                    set %in% c("HET_only","HOM_HET_intersect","HOM_only") ~ "compare within FastGxE",
                                    set %in% c("HET_all","HOM_all") ~ "FastGxE components"
        )) %>% mutate(odds_ratio = as.double(odds_ratio), conf_int.upper = as.double(conf_int.upper)) %>%
        filter(feature %in% major_features_bulk) %>% 
        mutate(feature = case_when(feature == "promoter" ~ "Promoter",
                                   feature == "enhancer" ~ "Enhancer",
                                   feature == "promoter_flanking_region" ~ "Promoter-flanking",
                                   feature == "CTCF_binding_site" ~ "CTCF binding site",
                                   feature == "TF_binding_site" ~ "TF binding site",
                                   feature == "5_prime_UTR_variant" ~ "5' UTR",
                                   feature == "3_prime_UTR_variant" ~ "3' UTR",
                                   feature == "upstream_gene_variant" ~ "Upstream gene",
                                   feature == "downstream_gene_variant" ~ "Downstream gene",
                                   feature == "missense_variant" ~ "Missense",
                                   feature == "synonymous_variant" ~ "Synonymous",
                                   feature == "intron_variant" ~ "Intron",
                                   feature == "non_coding_transcript_exon_variant" ~ "NC transcript Exon")) %>% 
        filter(set_type=="compare by intersecting methods") %>% 
        filter(set %in% c("FastGxE_only", "TBT_only"))  %>% 
        mutate(p.bh = p.adjust(p.value, method="BH")) %>%
        mutate(is_sign = case_when(p.bh <= 0.05 ~ "yes", T ~ "no")) 
      
      
      major_features_sc = c("promoter","enhancer","promoter_flanking","CTCF_binding_site",
                            "TF_binding_site","5_prime_UTR","3_prime_UTR","upstream_gene",
                            "downstream_gene","missense","synonymous",
                            "intron","non_coding_transcript_exon")
      
      ag_data_sc = read_csv("../Input_Files/Figure6_Enrichment/Enrichment.Tissue_Agnostic.VEP.fisher_results_fdr_hom_het_only.sc.csv") %>% 
        rename(feature = "desc") %>% filter(feature %in% major_features_sc) %>% 
        mutate(feature = case_when(feature == "promoter" ~ "Promoter",
                                   feature == "enhancer" ~ "Enhancer",
                                   feature == "promoter_flanking_region" ~ "Promoter-flanking",
                                   feature == "CTCF_binding_site" ~ "CTCF binding site",
                                   feature == "TF_binding_site" ~ "TF binding site",
                                   feature == "5_prime_UTR" ~ "5' UTR",
                                   feature == "3_prime_UTR" ~ "3' UTR",
                                   feature == "upstream_gene" ~ "Upstream gene",
                                   feature == "downstream_gene" ~ "Downstream gene",
                                   feature == "missense" ~ "Missense",
                                   feature == "synonymous" ~ "Synonymous",
                                   feature == "intron" ~ "Intron",
                                   feature == "non_coding_transcript_exon" ~ "NC transcript Exon")) %>% 
        filter(set %in% c("fastgxc_only", "TBT"))  %>% 
        mutate(p.bh = p.adjust(p.value, method="BH")) %>%
        mutate(is_sign = case_when(p.bh <= 0.05 ~ "yes", T ~ "no")) 
      
      
      fig_6A_right_combined = bind_rows(ag_data %>% dplyr::select(c("set", "feature", "odds_ratio", "conf_int.lower", "conf_int.upper", "p.bh", "is_sign")) %>%
                                          mutate(feature = factor(feature, levels = rev(major_features_label_bulk))) %>% 
                                          mutate(group = "Tissues") %>% mutate(set = case_when(set == "FastGxE_only" ~ "fastgxc_only",
                                                                                               set == "TBT_only" ~ "TBT",
                                                                                               TRUE ~ set)), 
                                        ag_data_sc %>% dplyr::select(c("set", "feature", "odds_ratio", "conf_int.lower", "conf_int.upper", "p.bh", "is_sign")) %>%
                                          mutate(feature = factor(feature, levels = rev(major_features_label_sc))) %>% 
                                          mutate(group = "PBMC")) %>% mutate(facet_lab = "FastGxC vs. CxC") %>% filter(feature %in% features_to_save) %>% 
        mutate(group = factor(group, levels = c("Tissues", "PBMC")))
      
      fig_6A_right_comb = fig_6A_right_combined %>% ggplot(aes(y=feature, x=odds_ratio, color=set, shape=is_sign)) +
        facet_grid(group~facet_lab) + 
        geom_errorbarh(aes(xmax = conf_int.lower, xmin = conf_int.upper), size = .5) + 
        geom_point(size=5) +
        geom_vline(xintercept=1, linetype="dashed", color = "black") + 
        
        scale_color_manual(values = c("fastgxc_only"="#579ABE",
                                      "TBT"="#c87e7e"),
                           labels = c("fastgxc_only"="FastGxC",
                                      "TBT"="CxC"),
                           name = "") +
        scale_shape_manual(values=c("yes" = 16, "no" = 1),
                           labels = c("yes" = "FDR <=0.05", "no" = "FDR >0.05"), name = "") +
        
        theme_bw() +
        theme(
          strip.text.x = element_text(size = 20),
          #legend.position=c(0.7,0.95),
          legend.box = 'horizontal',
          # legend.position=c(0.75,0.99),
          legend.position="top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # panel.background = element_blank(),
          # panel.border = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size=20,color="black"),
          axis.text = element_text(size=20,color="black"),
          plot.title = element_text(hjust = 0.5,size=15),
          legend.title=element_text(size=15),
          legend.text=element_text(size=15),
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(face = "bold", size = 20)
        ) +
        # guides(shape = guide_legend(order = 2),col = guide_legend(order = 1)) +
        xlab("") + ylab("") + 
        guides(colour = guide_legend(nrow = 2, byrow = T)) +
        guides(shape = guide_legend(nrow = 2, byrow = T))+
        scale_x_continuous(limits=c(0.7,1.4))
      
    }
    
    
    # combine left, right
    if(1){
      
      fig6A = plot_grid(
        plot_grid(
          fig_6A_left_comb, fig_6A_right_comb,
          rel_widths = c(1.35, 1)
        ),
        ggdraw() + 
          draw_label("Odds ratio", fontface = "bold", size = 20, x = 0.53, y = 0.95, hjust = 0) +
          theme(plot.margin = margin(0, 0, 0, 0)),
        
        ncol = 1,
        rel_heights = c(12,1)
      )
    }
    
  }
  
  # Fig6B: enrichment in ATAC
  if(1){
    
    atac_intersect_results = read_csv("../Input_Files/Figure6_Enrichment/final_atac_intersect_results.csv")
    encode_order=read_csv("../Input_Files/Figure6_Enrichment/Heatmap_ENCODE_TissueOrder.csv", col_types = cols()) %>% pull(ENCODE_order)
    gtex_order=read_csv("../Input_Files/Figure6_Enrichment/Heatmap_GTEx_TissueOrder.csv",  col_types = cols()) %>% pull(GTEx_order)
    
    atac_dat = atac_intersect_results %>% 
      filter(GTEx_tissue %in% gtex_order, ENCODE_tissue %in% encode_order) %>% 
      mutate(ENCODE_tissue = fct_relevel(ENCODE_tissue,encode_order)) %>% 
      mutate(GTEx_tissue = fct_relevel(GTEx_tissue,gtex_order))
    
    
  }
  
  # Fig6B: enrichment in ATAC single cell
  if(1){
    
    atac_intersect_results = read_csv("../Input_Files/Figure6_Enrichment/Enrichment.Tissue_Specific.ENCODE_ATAC_Intersect.fisher_results_fdr.sc.csv")
    atac_order=c("B", "T", "NK", "Mye", "Open")
    sc_order=c("B", "CD4", "CD8", "NK", "cMono", "ncMono", "pDC", "cDC")
    
    atac_dat_sc = atac_intersect_results %>% 
      filter(CLUES_tissue %in% sc_order, ATAC_tissue %in% atac_order) %>% 
      mutate(ATAC_tissue = fct_relevel(ATAC_tissue,atac_order)) %>% 
      mutate(CLUES_tissue = fct_relevel(CLUES_tissue,sc_order)) %>%
      filter(set %in% c("HET.single_tissue", "TBT.single_tissue")) %>% 
      mutate(is_matched = case_when(
        ATAC_tissue == "B" & CLUES_tissue == "B" ~ "yes",
        ATAC_tissue == "T" & CLUES_tissue %in% c("CD4", "CD8") ~ "yes",
        ATAC_tissue == "NK" & CLUES_tissue == "NK" ~"yes",
        ATAC_tissue == "Mye" & CLUES_tissue %in% c("ncMono", "cMono") ~ "yes",
        TRUE ~ "no"
      )) %>%
      mutate(is_sign = case_when(
        log2(odds_ratio) >= 0 & p.adjusted.BH <= 0.05 ~ 1,
        TRUE ~ 0
      ))
    
    gtex_tissues = unique(atac_dat$GTEx_tissue)
    gtex_encode = unique(atac_dat$ENCODE_tissue)
    sc_contexts = unique(atac_dat_sc$CLUES_tissue)
    sc_atac = unique(atac_dat_sc$ATAC_tissue)
    
    gtexOC_sc = expand.grid(gtex_encode, sc_contexts)
    names(gtexOC_sc) = c("ATAC_tissue", "tissue")
    gtexOC_sc$set = "HET.single_tissue"
    gtexOC_sc$odds_ratio = NA
    gtexOC_sc$is_sign = 0
    gtexOC_sc$is_matched = "no"
    gtexOC_sc_tbt = gtexOC_sc %>% mutate(set = "TBT.single_tissue")
    gtexOC_sc = rbind(gtexOC_sc, gtexOC_sc_tbt)
    scOC_gtex = expand.grid(sc_atac, gtex_tissues)
    names(scOC_gtex) = c("ATAC_tissue", "tissue")
    scOC_gtex$set = "HET.single_tissue"
    scOC_gtex$odds_ratio = NA
    scOC_gtex$is_sign = 0
    scOC_gtex$is_matched = "no"
    scOC_gtex_tbt = scOC_gtex %>% mutate(set = "TBT.single_tissue")
    scOC_gtex = rbind(scOC_gtex, scOC_gtex_tbt)
    final_na_df = rbind(gtexOC_sc, scOC_gtex)
    ## combine with real data
    
    gtex_atac_dat = atac_dat %>% dplyr::select(c("ENCODE_tissue", "GTEx_tissue", "set", "odds_ratio", "is_sign", "is_matched")) %>%
      mutate(set = case_when(set == "FastGxE ts-eQTL Variants" ~ "HET.single_tissue",
                             set == "TbT eQTL Variants" ~ "TBT.single_tissue",
                             TRUE ~ set)) %>% 
      rename(ATAC_tissue = "ENCODE_tissue", tissue = "GTEx_tissue")
    
    sc_atac_dat = atac_dat_sc %>% dplyr::select(c("ATAC_tissue", "CLUES_tissue", "set", "odds_ratio", "is_sign", "is_matched")) %>%
      rename(tissue = "CLUES_tissue")
    
    final_merged_df = rbind(gtex_atac_dat, final_na_df, sc_atac_dat)
    final_merged_df$ATAC_tissue = factor(final_merged_df$ATAC_tissue, levels = c(atac_order, encode_order))
    final_merged_df$tissue = factor(final_merged_df$tissue, levels = c(gtex_order, sc_order))
    facet_labs = c(HET.single_tissue = "FastGxC single-context sp-eQTLs", TBT.single_tissue="CxC single-context eQTLs")
    
    fig6B_combined = ggplot(final_merged_df, aes(x = tissue, y = ATAC_tissue)) + 
      facet_wrap(~set, labeller = labeller(set = facet_labs)) +
      geom_tile(aes(fill= log2(odds_ratio))) +
      scale_fill_continuous_divergingx(name="log2(OR)",palette = 'RdBu', mid = 0, rev=T, na.value = "#dddddd") + 
      # scale_fill_distiller(palette="RdBu", direction = -1, name = "Odds ratio (log2)", limits=c(-3.5,3.5)) +
      geom_tile(data = subset(final_merged_df, is_matched == "yes"), color = "gray40", size = 0.4, alpha=0, linetype="dotted") + # add "red" boxes
      geom_point(data = subset(final_merged_df, is_sign == 1), color = "gray40") + # add significance dots
      both_abb_scale + 
      xlab("eQTL") + ylab("Open Chromatin") +
      theme_bw() + 
      theme(legend.position="left",
            legend.title = element_text(size = 12),
            legend.text = element_text(size=12),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 14),
            axis.text.y = element_text(size=13),
            axis.title = element_text(face = "bold", size=23),
            strip.text = element_text(face = "bold", size = 20),
            strip.background = element_rect(fill="white"))+
      scale_size_area(max_size = 3) 
    
  }
  
  
  ### manuscript claims: mcnemar test
  
  mat = gtex_atac_dat %>%
    filter(is_matched == "yes") %>%
    dplyr::select(-is_matched, -odds_ratio) %>% 
    pivot_wider(names_from = set, values_from = is_sign) %>% 
    summarise(
      all_yes = sum(HET.single_tissue == 1 & TBT.single_tissue == 1, na.rm = TRUE),
      HET_no_TBT_yes = sum(HET.single_tissue == 0 & TBT.single_tissue == 1, na.rm = TRUE),
      HET_yes_TBT_no = sum(HET.single_tissue == 1 & TBT.single_tissue == 0, na.rm = TRUE),
      all_no = sum(HET.single_tissue == 0 & TBT.single_tissue == 0, na.rm = TRUE)
    ) %>% unlist %>% matrix(byrow = F,ncol = 2)
  
  
  format(mcnemar.test(x = mat)$p.val, scientific=T, digits=2)
  
  mat_sc = sc_atac_dat %>%
    filter(is_matched == "yes") %>%
    dplyr::select(-is_matched, -odds_ratio) %>% 
    pivot_wider(names_from = set, values_from = is_sign) %>% 
    summarise(
      all_yes = sum(HET.single_tissue == 1 & TBT.single_tissue == 1, na.rm = TRUE),
      HET_no_TBT_yes = sum(HET.single_tissue == 0 & TBT.single_tissue == 1, na.rm = TRUE),
      HET_yes_TBT_no = sum(HET.single_tissue == 1 & TBT.single_tissue == 0, na.rm = TRUE),
      all_no = sum(HET.single_tissue == 0 & TBT.single_tissue == 0, na.rm = TRUE)
    ) %>% unlist %>% matrix(byrow = F,ncol = 2)
  
  format(mcnemar.test(x = mat_sc)$p.val, scientific=T, digits=2)
  
  # Fig6: final, merge panels together
  if(1){
    
    ## no figure 6C and combined figure B
    Fig6_combined = plot_grid(
      fig6A,
      plot_grid(fig6B_combined,
                ncol = 1,
                rel_heights = c(3.25),
                labels = c("B"),
                label_size = 30),
      ncol = 1,
      rel_heights = c(1.5, 2.5),
      labels = c("A",""),
      label_size = 30
    )
    
    
    ggsave(plot = Fig6_combined,
           filename = 'Fig05_EnrichFuncAnnot_combined.pdf',
           dpi = 300,
           height = 22, width = 20,
           limitsize = F)
  }
  
}

### manuscript claims
## comparing two methods 

## median odds ratio of FastGxC vs CxC
final_merged_df %>% filter(is_matched == "yes") %>% dplyr::select(set, ATAC_tissue, tissue, odds_ratio) %>% 
  pivot_wider(names_from = set, values_from = odds_ratio) %>% 
  mutate(enrichment_ratio=exp(HET.single_tissue)/exp(TBT.single_tissue)) %>% 
  summarise(median(exp(HET.single_tissue)), median(enrichment_ratio))

numerator = final_merged_df %>% filter(is_matched == "yes" & is_sign == 1) %>% group_by(set) %>% summarise(n = n())
numerator$n[1]/numerator$n[2]

numerator = final_merged_df %>% filter(set == "HET.single_tissue" & is_matched == "yes") %>% group_by(tissue, set) %>% summarise(any(is_sign == 1)) %>% group_by(set) %>% summarize(n = sum(`any(is_sign == 1)`))
denominator = final_merged_df %>% filter(is_matched == "yes") %>% distinct(tissue)  %>% pull() %>% length()
numerator$n/denominator



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%% Figure 6 and and S16: GTEx eQTL examples
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(1){
  
  grabbed = read_csv("../Input_Files/Figure5_GTEx_eQTL_Examples/combined_grabbed_snp_gene_tissue_eqtl_stats.csv", col_types = cols()) %>% 
    mutate(is_shared = if_else(exp_type=="Sh-","shared","not_shared")) %>% 
    mutate(method = case_when(method == "TbT" ~ "CxC",
                              method == "FastGxE" ~ "FastGxC",
                              TRUE ~ method))  %>% 
    mutate(exp_type = case_when(exp_type == "TbT" ~ "CxC",
                                TRUE ~ exp_type)) %>% 
    mutate(is_sign = case_when(is_sign ~ "significant", !is_sign ~ "not significant"))
  
  pcs=4
  ## get legend = figure5_legend ----
  figure5_legend = get_legend(grabbed %>%
                                mutate(gene_snp = paste(gene, snp)) %>% 
                                filter(gene_snp %in% "SIGLEC14 rs872629") %>%
                                # mutate(gene_snp = fct_relevel(gene_snp, eqtls_to_plot)) %>% 
                                mutate(tissue = fct_reorder(tissue, desc(tissue))) %>% 
                                mutate(tissue = fct_relevel(tissue, "AverageTissue", after = Inf)) %>% 
                                mutate(method = fct_relevel(method, "CxC")) %>% 
                                mutate(is_shared = fct_relevel(is_shared, "shared")) %>%
                                ggplot() +
                                geom_point(aes(x = beta, y = tissue, color = tissue, shape = is_sign), size=pcs) +
                                geom_vline(xintercept=c(0), linetype="dotted") +
                                geom_errorbarh(aes(xmax = upper_int, xmin = lower_int, y = tissue), size = .5, height = .2, color = "gray50") +
                                gtex_colors_scale +
                                scale_shape_manual(values=c("significant" = 16, "not significant" = 1), 
                                                   labels=c("significant" = "< 0.05", "not significant" = "> 0.05"),
                                                   name="FDR p-value") +
                                theme_bw() +
                                theme(legend.spacing.x = unit(0.001, 'cm'),
                                      legend.title=element_text(size=20), 
                                      legend.text=element_text(size=19),
                                      legend.position = "top")+
                                xlab("eQTL effect size") + ylab("") +
                                scale_y_discrete(labels=gtex_abbrev.vec))
  
  
  
  ## main: plot eqtl examples and save ----
  fig5 = plot_grid(
    figure5_legend,
    plot_grid(plot_eqtl_example(plot_gene = "CBS", plot_snp = "rs234728", y_axis_tissue_lab = T),
              plot_eqtl_example(plot_gene = "SIGLEC14", plot_snp = "rs872629", y_axis_tissue_lab = F),
              plot_eqtl_example(plot_gene = "LDHC", plot_snp = "rs4757652", y_axis_tissue_lab = F),
              plot_eqtl_example(plot_gene = "MCPH1", plot_snp = "rs74607633", y_axis_tissue_lab = F),
              nrow=1,
              rel_widths = c(1, 0.8, 0.8),
              labels = c("A","B","C"),
              label_size = 25),
    ncol=1,
    rel_heights = c(1,20)
  )
  
  ggsave(plot = fig5,
         filename = 'Fig06_GTEx_eQTL_Examples_MCPH1.pdf',
         dpi = 300,
         height = 10, width = 16.5,
         limitsize = F)
  
  ## supps: additional examples + LDHC manhattan plots  ----
  
  ## A: additional examples  ----
  figsupp_examples = plot_grid(
    figure5_legend,
    plot_grid(plot_eqtl_example(plot_gene = "MUC20P1", plot_snp = "rs139637885", y_axis_tissue_lab = T),
              plot_eqtl_example(plot_gene = "GBP3", plot_snp = "rs7544740", y_axis_tissue_lab = F),
              plot_eqtl_example(plot_gene = "GSTT2", plot_snp = "rs369691", y_axis_tissue_lab = F),
              nrow=1,
              rel_widths = c(1, 0.8, 0.8),
              labels = c("B","C","D"),
              label_size = 25),
    ncol=1,
    rel_heights = c(1,20)
  )
  
  ## B: LDHC manhattan ----
  ldhc = read_csv("../Input_Files/Figure5_GTEx_eQTL_Examples/ldhc_manhattan.csv")
  
  ldhc_final = tissues2plot = c("ADPSBQ", "ADRNLG", "ARTTBL",  "BREAST", "BRNCHA", "CLNSGM", "CLNTRN", "ESPGEJ", "ESPMCS", "ESPMSL", "HRTAA", "HRTLV", "KDNCTX", "LCL", "LIVER", "LUNG", "MSCLSK",  "NERVET", "OVARY", "PNCREAS", "PRSTTE", "PTTARY", "SKINS", "SLVRYG", "SNTTRM", "SPLEEN", "WHLBLD", "TESTIS", "SHARED")
  
  labels4tissues=c(paste(tissues2plot[1:(length(tissues2plot)-1)],'CxC', sep = '-'),"TESTIS-FastGxC", "SHARED-FastGxC")
  
  ldhc_final = ldhc %>% mutate(method=factor(x = method, levels = c("FastGxE","CxC"), labels = c("FastGxC","CxC"))) %>%
    filter(tissue %in% tissues2plot) %>% 
    mutate(tissue = factor(x = paste(tissue,method, sep = '-'), levels = labels4tissues)) %>%
    ggplot(aes(x=SNP_pos,y=-log10(p.value),color=positive_beta)) + 
    facet_wrap(~tissue, scale="free_y") + 
    geom_point() +
    geom_vline(xintercept = c(18399690), linetype="dashed") +
    geom_hline(yintercept = c(-log10(1e-6)), linetype="dashed",color="gray") +
    scale_color_manual(name="Genetic effect",
                       values=c("Negative"="blue","Positive"="red")) +
    theme_bw() + 
    theme(legend.position = "none",
          strip.text = element_text(size = 20),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size=15),
          axis.title = element_text(size = 20)) + 
    xlab("SNP position") + 
    ylab(expression(-log[10](P))) + 
    scale_x_continuous(limits=c(18350000,18600000),
                       breaks=seq(18350000,18600000,100000)) + 
    theme(legend.position = "top",
          plot.title = element_text(vjust = -4))
  
  
  
  # Combined Figure ----  
  fig_supp_total = plot_grid(
    ldhc_final,
    figsupp_examples,
    ncol=1,
    rel_heights = c(1.5,1),
    labels=c("A",""),
    label_size = 25
  )
  
  
  ggsave(plot = fig_supp_total,
         filename = 'FigureS16_GTEx_eQTL_Examples_test.pdf',
         dpi = 300,
         height = 20, width = 18, limitsize = F)
  
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%% Figure 7 : GWAS
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(1){
  
  ## Data preprocess ----
  if(1){
    #  read enrichment results
    raw_gwas = read_csv("../Input_Files/Figure7_GWAS/Final_Merged.GWAS_enrichment_results.All_Tissues.SNPs_MAF_matched.with_additional_SNP_sets.FastGxE_TBT.by_mapped_and_parent_gwas_traits.csv", col_types=cols()) %>% 
      filter(ct3+ct4>=10) %>% 
      # add HOM to FastGxE and change f_in to "each_tissue"
      mutate(f_in = case_when(met == "HOM" ~ "each_tissue", TRUE ~ f_in)) %>% 
      mutate(met = case_when(met == "HET" ~ "FastGxE", met == "HOM" ~ "FastGxE", met == "TBT" ~ "TbT")) %>% 
      dplyr::select(met, f_in, cat, tis, tra, n_gwas_snps, OR, conf_int1, conf_int2, pval, ct1, ct2, ct3, ct4)
    
    #  perform heiachical multiple testing and get significant results 
    tree_sign_results = tribble(~met, ~f_in, ~cat, ~tis, ~tra)
    for (method in unique(raw_gwas$met)){
      for (found_in in unique(raw_gwas$f_in)){
        for (trait_category in unique(raw_gwas$cat)){
          
          tib=raw_gwas %>% filter(met == method, f_in == found_in, trait_category == cat) %>% dplyr::select(tis,tra,pval)
          
          # per tissue, get a simes p value
          tissue_simes = tribble(~tissue,~simes_p)
          for (t in unique(tib$tis))
            tissue_simes = bind_rows(tissue_simes, tribble(~tissue,~simes_p,
                                                           t,      simes.test(tib[tib$tis==t,"pval"])))
          # BH adjust 49x tissue p values
          tissue_simes = tissue_simes %>% mutate(simes_p_BH = p.adjust(simes_p,method="BH"))
          
          # signficant tissue: BH adjusted p value < .05
          thresh1=0.05
          tissue_simes_sign = tissue_simes[tissue_simes$simes_p_BH<thresh1,"tissue",drop=T]
          
          # for thresholding within tissue
          print(length(tissue_simes_sign))
          # print(nrow(tissue_simes))
          
          thresh2 = thresh1 * (length(tissue_simes_sign) / nrow(tissue_simes))
          for (sign_tissue in tissue_simes_sign){
            
            sign_tissue_tib = tib %>%
              # BH adjust p values across all traits
              filter(tis == sign_tissue) %>% dplyr::select(tra, pval) %>%
              mutate(tree_pval = p.adjust(pval, method = "BH")) %>%
              # filter BH adjusted p values by adjusted threshold 2
              filter(tree_pval <= thresh2) %>%
              # save results
              mutate(met = method, f_in = found_in, cat = trait_category, tis = sign_tissue) %>% 
              dplyr::select(met, f_in, cat, tis, tra)
            tree_sign_results = bind_rows(tree_sign_results, sign_tissue_tib)
          }
        } 
      }
    }
    
    tree_sign_results = tree_sign_results %>% mutate(tree_sign = "yes")
    
    # annotate gwas results as tree-adjusted significant
    all_gwas_enrich_sig_results = raw_gwas %>%  
      left_join(tree_sign_results, by = c("met","f_in","cat","tis","tra")) %>% 
      mutate(tree_sign = replace_na(tree_sign,"no")) %>% 
      filter(tree_sign == "yes") %>% 
      dplyr::select(met, f_in, cat, tis, tra, n_gwas_snps, OR, pval, conf_int1, conf_int2, ct1, ct2, ct3, ct4)
    
    all_49_tissues = all_gwas_enrich_sig_results %>% dplyr::select(tis) %>% distinct
    
    # remove average tissue 
    # all_gwas_enrich_sig_results = all_gwas_enrich_sig_results %>% filter(tis != "AverageTissue")
    
    #  grab, clean, and merge my old manual tissue-trait annotations
    relevant_tissues_annot = bind_rows(
      
      read_csv("../Input_Files/Figure7_GWAS/082020.known-tissue-trait-associations-table.lu_et_al.csv", col_types = cols()) %>%
        dplyr::select(MAPPED_TRAIT, Annotated_GTEx_Tissues) %>%
        arrange(MAPPED_TRAIT) %>%
        mutate(Annotated_GTEx_Tissues = gsub("[,]",", ",Annotated_GTEx_Tissues)) %>%
        rename(trait = MAPPED_TRAIT, relevant_tissues = Annotated_GTEx_Tissues),
      
      read_csv("../Input_Files/Figure7_GWAS/GWAS_cancer_traits_matched.csv", col_types = cols()) %>% 
        pivot_longer(-cancer_trait, names_to = "tissue_num", values_to = "tissue") %>% 
        filter(!is.na(tissue)) %>% 
        dplyr::select(cancer_trait, tissue) %>% 
        group_by(cancer_trait) %>% 
        group_modify(function(tib,key) tribble(~relevant_tissues, paste0(tib$tissue,collapse=", "))) %>% 
        ungroup %>% 
        rename(trait = cancer_trait) %>% 
        dplyr::select(trait, relevant_tissues),
      
      tribble(~trait, ~relevant_tissues,
              "body fat distribution", "Adipose_Visceral_Omentum, Adipose_Subcutaneous",
              "body fat percentage", "Adipose_Visceral_Omentum, Adipose_Subcutaneous",
              "schizophrenia", "Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Nucleus_accumbens_basal_ganglia") ## Added Schizophrenia Here
      
    ) %>% distinct %>% group_by(trait) %>% 
      summarise(relevant_tissues = paste(unique(unlist(strsplit(relevant_tissues, split = ", "))), collapse = ", ")) %>% 
      arrange(trait)
    
    # Uncomment this if you want to save the list of relevant tissues for each trait in a separate excel
     relevant_tissues_annot %>%  WriteXLS::WriteXLS("../SuppTables/TableS5_gwas_enrichment_results.xlsx")
    
    # annotate relevant tissue-trait pairs 
    relevant_tissue_trait_pairs=unlist(sapply(1:nrow(relevant_tissues_annot), function(i) paste(relevant_tissues_annot$trait[i],unlist(strsplit(relevant_tissues_annot$relevant_tissues[i], split = ', ')), sep = "-")))
    all_gwas_enrich_sig_results %<>% mutate(relevant_tissue_4_trait = paste(tra,tis, sep = '-') %in% relevant_tissue_trait_pairs)
    
    # True/False positive and True/False negative tissue - trait mapping
    cont_tables = all_gwas_enrich_sig_results %>% 
      filter(tra %in% unique(relevant_tissues_annot$trait)) %>% 
      filter(f_in == "each_tissue", cat=="MAPPED_TRAIT") %>% 
      dplyr::select(met, tis, tra, OR) %>% 
      filter(OR>1) %>%   #### only ORs > 1 = enriched!
      group_by(tra) %>% 
      group_modify(function(tib,key){
        
        # annotate list of 50 tissues with either the tissue is relevant for this trait or not
        rel_tissues = relevant_tissues_annot %>% filter(trait == key$tra) %>% pull(relevant_tissues)
        rel_tissues_tib = strsplit(gsub(" ", "", rel_tissues, fixed = TRUE), ",")[[1]] %>%
          as_tibble() %>% rename(tis = value) %>% mutate(is_relevant = T) %>% 
          dplyr::select(tis, is_relevant)
        all_49_tissues.rel = all_49_tissues %>% 
          left_join(rel_tissues_tib, by=c("tis")) %>% 
          mutate(is_relevant = case_when(is.na(is_relevant) ~ F, !is.na(is_relevant) ~ is_relevant)) %>% 
          dplyr::select(tis, is_relevant)
        
        tib_annotated_all = tib %>%
          group_by(met) %>%
          group_modify(function(tib2,key2){
            
            # classify each tissue as enriched/not for this method
            enriched_tissues_rank = tib2 %>% arrange(desc(OR)) %>% mutate(rank = as.character(1:n()))
            
            all_49_tissues.rel.enr = all_49_tissues.rel %>% 
              left_join(enriched_tissues_rank, by=c("tis")) %>% 
              mutate(is_enriched = !is.na(OR)) %>% 
              rename(enrich_OR = OR, enrich_rank = rank) %>% 
              dplyr::select(tis, is_enriched, enrich_OR, enrich_rank, is_relevant)
            
            # contingency table
            TP = nrow(all_49_tissues.rel.enr %>% filter(is_relevant,is_enriched))
            FP = nrow(all_49_tissues.rel.enr %>% filter(!is_relevant,is_enriched))
            FN = nrow(all_49_tissues.rel.enr %>% filter(is_relevant,!is_enriched))
            TN = nrow(all_49_tissues.rel.enr %>% filter(!is_relevant,!is_enriched))
            
            tribble(~TP, ~FP, ~FN, ~TN, TP, FP, FN, TN)
            
          })
        
        if(length(unique(tib_annotated_all$met))==2){
          return(tib_annotated_all)
        } else if(length(unique(tib_annotated_all$met))==1){
          if(unique(tib_annotated_all$met)=="FastGxE") method_to_make_NA = "TbT"
          if(unique(tib_annotated_all$met)=="TbT") method_to_make_NA = "FastGxE"
          bind_rows(
            tib_annotated_all,
            tribble(~met, ~TP, ~FP, ~FN, ~TN, 
                    method_to_make_NA, 
                    0, 0, 
                    nrow(rel_tissues_tib), 
                    nrow(all_49_tissues) - nrow(rel_tissues_tib))
          )
        }
        
        
      }) %>%
      ungroup %>% 
      mutate(met = case_when(met == "TbT" ~ "CxC", met == "FastGxE" ~ "FastGxE", T ~ "???"))
    
    # Precision, recall, and F1
    pre_rec = cont_tables %>% 
      mutate(precision = TP / (TP + FP)) %>% 
      mutate(recall = TP / (TP + FN)) %>% 
      mutate(f1 = 2 * precision * recall / (precision + recall)) %>% 
      filter(!is.nan(f1)) %>% 
      rename(Precision = precision, Recall = recall, F1 = f1) %>% 
      dplyr::select(met, tra, Precision, Recall, F1) %>% 
      pivot_longer(-c("met","tra"),names_to = "score_desc",values_to = "score") %>% 
      mutate(score_desc = fct_relevel(score_desc,"Precision","Recall")) %>% 
      mutate(met = if_else(met=="FastGxE","FastGxC",met))
    
    # Rank of all enriched tissues
    ranks = all_gwas_enrich_sig_results %>% 
      filter(tra %in% unique(relevant_tissues_annot$trait)) %>% 
      filter(f_in == "each_tissue", cat=="MAPPED_TRAIT")  %>% 
      dplyr::select(met, tra, tis, relevant_tissue_4_trait, OR, pval) %>%
      mutate(minusLog10P=-log10(pval)) %>% arrange(met,tra,relevant_tissue_4_trait) %>% 
      group_by(met, tra)  %>% 
      mutate(rank_or=order(OR, decreasing = T), rank_p=order(minusLog10P, decreasing = T)) %>% 
      ungroup %>% 
      mutate(met = case_when(met == "TbT" ~ "CxC", met == "FastGxE" ~ "FastGxE", T ~ "???"))
    
    #### add method to split shared and specific
    cont_tables_shared <- cont_tables %>%
      mutate(met = case_when(
        tra == "cancer" & met == "FastGxE" ~ "FastGxC Shared",
        TRUE ~ met  # retain original value otherwise
      ))
    
    # Precision, recall, and F1
    pre_rec_shared = cont_tables_shared %>% 
      mutate(precision = TP / (TP + FP)) %>% 
      mutate(recall = TP / (TP + FN)) %>% 
      mutate(f1 = 2 * precision * recall / (precision + recall)) %>% 
      filter(!is.nan(f1)) %>% 
      rename(Precision = precision, Recall = recall, F1 = f1) %>% 
      dplyr::select(met, tra, Precision, Recall, F1) %>% 
      pivot_longer(-c("met","tra"),names_to = "score_desc",values_to = "score") %>% 
      mutate(score_desc = fct_relevel(score_desc,"Precision","Recall")) %>% 
      mutate(met = if_else(met=="FastGxE","FastGxC Specific",met))
    
    # Rank of all enriched tissues
    ranks_shared = all_gwas_enrich_sig_results %>% 
      filter(tra %in% unique(relevant_tissues_annot$trait)) %>% 
      filter(f_in == "each_tissue", cat=="MAPPED_TRAIT")  %>% 
      dplyr::select(met, tra, tis, relevant_tissue_4_trait, OR, pval) %>%
      mutate(minusLog10P=-log10(pval)) %>% arrange(met,tra,relevant_tissue_4_trait) %>% 
      group_by(met, tra)  %>% 
      mutate(rank_or=order(OR, decreasing = T), rank_p=order(minusLog10P, decreasing = T)) %>% 
      ungroup %>% 
      mutate(met = case_when(met == "TbT" ~ "CxC", met == "FastGxE" ~ "FastGxE", T ~ "???"))
    
  }
  
  #### code to plot shared as it's own line
  FigA=ggplot() + 
    scale_fill_manual(name="",values = c("CxC"="#c87e7e","FastGxC Specific"="#56A3E9")) +
    theme_bw() + 
    theme(legend.position = "none",
          axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
          #axis.ticks.x = element_blank(),
          legend.text=element_text(size=20),
          axis.line = element_line(colour = "black"),
          strip.text.x = element_text(size=20),
          axis.text = element_text(size=15, color = 'black'),
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size=12)) 
  
  FigA1 = FigA + geom_boxplot(data = pre_rec_shared %>% filter(score_desc=="Precision"), mapping = aes(x = met, y = 100*score, fill= met)) + ggtitle("") +
    theme(plot.title = element_text(hjust=0.5, face = "bold", size = 12), axis.title = element_text(size = 10, face = "bold")) + ylab("Precision (%)")
  
  FigA2 = FigA + geom_boxplot(data = pre_rec_shared %>% filter(score_desc=="Recall"), mapping = aes(x = met, y = 100*score, fill= met)) + ylab("Recall (%)") + ggtitle("") + 
    theme(plot.title = element_text(hjust=0.5, face = "bold", size = 12), axis.title = element_text(size = 10, face = "bold"))
  dat <- ggplot_build(FigA2)$data[[1]] %>% mutate(fill=factor(x = fill, levels = c("#c87e7e","#56A3E9"), labels = c("CxC","FastGxC Specific")))
  FigA2 = FigA2 + geom_segment(data=dat, aes(x=xmin, xend=xmax,  y=middle, yend=middle, color=fill), size=4) + 
    scale_color_manual(name="",values = c("CxC"="#c87e7e","FastGxC Specific"="#56A3E9"))
  
  fig7A_shared=gridExtra::grid.arrange(FigA1,FigA2, nrow=1)
  ggsave(plot = fig7A_shared,
         filename = 'reviews/Figure07A_shared_separate_precision_recall.jpg', 
         width = 5,
         height = 4)
  
  
  # Figure A: precision, recall, rank relevant tissue, nr of enriched tissues  ----
  FigA=ggplot() + 
    scale_fill_manual(name="",values = c("CxC"="#c87e7e","FastGxC"="#56A3E9")) +
    theme_bw() + 
    theme(legend.position = "none",
          axis.text.x= element_blank(), #element_text(angle = 90, hjust = 1, vjust = .5),
          axis.ticks.x = element_blank(),
          legend.text=element_text(size=20),
          axis.line = element_line(colour = "black"),
          strip.text.x = element_text(size=20),
          axis.text = element_text(size=15, color = 'black'),
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size=12)) 
  
  FigA1 = FigA + geom_boxplot(data = pre_rec %>% filter(score_desc=="Precision"), mapping = aes(x = met, y = 100*score, fill= met)) + ggtitle("") +
    theme(plot.title = element_text(hjust=0.5, face = "bold", size = 12), axis.title = element_text(size = 10, face = "bold")) + ylab("Precision (%)")
  
  FigA2 = FigA + geom_boxplot(data = pre_rec %>% filter(score_desc=="Recall"), mapping = aes(x = met, y = 100*score, fill= met)) + ylab("Recall (%)") + ggtitle("") + 
    theme(plot.title = element_text(hjust=0.5, face = "bold", size = 12), axis.title = element_text(size = 10, face = "bold"))
  dat <- ggplot_build(FigA2)$data[[1]] %>% mutate(fill=factor(x = fill, levels = c("#c87e7e","#56A3E9"), labels = c("CxC","FastGxC")))
  FigA2 = FigA2 + geom_segment(data=dat, aes(x=xmin, xend=xmax,  y=middle, yend=middle, color=fill), size=4) + 
    scale_color_manual(name="",values = c("CxC"="#c87e7e","FastGxC"="#56A3E9"))
  
  FigA3 = FigA + geom_boxplot(data =   ranks %>% mutate(met = if_else(met=="FastGxE","FastGxC",met)) %>% 
                                filter(relevant_tissue_4_trait) %>% 
                                mutate(score_desc="Rank of relevant tissue", score=rank_or),
                              mapping = aes(x = met, y = score, fill= met)) + ylab("Relevant tissue rank") + ggtitle("") + 
    theme(plot.title = element_text(hjust=0.5, face = "bold", size = 12), axis.title = element_text(size = 10, face = "bold"))
  
  FigA4 = FigA + geom_boxplot(data =  cont_tables %>% 
                                mutate(met = if_else(met=="FastGxE","FastGxC",met)) %>% 
                                mutate(score_desc="Number of enriched tissues", score=FP + TP) 
                              , mapping = aes(x = met, y = score, fill= met)) + ylab("Nr enriched tissues") + ggtitle("")+
    theme(plot.title = element_text(hjust=0.5, face = "bold", size = 12), axis.title = element_text(size = 10, face = "bold"))
  
  
  fig7A=gridExtra::grid.arrange(FigA1,FigA2,FigA3,FigA4, nrow=1)
  #fig7A=gridExtra::grid.arrange(FigA1,FigA3,FigA4, nrow = 1)
  fig7A=gridExtra::grid.arrange(FigA1,FigA2,FigA3,FigA4, nrow=1)
  ggsave(plot = fig7A,
         filename = 'reviews/Figure07A_with_schiz.jpg', 
         width = 7,
         height = 4)
  
  
  
  # # manuscript claims : median precision per trait
  if(FALSE){
    pre_rec %>%  filter(score_desc=="Precision") %>% group_by(met)  %>% summarise(med=median(score))
    pre_rec %>%  filter(score_desc=="Precision") %>% group_by(tra) %>% summarise(ratio = score[1]/score[2]) %>% summarise(med=median(ratio, na.rm = T))
    
    ranks %>% group_by(met)  %>% summarise(med_rho=median(rank_or), med_pval=median(rank_p))
    ranks %>% group_by(tra,tis) %>% summarise(ratio = rank_or[2]/rank_or[1])%>% ungroup() %>% summarise(med=median(ratio, na.rm = T))
    
    cont_tables %>% 
      mutate(score=FP + TP) %>% 
      group_by(met)  %>% 
      summarise(med=median(score))
  }
  
  # Figure B: Table with most interesting gwas traits ----
  if(1){
    example_traits = c("breast carcinoma", "lung adenocarcinoma", "melanoma", "prostate carcinoma","Testicular Germ Cell Tumor", "cancer", "coronary artery disease","atrial fibrillation") #, "hypertension"
    
    table_dat = all_gwas_enrich_sig_results %>% 
      filter(tra %in% unique(relevant_tissues_annot$trait)) %>% 
      filter(tra %in% example_traits) %>%
      mutate(tra = factor(x = tra, levels = example_traits)) %>%
      filter(f_in == "each_tissue", cat=="MAPPED_TRAIT") %>% 
      dplyr::select(met, tis, tra, OR) %>% 
      filter(OR>1) %>% #### only ORs > 1 = enriched!
      group_by(tra) %>% 
      group_modify(function(tib,key){
        
        # print(key)
        # print(tib, n=1000)
        
        # annotate list of 50 tissues with either the tissue is relevant for this trait or not
        rel_tissues = relevant_tissues_annot %>% filter(trait == key$tra) %>% pull(relevant_tissues)
        rel_tissues_tib = strsplit(gsub(" ", "", rel_tissues, fixed = TRUE), ",")[[1]] %>%
          as_tibble() %>% rename(tis = value) %>% mutate(is_relevant = T) %>% 
          dplyr::select(tis, is_relevant)
        all_49_tissues.rel = all_49_tissues %>% 
          left_join(rel_tissues_tib, by=c("tis")) %>% 
          mutate(is_relevant = case_when(is.na(is_relevant) ~ F, !is.na(is_relevant) ~ is_relevant)) %>% 
          dplyr::select(tis, is_relevant)
        
        tib_annotated_all = tib %>%
          group_by(met) %>%
          group_modify(function(tib2,key2){
            
            # print("HEY!!")
            # print(key2)
            # print(tib2)
            
            # classify each tissue as enriched/not for this method
            enriched_tissues_rank = tib2 %>% arrange(desc(OR)) %>% mutate(rank = as.character(1:n()))
            
            all_49_tissues.rel.enr = all_49_tissues.rel %>% 
              left_join(enriched_tissues_rank, by=c("tis")) %>% 
              mutate(is_enriched = !is.na(OR)) %>% 
              rename(enrich_OR = OR, enrich_rank = rank) %>% 
              dplyr::select(tis, is_enriched, enrich_OR, enrich_rank, is_relevant)
            
            # convert names to abbreviations
            all_49_tissues.rel.enr = all_49_tissues.rel.enr %>% 
              inner_join(read_csv("../../manuscript/Input_Files/FigureS1_abbreviations_and_colors/GTEx_v8.AllTissues.Colors_Abbreviations.csv", col_types = cols()) %>% 
                           dplyr::select(tissue, abbreviation),
                         by = c("tis" = "tissue")) %>%
              dplyr::select(-tis) %>% 
              rename(tis = abbreviation) %>% 
              dplyr::select(tis, is_enriched, enrich_OR, enrich_rank, is_relevant)
            
            print(all_49_tissues.rel.enr)
            
            # Relevant Tissues (Rank, Odds Ratio)
            # Not relevant Tissues (Rank, Odds Ratio)
            
            rel_description = all_49_tissues.rel.enr %>% 
              filter(is_relevant, is_enriched) %>% 
              mutate(desc = paste0(tis," (",enrich_rank,")")) %>% 
              mutate(enrich_rank = as.numeric(enrich_rank)) %>% 
              arrange(enrich_rank) %>% 
              # mutate(desc = paste0(tis," (",enrich_rank,"|",round(enrich_OR,digits=1),")")) %>% 
              pull(desc) %>%
              paste0(collapse=", ")
            if(nchar(rel_description)==0) rel_description = "No enrichments"
            
            not_rel_description = all_49_tissues.rel.enr %>% 
              filter(!is_relevant, is_enriched) %>% 
              mutate(desc = paste0(tis," (",enrich_rank,")")) %>% 
              mutate(enrich_rank = as.numeric(enrich_rank)) %>% 
              arrange(enrich_rank) %>% 
              # mutate(desc = paste0(tis," (",enrich_rank,"|",round(enrich_OR,digits=1),")")) %>% 
              pull(desc) %>% 
              paste0(collapse=", ")
            if(nchar(not_rel_description)==0) not_rel_description = "No enrichments"
            
            tribble(~rel_desc, ~notrel_desc,
                    rel_description, not_rel_description)
            
          })
        
        if(length(unique(tib_annotated_all$met))==2){
          return(tib_annotated_all)
        } else if(length(unique(tib_annotated_all$met))==1){
          if(unique(tib_annotated_all$met)=="FastGxE") method_to_make_NA = "TbT"
          if(unique(tib_annotated_all$met)=="TbT") method_to_make_NA = "FastGxE"
          bind_rows(
            tib_annotated_all,
            tribble(~met, ~rel_desc, ~notrel_desc,
                    method_to_make_NA, "No enrichments", "No enrichments")
          )
        }
        
        
      }) %>%
      ungroup %>% 
      mutate(met = case_when(met == "TbT" ~ "CxC", met == "FastGxE" ~ "FastGxE", T ~ "???"))
    
    table_dat %<>% mutate(tra=ifelse(test = tra=="Testicular Germ Cell Tumor", yes = "testicular germ cell tumor", no = as.character(tra)))
    
    # rename columns
    table_dat = table_dat %>% 
      rename("GWAS Trait" = tra,
             "Method" = met,
             "Relevant Tissue(s) Enriched (Rank)" = rel_desc,
             "Other Enriched Tissue(s)" = notrel_desc)
    
    # remove duplicate GWAS trait names
    table_dat = bind_rows(lapply(unique(table_dat$`GWAS Trait`), function(g){
      ret = bind_rows(
        table_dat[table_dat$`GWAS Trait`==g & table_dat$Method == "FastGxE",],
        bind_cols(tribble(~"GWAS Trait",""),table_dat[table_dat$`GWAS Trait`==g & table_dat$Method == "CxC",-c(1)])
      )
      return(ret)
    }))
    
    table_dat = table_dat %>% 
      mutate(Method = if_else(Method=="FastGxE","FastGxC",Method)) 
    
    # Remove rank from other tissues, keep only top N_T other tissues, add nr of other tissues enriched
    nr_other_tissues=NULL
    all_other_tissues=NULL
    N_T=50
    for(i in 1:nrow(table_dat)){
      
      other_tissues=unlist(strsplit(table_dat$`Other Enriched Tissue(s)`[i], split = ', '))
      if(any(other_tissues=="No enrichments")) nr_other_tissues = c(nr_other_tissues,0) else nr_other_tissues = c(nr_other_tissues,length(other_tissues))
      if(length(other_tissues)>N_T) other_tissues=other_tissues[1:N_T] 
      if(any(other_tissues=="No enrichments")){ 
        all_other_tissues = c(all_other_tissues, other_tissues)
      } else {
        other_tissues = gsub(pattern = paste0("\\(|\\)|",paste(1:50, collapse = "|")), replacement = "", x = other_tissues)
        other_tissues = gsub(pattern = " ", replacement = "", x = other_tissues)
        all_other_tissues = c(all_other_tissues, paste(other_tissues, collapse = ", "))
      }
      
    }
    
    table_dat$`Other Enriched Tissue(s)` = all_other_tissues
    table_dat$`Nr Other Tissues Enriched` = nr_other_tissues
    
    table_dat=table_dat[,c(1:3,5,4)]
    
    table_dat$`GWAS Trait`<-
      str_wrap(table_dat$`GWAS Trait`, width = 50) 
    
    table_dat$`Relevant Tissue(s) Enriched (Rank)`<-
      str_wrap(table_dat$`Relevant Tissue(s) Enriched (Rank)`, width = 25) 
    
    table_dat$`Other Enriched Tissue(s)` <-
      str_wrap(table_dat$`Other Enriched Tissue(s)`, width = 80) 
    
    ncols=ncol(table_dat)
    nrows=nrow(table_dat) +1
    fig7B=ggtexttable(table_dat[,1:ncols], rows = NULL, theme = ttheme("blank")) %>%
      tab_add_hline(at.row = c(1:2,seq(4,nrows,2)), row.side = "top", linewidth = 2) %>%
      tab_add_hline(at.row = 1:nrows, row.side = "top", linewidth = 2, from.col = 2) %>%
      
      tab_add_hline(at.row = nrows, row.side = "bottom", linewidth = 2) %>%
      tab_add_vline(at.column = 2:ncols, column.side = "left", from.row = 1)
    # tab_add_vline(at.column = 2:3, column.side = "right", from.row = 1)
    
    fig7B = table_cell_bg(fig7B, row = seq(2,nrows,2), column = 2, fill="#56A3E9")
    fig7B = table_cell_bg(fig7B, row = seq(3,nrows,2), column = 2, fill="#c87e7e")
  }
  
  
  ###### COLOCALIZATION (FIG 7C) ###########
  
  ##### Create Colocalization tables for figure 7
  #clpp_threshold=.01
  clpp_mod_threshold = 0.5 ### like in paper
  
  coloc_GTEx=fread('../Input_Files/Figure7_GWAS/Updated_GTEx_FastGxC_cleaned.txt') %>% mutate(study="Tissues")
  coloc_PBMC=fread('../Input_Files/Figure7_GWAS/Updated_SingleCell_FastGxC_cleaned.txt') %>% mutate(study="PBMCs")
  all_coloc = rbind(coloc_GTEx,coloc_PBMC)  %>% filter(n_snps > 40)  %>% filter(clpp_mod>=clpp_mod_threshold) 
  all_coloc_no_filt = rbind(coloc_GTEx,coloc_PBMC) 
  
  
  #### remove traits 
  traits_remove = c("Type_2_Diabetes_Mahajan_2022_EUR",
                    "Triglycerides_Graham_2021_EUR",
                    "Handgrip_Strength_Tikkanen_2018",
                    "HDL_Graham_2021_EUR",
                    "Fasting_Insulin_Chen_2021_EUR",
                    "Fasting_Glucose_Chen_2021_EUR")
  all_coloc = all_coloc %>% filter(!gwas_trait %in% traits_remove)
  all_coloc_no_filt = all_coloc_no_filt%>%filter(!gwas_trait %in% traits_remove)
  
  
  # Fig 1 : FastGxC increases number of co-localization by % in tissues and % in single cell
  # Another figure would be to show % co-localization increase (y axis) as a function of the CLPP threshold (x-axis) adding CI around % increase
  
  fig7C = all_coloc %>% 
    group_by(lead_snp, gwas_trait, method,study) %>% 
    summarize(n=n()) %>% 
    ungroup() %>%
    group_by(lead_snp,gwas_trait, study) %>% 
    summarize(n_method = n_distinct(method),
              group = case_when(n_method == 2 ~ "Both",
                                n_method == 1 & method== "Regular" ~ "CxC",
                                n_method == 1 & method == "FastGxC" ~ "FastGxC",
                                TRUE ~ NA_character_)) %>%
    
    distinct() %>% 
    ungroup() %>%
    group_by(study, group) %>% 
    summarize(value = n()) %>% 
    group_by(study) %>% 
    mutate(percentage = value/sum(value)) %>% mutate(study = factor(study, levels=(c("Tissues", "PBMCs")))) %>%
    ggplot(mapping = aes(x = study, y = percentage, fill=factor(group, levels = c("FastGxC", "CxC", "Both")))) + 
    geom_bar(position="stack", stat="identity", color = "black", width = 0.7) +
    scale_fill_manual(values = c("#56A3E9", "#c87e7e", "grey")) + 
    scale_x_discrete(labels = c("GTEx" = "Tissues", 
                                "Single-Cell" = "PBMC"))+
    guides( pattern = guide_legend(),
            fill = guide_legend())+
    theme_bw() + 
    theme(legend.position="right",
          axis.title = element_text(face = "bold", size=12,color="black"),
          axis.text = element_text(size=12,color="black"),
          axis.text.x = element_text(face = "bold"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          plot.title = element_text(hjust = 0.5,size=12),
          legend.title=element_blank(),
          legend.text=element_text(size=10),
          strip.text.y = element_text(size = 0),
          strip.background = element_blank(),
          panel.spacing = unit(2, "lines")) + 
    xlab("") + ylab("Proportion of colocs")
  
  
  #### manuscript claims:
  ## we prioritized candidate causal genes for
  candidate_causal_genes = (all_coloc %>% dplyr::select(lead_snp, gwas_trait) %>% distinct() %>% dim())[1]
  total_gwas_variants = (all_coloc_no_filt %>% dplyr::select(lead_snp,gwas_trait) %>% distinct() %>% dim())[1]
  ## colocs identified in tissues
  tissues_colocs = (all_coloc %>% filter(study == "Tissues") %>% dplyr::select(lead_snp, gwas_trait) %>% distinct() %>% dim())[1]
  tissues_colocs_perc = tissues_colocs/candidate_causal_genes
  pbmcs_colocs = (all_coloc %>% filter(study == "PBMCs") %>% dplyr::select(lead_snp, gwas_trait) %>% distinct() %>% dim())[1]
  pbmcs_colocs_perc = pbmcs_colocs/candidate_causal_genes
  
  ## nr specific eQTLs that were tested for coloc
  specific_eQTLs_tested_total = (all_coloc_no_filt %>% filter(method == "FastGxC", tissue != "Shared") %>% dplyr::select(lead_snp,gwas_trait) %>% distinct() %>% dim())[1]
  specific_eQTLs_tested_tiss = (all_coloc_no_filt %>% filter(method == "FastGxC", study == "Tissues", tissue != "Shared") %>% dplyr::select(lead_snp,gwas_trait) %>% distinct() %>% dim())[1]
  specific_eQTLs_tested_celltype = (all_coloc_no_filt %>% filter(method == "FastGxC", study == "PBMCs", tissue != "Shared") %>% dplyr::select(lead_snp,gwas_trait) %>% distinct() %>% dim())[1]
  ## nr specific eQTLs that were tested for coloc with clpp > 50%
  specific_eQTLs_tested_total_filt = (all_coloc %>% filter(method == "FastGxC", tissue != "Shared") %>% dplyr::select(lead_snp,gwas_trait) %>% distinct() %>% dim())[1]
  specific_eQTLs_tested_tiss_filt = (all_coloc %>% filter(method == "FastGxC", study == "Tissues", tissue != "Shared") %>% dplyr::select(lead_snp,gwas_trait) %>% distinct() %>% dim())[1]
  specific_eQTLs_tested_celltype_filt = (all_coloc %>% filter(method == "FastGxC", study == "PBMCs", tissue != "Shared") %>% dplyr::select(lead_snp,gwas_trait) %>% distinct() %>% dim())[1]
  
  ## nr shared eQTLs that were tested for coloc
  shared_eQTLs_tested_total = (all_coloc_no_filt %>% filter(method == "FastGxC", tissue == "Shared") %>% dplyr::select(lead_snp,gwas_trait) %>% distinct() %>% dim())[1]
  shared_eQTLs_tested_tiss = (all_coloc_no_filt %>% filter(method == "FastGxC", study == "Tissues", tissue == "Shared") %>% dplyr::select(lead_snp,gwas_trait) %>% distinct() %>% dim())[1]
  shared_eQTLs_tested_celltype = (all_coloc_no_filt %>% filter(method == "FastGxC", study == "PBMCs", tissue == "Shared") %>% dplyr::select(lead_snp,gwas_trait) %>% distinct() %>% dim())[1]
  ## nr shared eQTLs that were tested for coloc with clpp > 50%
  shared_eQTLs_tested_total_filt = (all_coloc %>% filter(method == "FastGxC", tissue == "Shared") %>% dplyr::select(lead_snp,gwas_trait) %>% distinct() %>% dim())[1]
  shared_eQTLs_tested_tiss_filt = (all_coloc %>% filter(method == "FastGxC", study == "Tissues", tissue == "Shared") %>% dplyr::select(lead_snp,gwas_trait) %>% distinct() %>% dim())[1]
  shared_eQTLs_tested_celltype_filt = (all_coloc %>% filter(method == "FastGxC", study == "PBMCs", tissue == "Shared") %>% dplyr::select(lead_snp,gwas_trait) %>% distinct() %>% dim())[1]
  
  ## dataframe to show that specific eQTLs are more likely to colocalize than shared
  totals = all_coloc_no_filt %>% filter(method == "FastGxC") %>% 
    mutate(component=ifelse(tissue=="Shared",yes = "Shared",no = "Specific")) %>% dplyr::select(lead_snp, gwas_trait,component, study) %>% group_by(component,study) %>% summarize(n = n_distinct(lead_snp, gwas_trait))
  specific_shared_coloc = all_coloc %>% filter(method == "FastGxC") %>% 
    mutate(component=ifelse(tissue=="Shared",yes = "Shared",no = "Specific")) %>% dplyr::select(lead_snp, gwas_trait,component, study) %>% group_by(component,study) %>% summarize(n = n_distinct(lead_snp, gwas_trait))
  specific_shared_coloc$totals = totals$n
  specific_shared_coloc$percentage = specific_shared_coloc$n/specific_shared_coloc$totals
  
  immune_traits <- c(
    "Systemic_Lupus_Erythematosus_Bentham_2015",
    "Rheumatoid_Arthritis_Ishigaki_2022",
    "Ulcerative_Colitis_deLange_2017",
    "Crohns_Disease_deLange_2017",
    "Psoriasis_Tsoi_2012",
    "Primary_Biliary_Cirrhosis_Cordell_2021",
    "Asthma_Tsuo_2022",
    "Multiple_Sclerosis_Andlauer_2016",
    "Hypothyroidism_Mathieu_2022",
    "Type_1_Diabetes_Chiou_2021"
  )
  totals_immune = all_coloc_no_filt %>% filter(method == "FastGxC" & gwas_trait %in% immune_traits) %>% 
    mutate(component=ifelse(tissue=="Shared",yes = "Shared",no = "Specific")) %>% dplyr::select(lead_snp, gwas_trait,component, study) %>% group_by(component,study) %>% summarize(n = n_distinct(lead_snp, gwas_trait))
  specific_shared_coloc_immune = all_coloc %>% filter(method == "FastGxC" & gwas_trait %in% immune_traits) %>% 
    mutate(component=ifelse(tissue=="Shared",yes = "Shared",no = "Specific")) %>% dplyr::select(lead_snp, gwas_trait,component, study) %>% group_by(component,study) %>% summarize(n = n_distinct(lead_snp, gwas_trait))
  specific_shared_coloc_immune$totals = totals_immune$n
  specific_shared_coloc_immune$percentage = specific_shared_coloc_immune$n/specific_shared_coloc_immune$totals
  
  final_specific_shared_coloc = bind_rows(specific_shared_coloc %>% filter(study == "Tissues"), specific_shared_coloc_immune %>% filter(study %in% "PBMCs"))
  factor(final_specific_shared_coloc$study, levels = c("Tissues", "PBMCs"))
  
  binom.test(final_specific_shared_coloc %>% filter(component == "Specific" & study == "Tissues") %>% 
               ungroup() %>% dplyr::select(n, totals) %>% unlist() %>% unname(), p = final_specific_shared_coloc %>% filter(component == "Shared" & study == "Tissues") %>% 
               ungroup() %>% dplyr::select(percentage) %>% unlist() %>% unname(), alternative = c("greater"))
  
  
  ## temp Fig 7c_b
  fig7C_b_temp = final_specific_shared_coloc  %>%
    ggplot(mapping = aes(x = factor(study, levels = c("Tissues", "PBMCs")), y = 100*percentage, alpha = factor(component, levels = c("Specific", "Shared")), fill=factor(component, levels = c("Specific", "Shared")))) + 
    geom_bar(position="dodge", stat="identity", color = "black", width = 0.7) +
    scale_fill_manual(values = c("#0f4271", "#56A3E9")) + 
    scale_alpha_manual(values = c(1,0.3, 1))+
    #facet_wrap(~factor(study, levels = c("Tissues", "PBMCs")), scales = "free")+
    #scale_y_continuous(limits = c(0,1))+
    theme_bw() + 
    theme(legend.position="top",
          axis.title = element_text(face = "bold", size=12,color="black"),
          axis.text = element_text(size=12,color="black"),
          axis.text.x = element_text(face = "bold"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          plot.title = element_text(hjust = 0.5,size=12),
          legend.title=element_blank(),
          legend.text=element_text(size=10),
          strip.text.y = element_text(size = 0),
          strip.text.x = element_blank(),
          strip.background = element_blank(),
          panel.spacing = unit(2, "lines")) + 
    xlab("") + ylab("eQTLs Colocalizing (%)")
  
  
  ## colocs identified in tissues
  coloc_table = all_coloc %>% 
    group_by(lead_snp, gwas_trait, method,study) %>% 
    summarize(n=n()) %>% 
    ungroup() %>%
    group_by(lead_snp,gwas_trait, study) %>% 
    summarize(n_method = n_distinct(method),
              group = case_when(n_method == 2 ~ "Both",
                                n_method == 1 & method== "Regular" ~ "CxC",
                                n_method == 1 & method == "FastGxC" ~ "FastGxC",
                                TRUE ~ NA_character_)) %>%
    
    distinct() %>% 
    ungroup() %>%
    group_by(study, group) %>% 
    summarize(value = n()) %>% 
    group_by(study) %>% 
    mutate(percentage = value/sum(value)) %>% mutate(study = factor(study, levels=(c("Tissues", "PBMCs"))))
  
  ## colocs identified by both methods: 
  ## PBMCs - 
  tissues_both = 
    unname(unlist(coloc_table[coloc_table$group == "Both" & coloc_table$study == "Tissues", "value"]))/sum(coloc_table[coloc_table$study == "Tissues", "value"])
  pbmcs_both = 
    unname(unlist(coloc_table[coloc_table$group == "Both" & coloc_table$study == "PBMCs", "value"]))/sum(coloc_table[coloc_table$study == "PBMCs", "value"])
  
  tissues_fastgxc = 
    unname(unlist(coloc_table[coloc_table$group == "FastGxC" & coloc_table$study == "Tissues", "value"]))/sum(coloc_table[coloc_table$study == "Tissues", "value"])
  pbmcs_fastgxc = 
    unname(unlist(coloc_table[coloc_table$group == "FastGxC" & coloc_table$study == "PBMCs", "value"]))/sum(coloc_table[coloc_table$study == "PBMCs", "value"])
  
  ## increase:
  tissues_increase = 
    unname(unlist(coloc_table[coloc_table$group == "FastGxC" & coloc_table$study == "Tissues", "value"]))/(566+4711)
  pbmcs_increase = 
    unname(unlist(coloc_table[coloc_table$group == "FastGxC" & coloc_table$study == "PBMCs", "value"]))/(817+239)
  
  ## contribution of shared and specific eQTLs to colocs identified by FastGxC
  ## from Figure 7C_b
  sh_sp_coloc_table = all_coloc %>% filter(method=="FastGxC") %>% 
    mutate(component=ifelse(tissue=="Shared",yes = "Shared",no = "Specific")) %>% 
    
    
    group_by(lead_snp, component,study) %>% 
    summarize(n=n()) %>% 
    
    group_by(lead_snp,study) %>% 
    summarize(n_method = n_distinct(component),
              group = case_when(n_method == 2 ~ "Both",
                                n_method == 1 & component== "Specific" ~ "Specific",
                                n_method == 1 & component == "Shared" ~ "Shared",
                                TRUE ~ NA_character_)) %>%
    distinct() %>% 
    ungroup() %>%
    group_by(study, group) %>% 
    summarize(value = n()) %>% 
    group_by(study) %>% 
    mutate(percentage = value/sum(value)) %>% mutate(study = factor(study, levels=(c("Tissues", "PBMCs"))))
  
  tissues_sh_sp = 2461/(2461+542+1533)
  pbmcs_sh = 404/(322+405+260)
  pbmcs_sh_sp = 322/(322+405+260)
  
  #### Figure S18: 
  if(1){
    # Fig 2 : Context-specific eQTLs increase co-localization over shared by % in tissues and % in cell types
    # Another figure would be to show % co-localization increase (y axis) as a function of the CLPP threshold (x-axis)
    fig7C_b = all_coloc %>% filter(method=="FastGxC") %>% 
      mutate(component=ifelse(tissue=="Shared",yes = "Shared",no = "Specific")) %>% 
      
      
      group_by(lead_snp, component,study) %>% 
      summarize(n=n()) %>% 
      
      group_by(lead_snp,study) %>% 
      summarize(n_method = n_distinct(component),
                group = case_when(n_method == 2 ~ "Both",
                                  n_method == 1 & component== "Specific" ~ "Specific",
                                  n_method == 1 & component == "Shared" ~ "Shared",
                                  TRUE ~ NA_character_)) %>%
      distinct() %>% 
      ungroup() %>%
      group_by(study, group) %>% 
      summarize(value = n()) %>% 
      group_by(study) %>% 
      mutate(percentage = value/sum(value)) %>% mutate(study = factor(study, levels=(c("Tissues", "PBMCs")))) %>%
      ggplot(mapping = aes(x = study, y = percentage, alpha = factor(group, levels = c("Specific", "Shared", "Both")), fill=factor(group, levels = c("Specific", "Shared", "Both")))) + 
      geom_bar(position="stack", stat="identity", color = "black", width = 0.7) +
      scale_fill_manual(values = c("#56A3E9", "#56A3E9", "white")) + 
      scale_alpha_manual(values = c(1,0.3, 1))+
      scale_x_discrete(labels = c("GTEx" = "Tissues", 
                                  "Single-Cell" = "PBMC"))+
      guides( pattern = guide_legend(),
              fill = guide_legend())+
      theme_bw() + 
      theme(legend.position="right",
            axis.title = element_text(face = "bold", size=12,color="black"),
            axis.text = element_text(size=12,color="black"),
            axis.text.x = element_text(face = "bold"),
            axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            plot.title = element_text(hjust = 0.5,size=12),
            legend.title=element_blank(),
            legend.text=element_text(size=10),
            strip.text.y = element_text(size = 0),
            strip.background = element_blank(),
            panel.spacing = unit(2, "lines")) + 
      xlab("") + ylab("Proportion of FastGxC Colocs")
    
    ggsave(plot = fig7C_b,
           filename = 'FigureS18_FastGxC_coloc.jpg', 
           width = 6,
           height = 7)
  }
  
  ### of the FastGxC unique colocs, what is the sharing and specificity?
  fig7C_c = all_coloc %>% 
    group_by(lead_snp, gwas_trait, method,study) %>% 
    mutate(n=n()) %>% 
    ungroup() %>%
    group_by(lead_snp,gwas_trait, study) %>% 
    mutate(n_method = n_distinct(method),
           group = case_when(n_method == 2 ~ "Both",
                             n_method == 1 & method== "Regular" ~ "CxC",
                             n_method == 1 & method == "FastGxC" ~ "FastGxC",
                             TRUE ~ NA_character_)) %>% filter(group == "FastGxC") %>% 
    dplyr::select(lead_snp, gwas_trait, study, group, tissue) %>% 
    mutate(component=ifelse(tissue=="Shared",yes = "Shared",no = "Specific")) %>% 
    
    group_by(lead_snp,gwas_trait,study) %>% 
    summarize(n_method = n_distinct(component),
              group = case_when(n_method == 2 ~ "Both",
                                n_method == 1 & component== "Specific" ~ "Specific",
                                n_method == 1 & component == "Shared" ~ "Shared",
                                TRUE ~ NA_character_)) %>%
    distinct() %>% 
    ungroup() %>%
    group_by(study, group) %>% 
    summarize(value = n()) %>% 
    group_by(study) %>% 
    mutate(percentage = value/sum(value)) %>% mutate(study = factor(study, levels=(c("Tissues", "PBMCs")))) %>%
    ggplot(mapping = aes(x = study, y = percentage, alpha = factor(group, levels = c("Specific", "Shared", "Both")), fill=factor(group, levels = c("Specific", "Shared", "Both")))) + 
    geom_bar(position="stack", stat="identity", color = "black", width = 0.7) +
    scale_fill_manual(values = c("#56A3E9", "#56A3E9", "white")) + 
    scale_alpha_manual(values = c(1,0.3, 1))+
    guides( pattern = guide_legend(),
            fill = guide_legend())+
    theme_bw() + 
    theme(legend.position="right",
          axis.title = element_text(face = "bold", size=12,color="black"),
          axis.text = element_text(size=12,color="black"),
          axis.text.x = element_text(face = "bold"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          plot.title = element_text(hjust = 0.5,size=12),
          legend.title=element_blank(),
          legend.text=element_text(size=10),
          strip.text.y = element_text(size = 0),
          strip.background = element_blank(),
          panel.spacing = unit(2, "lines")) + 
    xlab("") + ylab("Proportion of unique FastGxC colocs")
  
  
  
  legend = get_legend(fig7C + theme(legend.position = "top", legend.text = element_text(size = 20)))
  legend_sh_sp = get_legend(fig7C_b_temp + theme(legend.position = "top", legend.text = element_text(size = 20)))
  
  combineLegend = ggdraw() +
    draw_plot(legend, x = 0.25, y = 0.3, width = 0.3, height = 0.15) +
    draw_plot(legend_sh_sp, x = 0.50, y = 0.3, width = 0.3, height = 0.15)
  
  fig7 = plot_grid(
    plot_grid(
      combineLegend 
    ),
    ggdraw() + 
      draw_plot(plot = fig7A, x = 0.01, y = 0.77, height = 0.2, width = 0.59) +
      draw_plot(plot = fig7B, x = 0.12, y = 0, width = 0.8) +
      draw_plot(plot = fig7C + theme(legend.position = "none"), x = 0.62, y = 0.75, width = 0.15, height = 0.2)+
      draw_plot(plot = fig7C_b_temp + theme(legend.position = "none"), x = 0.77, y = 0.75, width = 0.22, height = 0.2)+
      draw_plot_label(label = c('A','B', 'C'),  x = c(0,0, 0.6), y = c(0.98,0.76,0.98), size = 30),
    
    ncol = 1,
    rel_heights = c(0.3, 12)
  )
  
  #fig7 = plot_grid(fig7, barplot_legend, ncol = 1, rel_heights = c(1, .05))
  ggsave(filename = 'reviews/Fig07_final_with_schiz.pdf', plot = fig7, width = 15, height = 15, limitsize = F)
  #ggsave(filename = 'Fig07_GWAS.pdf', plot = fig7, width = 15, height = 15, limitsize = F)
  
}  

#### Figure S17: FastGxC increase across CLPP thresholds
if(1){
  clpp_mod_threshold = 0.5 ### like in paper
  
  coloc_GTEx=fread('../Input_Files/Figure7_GWAS/Updated_GTEx_FastGxC_cleaned.txt') %>% mutate(study="Tissues")
  coloc_PBMC=fread('../Input_Files/Figure7_GWAS/Updated_SingleCell_FastGxC_cleaned.txt') %>% mutate(study="PBMCs")
  all_coloc = rbind(coloc_GTEx,coloc_PBMC)  %>% filter(n_snps > 40)  %>% filter(clpp_mod>=clpp_mod_threshold) 
  all_coloc_no_filt = rbind(coloc_GTEx,coloc_PBMC) %>% filter(n_snps > 40)
  
  
  get_bin = function(df, threshold, z = zval){
    cur_df = df %>% filter(clpp_mod >= threshold) %>% group_by(study) %>% mutate(n = n_distinct(lead_snp, gwas_trait)) %>%
      ungroup() %>% group_by(lead_snp, gwas_trait, study) %>%
      mutate(n_method = n_distinct(method),
             group = case_when(n_method == 2 ~ "Both",
                               n_method == 1 & method== "Regular" ~ "CxC",
                               n_method == 1 & method == "FastGxC" ~ "FastGxC",
                               TRUE ~ NA_character_)) %>% dplyr::select(c("lead_snp", "gwas_trait", "group", "study", "n")) %>%
      distinct() %>% ungroup() %>%
      group_by(study, group) %>%
      mutate(value = n()) %>% dplyr::select(c("study", "group", "value", "n")) %>% distinct() %>%
      group_by(study) %>% mutate(increase = value/(sum(value) - value)) %>%
      filter(group == "FastGxC") %>% mutate(clpp = threshold) %>% ungroup() %>% group_by(study) %>%
      mutate(upper_ci = (increase+(z*sqrt((1-increase)/n))), lower_ci = (increase-(z*sqrt((1-increase)/n))))
    return(cur_df)
  }
  
  clpp_thresholds = seq(0.1,0.95, by = 0.05)
  final_df = data.frame()
  for(cur_clpp in clpp_thresholds){
    cur_df = get_bin(all_coloc_no_filt, cur_clpp)
    final_df = rbind(final_df, cur_df)
  }
  
  plot = ggplot(final_df, aes(x = clpp, y = increase, color = study, alpha = study)) + geom_point() +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), size = .5) +
    theme_bw() +
    scale_color_manual(values = c("purple", "purple"))+
    scale_alpha_manual(values = c(0.2, 1))+
    theme(legend.position = "top",
          axis.title = element_text(face = "bold", size=25,color="black"),
          axis.text = element_text(size=21,color="black"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          plot.title = element_text(hjust = 0.5,size=15),
          legend.title=element_blank(),
          legend.text=element_text(size=20),
          plot.margin = unit(c(1,1,1,1), "cm"))+
    #strip.background = element_blank())+
    xlab("CLPP_mod threshold") + ylab("Proportion Increase in Colocalization")
  
  ggsave(filename = 'FigureS17_clpp_mod_thresh.pdf', plot = plot, width = 10, height = 8)
}

####################################################################
########################## Figure S13 ###############################
### Correlation of expression PCs with covariates in GTEx. 
####################################################################
if(1){
  
  # read PCs on each expression component, i.e. full, tissue-shared, and tissue-specific
  PCs=data.frame(read.table(file = "../Input_Files/SuppFigure13_PCA/PCs.v8.EUR.normalized_expression_full.txt", header = T, sep = '\t'), row.names = 1)
  
  bPCs=data.frame(read.table(file = "../Input_Files/SuppFigure13_PCA/PCs.v8.EUR.normalized_expression_homogeneous.txt", header = T, sep = '\t'), row.names = 1)
  bPCs=bPCs[matrix(data = unlist(strsplit(rownames(PCs), split = " - ")), ncol = 2, byrow = T)[,1],]
  rownames(bPCs)=rownames(PCs)
  
  wPCs=data.frame(read.table(file ="../Input_Files/SuppFigure13_PCA/PCs.v8.EUR.normalized_expression_heterogeneous.txt", header = T, sep = '\t'), row.names = 1)
  
  # read covariates 
  covs=read.table(file = "../Input_Files/SuppFigure13_PCA/GTEx_v8_all_covariates.txt", header = T, sep = '\t')
  cov_attr=read.table(file = "../Input_Files/SuppFigure13_PCA/GTEx_v8_all_covariates_description.txt", header = T, sep = '\t', stringsAsFactors = F)
  old_names=c(colnames(covs)[colnames(covs) %in% row.names(cov_attr)],"pcr","COHORT","SEX","AGE","ETHNCTY","DTHHRDY","HGHT","WGHT")
  new_names=c(cov_attr[colnames(covs)[colnames(covs) %in% row.names(cov_attr)],1],"PCR", "Cohort", "Sex", "Age", "Ethnicity", "Death cause","Height", "Weight")
  names(new_names)=old_names
  
  C=read.table(file = "../Input_Files/SuppFigure13_PCA/corr_PCs_covs.txt", header = T, sep = "\t")[-55,-55]
  C=C[rownames(C)!="SMGEBTCH",colnames(C)!="SMGEBTCH"]
  rownames(C)[rownames(C) %in% old_names] = new_names[rownames(C)[rownames(C) %in% old_names] ]
  
  bC=read.table(file = "../Input_Files/SuppFigure13_PCA/corr_bPCs_covs.txt", header = T, sep = "\t")[-55,-55]
  bC=bC[rownames(bC)!="SMGEBTCH",colnames(bC)!="SMGEBTCH"]
  rownames(bC)[rownames(bC) %in% old_names] = new_names[rownames(bC)[rownames(bC) %in% old_names] ]
  
  wC=read.table(file = "../Input_Files/SuppFigure13_PCA/corr_wPCs_covs.txt", header = T, sep = "\t")[-55,-55]
  wC=wC[rownames(wC)!="SMGEBTCH",colnames(wC)!="SMGEBTCH"]
  rownames(wC)[rownames(wC) %in% old_names] = new_names[rownames(wC)[rownames(wC) %in% old_names] ]
  
  vars=c("gPC1", "gPC2", "gPC3", "Cohort", "Sex",  "Age", "Ethnicity", "Height", "Weight", "BMI", "Death cause",  "Ischemic time", "RIN Number",  "Intragenic Rate", "Total Mapping Rate", "Exonic Mapping Rate",  "rRNA Rate") #, "gPC4", "gPC5", "PCR", "Uberon ID", "Tissue Type", "Tissue"
  
  n_PCs_to_show=10
  
  # old version
  if(0){
    # p1=ggcorrplot(corr = C[vars,colnames(PCs)[1:n_PCs_to_show]], show.legend = F) + 
    #   theme(plot.title = element_text(hjust = .5, size = 15),
    #         axis.text.x = element_text(angle=90, vjust = .5, hjust = 1), 
    #         axis.text = element_text(color = "black")) + 
    #   ggtitle("Tissue Expression") 
    # p2=ggcorrplot(corr = bC[vars,colnames(bPCs)[1:n_PCs_to_show]], show.legend = F)+ 
    #   theme(plot.title = element_text(hjust = .5, size = 15), 
    #         axis.text = element_text(color = "black"),
    #         axis.text.x = element_text(angle=90, vjust = .5, hjust = 1)) + 
    #   ggtitle("Tissue-Shared Expression")
    # 
    # p3=ggcorrplot(corr = wC[vars,colnames(wPCs)[1:n_PCs_to_show]], show.legend = T) + 
    #   theme(plot.title = element_text(hjust = .5, size = 15), 
    #         axis.text = element_text(color = "black"),
    #         axis.text.x = element_text(angle=90, vjust = .5, hjust = 1)) + 
    #   ggtitle("Tissue-Specific Expression") 
    
    # fig3A= ggdraw() +
    #   draw_plot(plot = p1, x =.0,  y = 0, width = .32, height = 1) +
    #   draw_plot(plot = p2, x =.32,  y = 0, width = .32, height = 1) +
    #   draw_plot(plot = p3, x =.64,  y = 0,  width = .37, height = 1)
  }
  
  breaksList = seq(0, 1, by = 0.00001)
  
  p1 = pheatmap(t(C[vars,colnames(PCs)[n_PCs_to_show:1]]),
                cluster_rows = F,
                cluster_cols = F,
                color = colorRampPalette(c("white",brewer.pal(n = 6, name = "Reds")))(length(breaksList)),
                breaks = breaksList,
                legend = F, 
                silent = T,
                show_colnames = F,
                fontsize = 20, fontsize_row = 10, fontsize_col = 10, angle_col = 45, 
                main = "Tissue Expression")
  
  p2 = pheatmap(t(bC[vars,colnames(bPCs)[n_PCs_to_show:1]]),
                cluster_rows = F,
                cluster_cols = F,
                color = colorRampPalette(c("white",brewer.pal(n = 6, name = "Reds")))(length(breaksList)),
                breaks = breaksList,
                legend = F,
                silent = T,
                show_colnames = F,
                fontsize = 20, fontsize_row = 10, fontsize_col = 10, angle_col = 45, 
                main = "Tissue-Shared Expression")
  
  p3 = pheatmap(t(wC[vars,colnames(wPCs)[n_PCs_to_show:1]]),
                cluster_rows = F,
                cluster_cols = F,
                color = colorRampPalette(c("white",brewer.pal(n = 6, name = "Reds")))(length(breaksList)),
                breaks = breaksList,
                legend = F,
                silent = T,
                fontsize = 20, fontsize_row = 10, fontsize_col = 10, angle_col = 45, 
                main = "Tissue-Specific Expression")
  
  # fig3A = plot_grid(p1$gtable, p2$gtable, p3$gtable, rel_widths = c(1,1,1.1),nrow = 1)
  figS1 = plot_grid(p1$gtable, NULL, p2$gtable, NULL, p3$gtable, 
                    ncol = 1,
                    rel_heights = c(1, 0.12, 1, 0.12, 1.32))
  
  # This code can be used to grab the legend for p3 and make it's own plot if we want to put it somewhere else
  # library(ggpubr)
  # legend=get_legend(p3)
  # as_ggplot(legend)
  
  
  ggsave(plot = figS1,
         filename = 'FigureS13_heatmaps.jpg', 
         width = 6,
         height = 12)
  
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Table S3: gtex color + abbreviation
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(1){
  p = inner_join(gtex_colors.tib, gtex_abbrev.tib, by=c("tissue")) %>% 
    mutate(tissue = if_else(tissue=="AverageTissue","Shared",tissue)) %>% 
    formattable(list("color_hex" = custom_tis_color()))
  sc = inner_join(sc_colors.tib, sc_abbrev.tib, by=c("tissue")) %>% rename(cell_type = "tissue") %>%
    mutate(cell_type = if_else(cell_type=="AverageContext","Shared",cell_type)) %>% 
    formattable(list("color_hex" = custom_tis_color()))
  
  export_formattable(f = p,file = "../SuppTables/TableS3_GTEx_colors.png")
  export_formattable(f = sc, width = "50%", file = "../SuppTables/TableS3_sc_colors.png")
}

### Supp Table S4
if(1){
  outdir = "../SuppTables/"
  egenes_gtex = read_csv("../Input_Files/Figure3_Performance/eGenes.v8.EUR.all_tissues.residualized_exp_types.txt") %>%
    mutate(exp_type = case_when(exp_type == "normalized_and_residualized_expression" ~ "CxC eGene",
                                exp_type == "normalized_and_residualized_expression_heterogeneous" ~ "FastGxC sp-eGene",
                                exp_type == "normalized_and_residualized_expression_homogeneous" ~ "FastGxC sh-eGene",
                                TRUE ~ exp_type)) %>%
    rename(eGene_type = exp_type, eGene = gene) 
  
  egenes_sc = read_csv("../Input_Files/Figure3_Performance/eGenes.scMeta.all_contexts.residualized_exp_types.txt") %>%
    rename(eGene_type = exp_type, celltype = tissue, eGene = gene) %>% 
    mutate(eGene_type = case_when(eGene_type == "mean_norm_res_exp.specific" ~ "FastGxC sp-eGene",
                                  eGene_type == "mean_norm_res_exp.shared" ~ "FastGxC sh-eGene",
                                  eGene_type == "mean_norm_res_exp" ~ "CxC eGene"))
  x = c("egenes_gtex", "egenes_sc")
  WriteXLS(x, ExcelFileName = paste0(outdir, "TableS4_FastGxC_CxC_eGenes_GTEx_and_singlecell.xlsx"), SheetNames = c("GTEx", "Single-Cell"), row.names = FALSE, col.names = TRUE)
}

### Supp Table S6
if(1){
  outdir = "../SuppTables/"
  coloc_GTEx=fread('../Input_Files/Figure7_GWAS/Updated_GTEx_FastGxC_cleaned.txt') %>% mutate(study="Tissues")
  coloc_PBMC=fread('../Input_Files/Figure7_GWAS/Updated_SingleCell_FastGxC_cleaned.txt') %>% mutate(study="PBMCs")
  
  traits_remove = c("Type_2_Diabetes_Mahajan_2022_EUR",
                    "Triglycerides_Graham_2021_EUR",
                    "Handgrip_Strength_Tikkanen_2018",
                    "HDL_Graham_2021_EUR",
                    "Fasting_Insulin_Chen_2021_EUR",
                    "Fasting_Glucose_Chen_2021_EUR")
  coloc_GTEx = coloc_GTEx %>% filter(!gwas_trait %in% traits_remove)
  coloc_PBMC = coloc_PBMC %>%filter(!gwas_trait %in% traits_remove)
  
  coloc_GTEx = coloc_GTEx %>% mutate(method = case_when(method == "Regular" ~ "CxC",
                                                        TRUE ~ method)) %>%
    rename(gene = feature) %>% select(!c("neg_log_gwas_pval", "neg_log_qtl_pval"))
  
  coloc_PBMC = coloc_PBMC %>% mutate(method = case_when(method == "Regular" ~ "CxC",
                                                        TRUE ~ method)) %>%
    rename(gene = feature, celltype = tissue) %>% dplyr::select(!c("neg_log_gwas_pval", "neg_log_qtl_pval"))
  datasets = list("GTEx" = coloc_GTEx, "Single-Cell" = coloc_PBMC)
  write.xlsx(datasets, file = paste0(outdir, "TableS6_FastGxC_CxC_coloc_GTEx_and_singlecell.xlsx"))
}
