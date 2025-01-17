setwd("/Users/ejagoda/Documents/HGRM/Crispri/dc_tap_paper_uploads/")
library(stringr)
library(ggplot2)
library (reshape2)
library(gridExtra)
library(grid)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(cowplot)
library(RColorBrewer)

mytheme <- theme_classic() + theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20))
#Fig 1c left crop vs sg-opti-cs
#realized we can just put these directly together
#okay this is gonna be umis per guide per 1000 reads

p10_guide_tab = read.table("figures/inputs_for_fig1/P10_tap_all_guide_info.txt",header=T,sep = '\t')
p10_reads = 1677.587 #reads per cell

#33606943, total reads per cell  #1677.587, 

lsg_reads = 1091.170592310356772 #reads per cell
lsg_guides_per_cell = read.table("figures/inputs_for_fig1/protospacer_calls_per_cell_052322_analysis_L-SG-TAP-Whole.csv",header=T,sep = ',')
total_lsg_cells = 17322 
lsg_n_cells_no_guide = total_lsg_cells-nrow(lsg_guides_per_cell)

#for p10 we don't have every guide to so to ensure that it's reads per guide doing max_umis x fraction_top_guide
p10_guide_tab$umis_per_top_guide = p10_guide_tab$max_umis * p10_guide_tab$fraction_top_guide
p10_guide_tab$umis_per_top_guide[p10_guide_tab$umis_per_top_guide == "NaN"] <- 0

combined_tab = data.frame(p10_guide_tab[,"umis_per_top_guide"])
colnames(combined_tab) = "total_guide_umis"
combined_tab$experiment = "Crop"
combined_tab$num_features = 1

lsg_guides_per_cell$experiment = "hU6"
colnames(lsg_guides_per_cell)[4] = "total_guide_umis"
lsg_guides_per_cell_small = lsg_guides_per_cell[,c("total_guide_umis","experiment","num_features")]
lsg_nulls = data.frame(cbind(rep(0,lsg_n_cells_no_guide),rep("hU6",lsg_n_cells_no_guide),rep(0,lsg_n_cells_no_guide)))
colnames(lsg_nulls)= c("total_guide_umis","experiment","num_features")
lsg_guides_per_cell_small = rbind(lsg_nulls,lsg_guides_per_cell_small)

dc_file <- "figures/inputs_for_fig1/protospacer_calls_per_cell_052322_analysis_H-DC-TAP-Whole.csv"
dc_data <- read.csv(dc_file)
dc_mean_reads_per_cell <- 1363
dc_total_cells = 25669
dc_n_cells_no_guide = dc_total_cells-nrow(dc_data)

dc_data$experiment = "mU6"
colnames(dc_data)[4] = "total_guide_umis"
dc_data_small = dc_data[,c("total_guide_umis","experiment","num_features")]
dc_data_nulls = data.frame(cbind(rep(0,dc_n_cells_no_guide),rep("mU6",dc_n_cells_no_guide),rep(0,dc_n_cells_no_guide)))
colnames(dc_data_nulls)= c("total_guide_umis","experiment","num_features")
dc_guides_per_cell_small = rbind(dc_data_nulls,dc_data_small)


## this tab has total per guide cell guide umis 
combined_tab = data.frame(rbind(combined_tab,lsg_guides_per_cell_small,dc_guides_per_cell_small))
#nom

#weight reads to per 1000 reads per cell to make
#need to sum umis across guides for cells with multiplue guides #divide by n guides?
combined_tab$total_guide_umis_real = "x"
combined_tab$total_guide_umis_real_weighted = "x"

for (i in 1:nrow(combined_tab)){
  if (grepl(combined_tab$total_guide_umis[i],pattern = "\\|")){
    total_umis = sum(as.numeric(paste0(str_split(combined_tab$total_guide_umis[i],"\\|")[[1]])))
    umis = total_umis/as.numeric(paste0(combined_tab$num_features[i]))
    #umis = total_umis/co
    combined_tab$total_guide_umis_real[i] = umis
  }else{
    combined_tab$total_guide_umis_real[i] = combined_tab$total_guide_umis[i]
  }
  if (combined_tab$experiment[i] == "Crop"){
    combined_tab$total_guide_umis_real_weighted[i] = (as.numeric(paste0(combined_tab$total_guide_umis_real[i])) / p10_reads) * 1000 #convert to per 1000 reads, by dividing by reads per cell * 1000 reads
  } else if (combined_tab$experiment[i] == "hU6"){
    combined_tab$total_guide_umis_real_weighted[i] = (as.numeric(paste0(combined_tab$total_guide_umis_real[i])) / lsg_reads) * 1000
  }else if (combined_tab$experiment[i] == "mU6"){
    combined_tab$total_guide_umis_real_weighted[i] = (as.numeric(paste0(combined_tab$total_guide_umis_real[i])) / dc_mean_reads_per_cell) * 1000
  }
}

combined_tab$experiment = factor(combined_tab$experiment, levels = c("Crop","mU6","hU6"))

p_guides = ggplot(combined_tab, aes(y = log10(as.numeric(paste0(total_guide_umis_real_weighted))), x = experiment,fill = experiment)) +
  geom_boxplot()+
  labs(y = "log10(Guide UMIs Guide Per Cell per 1k reads)")+
  mytheme

means = by(as.numeric(paste0(na.omit(combined_tab$total_guide_umis_real_weighted))),combined_tab$experiment,mean)


#combined_tab$experiment: Crop
#[1] 12.90019
#------------------------------------------------------------------------------------------ 
#  combined_tab$experiment: mU6
#[1] 55.2867
#------------------------------------------------------------------------------------------ 
#  combined_tab$experiment: hU6
#[1] 88.04974


png("figures/Fig1c.png",height = 400, width = 400)
print(p_guides)
dev.off()

pdf("figures/Fig1c.pdf",height = 4, width = 4)
print(p_guides)
dev.off()

#multi moi figure 1d
#k562
tab = read.table("figures/inputs_for_fig1/k562_moi_merged_all_experiments_GATA1_data.txt",header = T,sep = '\t')
tab = tab[tab$experiment != "moi5_multi_moi" & tab$experiment != "moi3_multi_moi",]
exp_order = c("moi1_lsg","moi1_multi_moi","moi3_multi_moi","moi5_multi_moi", "moi6_multi_moi", "moi9_multi_moi", "moi14_multi_moi","moi20_hsg")


merge = tab
merge$scaled_to_per_exp_moi1 = "x"
merge$scaled_to_mean_exp_moi1 = "x"
merge$ci_high_scaled_to_per_exp_moi1 = "x"
merge$ci_low_scaled_to_per_exp_moi1 = "x"
merge$scaled_to_mean_exp_moi1 = "x"
merge$ci_high_to_mean_exp_moi1 = "x"
merge$ci_low_to_mean_exp_moi1 = "x"
merge = unique(merge)
for (i in 1:nrow(merge)){
  if (merge[i,"experiment"] == "moi1_lsg" | merge[i,"experiment"] == "moi1_multi_moi"){
    merge[i,"scaled_to_per_exp_moi1"] = 1
    scale_factor_per_mean_exp_moi1 = mean(as.numeric(paste0(merge[(merge$experiment == "moi1_lsg" | merge$experiment == "moi1_multi_moi")  & merge$gene_pert == merge$gene_pert[i],"output_lfcs"])))
    merge[i,"scaled_to_mean_exp_moi1"] = as.numeric(paste0(merge[i,"output_lfcs"]))/scale_factor_per_mean_exp_moi1 
    merge[i,"ci_high_to_mean_exp_moi1"] =  as.numeric(paste0(merge[i,"ci_highs"]))/scale_factor_per_mean_exp_moi1 
    merge[i,"ci_low_to_mean_exp_moi1"] =  as.numeric(paste0(merge[i,"ci_lows"]))/scale_factor_per_mean_exp_moi1 
  }
  else if (merge[i,"experiment"] == "moi20_hsg") {
    merge[i,"scaled_to_per_exp_moi1"] = as.numeric(paste0(merge[i,"output_lfcs"]))/as.numeric(paste0(merge[merge$experiment == "moi1_lsg" & merge$gene_pert == merge$gene_pert[i],"output_lfcs"]))
    scale_factor_per_mean_exp_moi1 = mean(as.numeric(paste0(merge[(merge$experiment == "moi1_lsg" | merge$experiment == "moi1_multi_moi")  & merge$gene_pert == merge$gene_pert[i],"output_lfcs"])))
    merge[i,"scaled_to_mean_exp_moi1"] = as.numeric(paste0(merge[i,"output_lfcs"]))/scale_factor_per_mean_exp_moi1
    merge[i,"ci_high_to_mean_exp_moi1"] =  as.numeric(paste0(merge[i,"ci_highs"]))/scale_factor_per_mean_exp_moi1 
    merge[i,"ci_low_to_mean_exp_moi1"] =  as.numeric(paste0(merge[i,"ci_lows"]))/scale_factor_per_mean_exp_moi1 
    
  }
  else if (grepl(merge[i,"experiment"],pattern =  "multi_moi")) {
    merge[i,"scaled_to_per_exp_moi1"] = as.numeric(paste0(merge[i,"output_lfcs"]))/as.numeric(paste0(merge[merge$experiment == "moi1_multi_moi" & merge$gene_pert == merge$gene_pert[i],"output_lfcs"]))
    scale_factor_per_mean_exp_moi1 = mean(as.numeric(paste0(merge[(merge$experiment == "moi1_lsg" | merge$experiment == "moi1_multi_moi")  & merge$gene_pert == merge$gene_pert[i],"output_lfcs"])))
    merge[i,"scaled_to_mean_exp_moi1"] = as.numeric(paste0(merge[i,"output_lfcs"]))/scale_factor_per_mean_exp_moi1
    merge[i,"ci_high_to_mean_exp_moi1"] =  as.numeric(paste0(merge[i,"ci_highs"]))/scale_factor_per_mean_exp_moi1 
    merge[i,"ci_low_to_mean_exp_moi1"] =  as.numeric(paste0(merge[i,"ci_lows"]))/scale_factor_per_mean_exp_moi1 
  }
}

tab2 = tab
tab2$experiment_use = tab2$experiment
for (i in 1:nrow(tab2)){
  if (tab2$experiment_use[i] == "moi1_multi_moi" | tab2$experiment_use[i] == "moi1_lsg"){
    tab2$experiment_use[i] = "moi1_combined_reps"
  }
}


exp_order2 = c("moi1_combined_reps","moi3_multi_moi","moi5_multi_moi", "moi6_multi_moi", "moi9_multi_moi", "moi14_multi_moi","moi20_hsg")

p_v2 = ggplot(tab2,aes(x = factor(experiment_use,levels = exp_order2), y = as.numeric(paste0(scaled_to_mean_exp_moi1)),fill = factor(experiment_use,levels = exp_order2)))+
  geom_boxplot()+
  scale_x_discrete(labels = c("moi1*","moi6","moi9","moi14","moi20"))+
  mytheme+
  geom_hline(yintercept = 1, linetype = "dashed")+
  scale_fill_manual(labels = c("moi1*","moi6","moi9","moi14","moi20"),values = brewer.pal(5,"YlOrRd"))+
  labs(fill = "MOI",y = "KD Normalized to mean MOI1 KD", x = "")+
  ylim(0,1.75)

png("figures/fig1_d_k562_multi_moi_use.png")
print(p_v2)
dev.off()

pdf("figures/fig1_d_k562_multi_moi_use.pdf")
print(p_v2)
dev.off()

#wtc11 multi_moi_s_fig
tab = read.table("figures/inputs_for_fig1/wtc11_all_moi_mast_tab.txt",header = T,sep = '\t')
tab_tss = tab[grepl(tab$perturbation,pattern =  "TSS"), ]
keep_rows = c()
for (i in 1:nrow(tab_tss)){
  gene = tab_tss[i,"gene"]
  if (grepl(tab_tss$perturbation[i],pattern = gene)){
    keep_rows = c(keep_rows,i)
  }
}

tss_tab1 = tab_tss[keep_rows,]
tss_tab1_use = tss_tab1[tss_tab1$pval_adj_moi1 < 0.01,]

wtc11_tab = melt(data =tss_tab1_use,id.vars = c("perturbation","gene","logFC_moi1"),measure.vars = c("logFC_moi1", "logFC_moi2","logFC_moi3",
                                                                                                     "logFC_moi6","logFC_moi10"))
wtc11_tab$normalized_to_moi1 = wtc11_tab$value / wtc11_tab$logFC_moi1



p_v3 = ggplot(wtc11_tab,aes(x = variable, y = as.numeric(paste0(normalized_to_moi1)),fill = variable))+
  geom_boxplot()+
  scale_x_discrete(labels = c("moi1","moi2","moi3","moi6","moi10"))+
  mytheme+
  geom_hline(yintercept = 1, linetype = "dashed")+
  scale_fill_manual(labels = c("moi1","moi2","moi3","moi6","moi10"),values =c("#FFFFB2","#FFEDA0","#FED976","#FECC5C","#FC4E2A"))+
  labs(fill = "MOI",y = "KD Normalized to mean MOI1 KD", x = "") +
  ylim(0,2)



png("figures/Fig_SX_WTC11_multi_moi_wtc11_multi_moi_use.png")
print(p_v3)
dev.off()

pdf("figures/Fig_SX_WTC11_multi_moi_wtc11_multi_moi_usewtc11_multi_moi_use.pdf")
print(p_v3)
dev.off()