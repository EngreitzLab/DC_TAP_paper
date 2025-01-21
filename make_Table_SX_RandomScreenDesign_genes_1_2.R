##for the gene tables
#Table SX_RandomScreenDesign_1
library(stringr)
library(biomaRt)



master_path =  "/Users/ejagoda/Documents/HGRM/Crispri/dc_tap_paper_uploads/tables/inputs_for_Table_SX_RandomScreenDesign_1"

setwd(master_path)

cell_type = "k562"
design_output_file_name = "K562.rep6.GenesForTAPseq_wV29_genenames_controls_fixed.txt"
primer_output_list_o_name = "K562-ENCODE-Primers.xlsx\ -\ outer_primers_replaced.tsv"
primer_output_list_i_name = "K562-ENCODE-Primers.xlsx\ -\ inner_primers_replaced.tsv"
locus_list_tab = read.table("rep_6_loci_with_3_other_50tpm_genes_w_g_40_w_g_30.txt",header=T,sep = '\t')
locus_coordinates = read.table("loci_only_bed_version_rep6.bed")
colnames(locus_coordinates) = c("chr","start","end","locus")
locus_tab_all = merge(locus_coordinates,locus_list_name)
outer_controls = read.table("control_primer_info_outer.txt",header = T,sep = '\t')
tpm_file = read.table("k562_tpm.csv",sep = ",",header=T)

tab = read.table(design_output_file_name,header=T,sep = '\t')
primer_output_list_o = read.table(primer_output_list_o_name,header = T,sep = '\t')
primer_output_list_i = read.table(primer_output_list_i_name,header = T,sep = '\t')

tab_use = data.frame(
  gene = character(),
  ensemble_id = character(),
  tpm = numeric(),
  type = character(),
  locus = character(),
  outer_primer = character(),
  inner_primer = character()
)

for (i in 1:nrow(primer_output_list_o)){
  name = str_split(primer_output_list_o$primer_id[i],pattern = "\\.")[[1]][1]
  primer_output_list_o$gene_name[i] = name
  g_tab = tab[tab$gene_name_v29 == name,]
  locus_tab = locus_tab_all[grepl(locus_tab_all$genes_g_30tpm,pattern = name),]
  inner_primer_tab = primer_output_list_i[grepl(primer_output_list_i$primer_id,pattern = name),]
  inner_primer = inner_primer_tab[1,"sequence"]
  if (nrow(g_tab) == 1){
    ensemble_id = g_tab[1,"gene"]
    tpm = as.numeric(paste0(g_tab[1,"tpm"]))
    type = g_tab[1,"type"]
  }else{
    control_tab = outer_controls[grepl(outer_controls$gene,pattern = name),]
    if (nrow(control_tab) > 0){
      type = paste0(str_split(control_tab$type[1],pattern = "-")[[1]][1],"_control")
      ensemble_id = "x"
      tpm = "x"
    }else{
      type = "x"
      ensemble_id = "x"
      tpm = "x"
    }
  }
  if (nrow(locus_tab) == 1){
    locus = paste0(locus_tab[1,"chr"],":",locus_tab[1,"start"],"-",locus_tab[1,"end"])
  }else{
    locus = "x"
  }
  tab_use = data.frame(rbind(tab_use,c(name,ensemble_id,tpm,type,locus,primer_output_list_o$sequence[i],inner_primer)))
}

colnames(tab_use) = c("gene","ensemble_id","tpm", "type","locus","outer_primer","inner_primer")

#double checked that the two with the missing type are locus 
for (i in 1:nrow(tab_use)){
  if (tab_use$tpm[i] == "x"){
  gene = tab_use$gene[i]
  tpm_tab = tpm_file[tpm_file$gene_name == gene,]
  tab_use$ensemble_id[i] = tpm_tab$gene[1]
  tab_use$tpm[i] = tpm_tab$tpm[1]
  }
  if (tab_use$type[i] == "x" & tab_use$locus[i] != "x"){
    tab_use$type[i] = "locus_gene"
  }
}
#figure out why those 2 locus genes don't ahve tpm and you'll be fineee

write.table(tab_use,paste0(cell_type,"_gene_table.txt"),quote = F,sep = '\t',row.names = F)



#################

cell_type = "wtc11"
design_output_file_name = "wtc11_full_gene_set_use_real_w_ensemble_id.txt"
primer_output_list_o_name = "WTC11-ENCODE_primer_list-ordered - Outer.tsv"
primer_output_list_i_name = "WTC11-ENCODE_primer_list-ordered - Inner.tsv"
locus_list_tab = read.table("rep_106_loci_with_3_other_50tpm_genes_w_g_40_w_g_30.txt",header=T,sep = '\t')
locus_tab_all = locus_list_tab
outer_controls = read.table("control_primer_info_outer.txt",header = T,sep = '\t')
tpm_file = read.table("wtc11_tpm.out",sep = "\t",header=T,fill = T)


tab = read.table(design_output_file_name,header=T,sep = '\t')
primer_output_list_o = read.table(primer_output_list_o_name,header = T,sep = '\t')
primer_output_list_i = read.table(primer_output_list_i_name,header = T,sep = '\t')

tab_use = data.frame(
  gene = character(),
  alt_name = character(),
  ensemble_id = character(),
  tpm = numeric(),
  type = character(),
  locus = character(),
  outer_primer = character(),
  inner_primer = character()
)

for (i in 1:nrow(primer_output_list_o)){
  name = primer_output_list_o$Name[i]
  primer_output_list_o$gene_name[i] = name
  g_tab = tab[tab$v29_gene_name == name | tab$v26_gene_name == name,]
  v26_name = g_tab$v26_gene_name[1] 
  v29_name = g_tab$v29_gene_name[1]
  locus_tab = locus_tab_all[grepl(locus_tab_all$genes_g_30tpm,pattern = v26_name) | grepl(locus_tab_all$genes_g_30tpm,pattern = v29_name),]
  inner_primer_tab = primer_output_list_i[primer_output_list_i$Name == name,]
  inner_primer = inner_primer_tab[1,"Sequence"]
  if (nrow(g_tab) == 1){
    ensemble_id = g_tab[1,"ensemble_id"]
    tpm = as.numeric(paste0(g_tab[1,"tpm"]))
    type = g_tab[1,"type"]
    gene_names = unique(c(name, v26_name,v29_name))
    if (length(gene_names) == 1){
      alt_name = "x"
    }else{
      alt_name = gene_names[which(gene_names != name)]
    }
  }else{
    control_tab = outer_controls[grepl(outer_controls$gene,pattern = name),]
    alt_name ="x"
    if (nrow(control_tab) > 0){
      type = paste0(str_split(control_tab$type[1],pattern = "-")[[1]][1],"_control")
      ensemble_id = "x"
      tpm = "x"
    }else{
      type = "x"
      ensemble_id = "x"
      tpm = "x"
    }
  }
  if (nrow(locus_tab) == 1){
    locus = paste0(locus_tab[1,"chr"],":",locus_tab[1,"start"],"-",locus_tab[1,"end"])
  }else{
    locus = "x"
  }

  tab_use = data.frame(rbind(tab_use,c(name,alt_name,ensemble_id,tpm,type,locus,primer_output_list_o$Sequence[i],inner_primer)))
}

colnames(tab_use) = c("gene","alt_name", "ensemble_id","tpm", "type","locus","outer_primer","inner_primer")

#double checked that the two with the missing type are locus 
for (i in 1:nrow(tab_use)){
  if (tab_use$tpm[i] == "x"){
    gene = tab_use$gene[i]
    alt = tab_use$alt_name[i]
    tpm_tab = tpm_file[tpm_file$Gene_Id == gene | tpm_file$Gene_Id == alt_name,]
    #tab_use$ensemble_id[i] = tpm_tab$gene[1]
    tab_use$tpm[i] = tpm_tab$TPM[1]
  }
  if (tab_use$type[i] == "x" & tab_use$locus[i] != "x"){
    tab_use$type[i] = "locus_gene"
  }
  if (gene == "TFC4"){
    tab_use$ensemble_id[i] = "ENSG00000196628"
    tab_use$type[i] = "control"
  }else if (gene == "YBX1"){
    tab_use$ensemble_id[i] = "ENSG00000065978"
    tab_use$type[i] = "control"
  }
}
#figure out why those 2 locus genes don't ahve tpm and you'll be fineee

#TCF4 and YBX1 are controls ?
#ENSG00000196628
#ENSG00000065978

write.table(tab_use,paste0(cell_type,"_gene_table.txt"),quote = F,sep = '\t',row.names = F)





#Figure out what needs to be uploaded for this

