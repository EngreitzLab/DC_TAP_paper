#List of all the guides in each screen, noting whether it is an experimental element, 
#positive or negative control, target region, and the coordinates of the 2MB in which is located.

setwd("/Users/ejagoda/Documents/HGRM/Crispri/dc_tap_paper_uploads/tables/inputs_for_Table_SX_RandomScreenDesign_2")

cell_type = "k562"
grna_table_sceptre = read.table("new_merged_k562_grna_groups_table_fixed.txt",header=T,sep = '\t')
cell_ranger_input = read.table("K562_random_screen_tap_guides_feature_reference_FIXED_USE.csv",header = T,sep = ',')


colnames(grna_table_sceptre)[2] = "name"
cell_ranger_input_small = cell_ranger_input[,c("name","sequence")]

merge1 = merge(grna_table_sceptre,cell_ranger_input_small)
colnames(merge1)[4] = "alt_name"

done_tab = merge1[,c("name","grna_target","gRNA_type","sequence","alt_name")]
write.table(done_tab,paste0(cell_type,"guide_tab.txt"),quote = F,sep = '\t',row.names = F)



cell_type = "wtc11"
grna_table_sceptre = read.table("new_merged_wtc11_grna_groups_table_fixed.txt",header=T,sep = '\t')
cell_ranger_input = read.table("wtc11_encode_feature_reference_file_dups_removed_FIXED.csv",header = T,sep = ',')


colnames(grna_table_sceptre)[2] = "name"
cell_ranger_input_small = cell_ranger_input[,c("name","sequence")]

merge1 = merge(grna_table_sceptre,cell_ranger_input_small)
colnames(merge1)[4] = "alt_name"

done_tab = merge1[,c("name","grna_target","gRNA_type","sequence","alt_name")]
done_tab$gRNA_type = ifelse(done_tab$gRNA_type == "DE","pos_control_enh",done_tab$gRNA_type)
write.table(done_tab,paste0(cell_type,"guide_tab.txt"),quote = F,sep = '\t',row.names = F)

