# Key Paths

## Cell Ranger Outputs:
- **K562**: /oak/stanford/groups/engreitz/Users/ejagoda/230327_Encode_K562_Tap_seq_full_seq/230327_Encode_Tap_all_moi5_w_fresh_aggregate_no_norm
- **WTC11**: /oak/stanford/groups/engreitz/Users/ejagoda/231114-WTC11-ENCODE-DC-TAP-MOI5-Frozen-16-Lanes/count/

### Genome Build to Use:
- v32 (v32lift37) when needed for hg19-based analyses — original target/guide code is hg19.

## Differential Expression Outputs:
- **K562**: /oak/stanford/groups/engreitz/Users/jgalante/ENCODE_Sceptre_Analysis/results/process_validation_datasets/K562_DC_TAP_Seq
- **WTC11**: /oak/stanford/groups/engreitz/Users/jgalante/ENCODE_Sceptre_Analysis/results/process_validation_datasets/WTC11_DC_TAP_Seq

## Pilot Studies:
- **1. TAP-seq pilot assay (JR06_21-11-03)**:
  - Direct capture 10x output: /oak/stanford/groups/engreitz/Users/kangh/process_sequencing_data/220210_CRISPRi_direct_capture/CRISPRiGuideCaptureTAP-seqUpdatedFeatureReference
  - Path to TAPseq_workflow to generate dge.txt for all and guide assignment for the TAP-seq (crop-seq vector): /oak/stanford/groups/engreitz/Users/ejagoda/dc_tap_seq/from_Helen_scratch/TAPseq_workflow/
  -  DE files:
      - /oak/stanford/groups/engreitz/Users/ejagoda/reorganization_2025/crispri_analysis_from_evvies_local/JR06_21-11-03_pilot_tap_seq_ampli_dc/v3/fromHelen_P10/summary_table_P10_tap_kd_mean_fc_all_genes.txt
      - /oak/stanford/groups/engreitz/Users/ejagoda/reorganization_2025/crispri_analysis_from_evvies_local/JR06_21-11-03_pilot_tap_seq_ampli_dc/v3/summary_table_kd_mean_fc_all_ampli_tap_samples_all_genes.txt
  - Analysis scripts:
      - /oak/stanford/groups/engreitz/Users/ejagoda/reorganization_2025/crispri_analysis_from_evvies_local/JR06_21-11-03_pilot_tap_seq_ampli_dc/v3/fromHelen_P10/P10_tap_analysis.R
      -  /oak/stanford/groups/engreitz/Users/ejagoda/reorganization_2025/crispri_analysis_from_evvies_local/JR06_21-11-03_pilot_tap_seq_ampli_dc/v3/from_Helen_ampliTap/ampli_tap_analysis.R
      -  /oak/stanford/groups/engreitz/Users/ejagoda/reorganization_2025/crispri_analysis_from_evvies_local/JR06_21-11-03_pilot_tap_seq_ampli_dc/v3/fromHelen_DC/Analyze_gene_expression_data_dc_10x.R

- **2. TAP-direct capture test (Low vs High MOI) (220523_high_moi_tap_dc_ampli_tap_pilot) (JR18_22-05-13)**:
  - Path to 10x results: /oak/stanford/groups/engreitz/Users/ejagoda/hgrm/220523_high_moi_tap_dc_ampli_tap_pilot/10x_analysis/
  - Path to differential expression outputs: /oak/stanford/groups/engreitz/Users/ejagoda/reorganization_2025/crispri_analysis_from_evvies_local/JR18_22-05-13_220523_high_moi_tap_dc_ampli_tap_pilot/220523_high_moi_tap_dc_ampli_tap_pilot/10x_analysis/052322_analysis_{FULL_SAMPLE_NAME}/{short_sample_name}_w_all_pertsdifferential_expression_output_tab_WITHLOCUSFILTERS_fixed.tx

- **3. K562 multi moi  (JR28_22-09-13)**
  - Path to 10x results: /oak/stanford/groups/engreitz/Users/ejagoda/multi_moi/10x_analysis
  - Path to  differential expression results: /oak/stanford/groups/engreitz/Users/ejagoda/Encode_K562_Random_TAP_screen/CRISPRiScreen-example/results/multi_moi{1,3,5,6,9,14}
  - Post MAST processing script: /oak/stanford/groups/engreitz/Users/ejagoda/reorganization_2025/crispri_analysis_from_evvies_local/JR28_22-09-13_k562_multi_moi/multi_moi/analyze_mast_multi_moi.R

- **4. wtc11 multi moi (JR38_23-02-10)**
  - Path to 10x results: /oak/stanford/groups/engreitz/Users/ejagoda/230309_WTC11_Multi_MOI_TAP_DC_Pilot_FIXED/{moi}_230309_WTC11_Multi_MOI_TAP_DC_Pilot_new_seq
  - Path to differential expression results: /oak/stanford/groups/engreitz/Users/ejagoda/Encode_K562_Random_TAP_screen/CRISPRiScreen-example/results/230309_WTC11_Multi_MOI_TAP_DC_Pilot_FIXED_{moi}
  - Post MAST processing script: oak/stanford/groups/engreitz/Users/ejagoda/reorganization_2025/crispri_analysis_from_evvies_local/JR38_23-02-10_230309_WTC11_Multi_MOI_TAP_DC_Pilot/230309_WTC11_Multi_MOI_TAP_DC_Pilot]$ more analyze_snakemake_mast_230309_WTC11_Multi_MOI_TAP_DC_Pilot.R
    
- **Shared Scripts for processing 2,3,4**
  - Make dge file: /oak/stanford/groups/engreitz/Users/ejagoda/230327_Encode_K562_Tap_seq_full_seq/process_gene_expression_and_make_dge_file_230327_Encode_K562_Tap_seq_full_seq_for_server_w_untar.sh
  - Make perturb_status_file:
  - 1. /oak/stanford/groups/engreitz/Users/ejagoda/Encode_K562_Random_TAP_screen/scripts/crispri_get_perturbation_status_faster.sh
  - 2. /oak/stanford/groups/engreitz/Users/ejagoda/Encode_K562_Random_TAP_screen/scripts/transpose_attempt.R 
- Differential expression analysis for experiments 2,3,4 done via MAST here: /oak/stanford/groups/engreitz/Users/ejagoda/Encode_K562_Random_TAP_screen/CRISPRiScreen-example/
 
