# U-rich elements drive pervasive cryptic splicing in 3’ UTR massively parallel reporter assays 

Last updated: 02/18/2025


### Analysis codes
  - `F0_process_raw_count_to_splicing_efficiency.ipynb`
      - `Input`: raw count by categories for each reporters (ptreseq_raw_count/). See **Fig. S1D** and **Method** for a description of each category. Codes for generating raw count are available upon request.
      - `Output`: median observed splicing fraction (ptreseq_splicing_quantification/)
  - `F1_splicing_quantification.ipynb`: Figure 1
  - `F1S_splicing_reproducibility.ipynb`: Figure S1
  - `F2_donor_acceptor_identification.ipynb`: Figure 2 and S2
  - `F3_AU_rich_element_splicing.ipynb`: Figure 3 and S3
  - `F4_modelling_splicing_impact.ipynb`: Figure 4 (A-C) and S4
  - `F4_ptreseq_splicing_reanalysis.ipynb`: Figure 4D and S5
  - `F5_cryptic_splicing_mpra.ipynb`: Figure 5, S6, and S7
     
### additional data 
  - Tape station data comparing DNA and RNA sequencing library
        - _2023-08-04 - 12-17-54-D1000_Electropherogram.csv_
  - Nucleotide sequences at boundaries of internal deletions in HeLa RNA sample
        - _HELA-1_gap_sequences.txt.gz_
  - Match between 9-nt barcode sequences and reporter name
      - _barcode_indexed.txt_
  - Single reporter validation with westernblot and RT-PCR
      - _single_reporter_validation_quantification.xlsx_
  - Maximum Entropy score of annotated human splice sites (Ensembl hg38 assembly)
      - _human_5UTR_5ss_MAXENT.txt_
      - _human_5UTR_3ss_MAXENT.txt_
      - _human_CDS_5ss_MAXENT.txt_
      - _human_CDS_3ss_MAXENT.txt_
      - _human_3UTR_5ss_MAXENT.txt_
      - _human_3UTR_3ss_MAXENT.txt_
      - _human_5ss_sequence.fa_ and _human_3ss_sequence.fa_ (sequence for all annotated human splice sites)
  - Measurements from original PTRE-seq publication
      - _sup2_readscount_plasmid.xlsx_ - read count for input DNA plasmid library
      - _sup3_readscount_HeLa_total.xlsx_ - read count for total RNA
      - _sup4_readscount_HeLa_polysome.xlsx_ - read count for polysome-associated RNA
  - Expression measurements from 3'UTR MPRA
      - _griesemer_expression_measurements.xlsx_
      - _siegel_expression_measurements_jurkat.csv_
      - _zhao_expression_measurements.xls_
      - _Fu_MPRA_expression_measurements.xlsx_
  - SpliceAI prediction for 3'UTR MPRAs
      - _griesemer_spliceai_full_transcript_prediction_best_sites.txt_
      - _zhao_spliceai_full_transcript_prediction_best_sites.txt_
      - _siegel_spliceai_full_transcript_prediction_best_sites.txt_
  - XSTREME motif analysis for spliced reporters in 3'UTR MPRAs
      - _griesemer_motif_analysis.txt_ and _griesemer_motif_analysis.html_
      - _siegel_motif_analysis.txt_ and _siegel_motif_analysis.html_
      - _zhao_motif_analysis.txt_ and _zhao_motif_analysis.html_
- **PLOTS**: contain all figures saved as PDF


## Changelog
07/30/24
  - Intial commit: add analysis codes and data files to generate all figures for publication
02/18/25
  - Update code following revision
