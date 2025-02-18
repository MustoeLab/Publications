# Cryptic splicing in massively parallel reporter assays (MPRAs)

Publication: _U-rich elements drive pervasive cryptic splicing in 3'UTR massively parallel reporter assays_

Last updated: 07/30/24

## Analysis codes
  - _F0_process_raw_count_to_splicing_efficiency.ipynb_
      - Use categorized raw count (ptreseq_raw_count/) to compute splicing efficiency (ptreseq_splicing_quantification/)
      - Reported value is a median over biological replicates and internal barcode replicates
  - _F1_splicing_quantification.ipynb_
      - Codes to create figure 1
  - _F1S_splicing_reproducibility.ipynb_
      - Codes to create supplementary figure S1
  - _F2_donor_acceptor_identification.ipynb_
      - Codes to create figure 2 and supplementary figure S2
  - _F3_AU_rich_element_splicing.ipynb_
      - Codes to create figure 3 and supplementary figure S3
  - _F4_modelling_splicing_impact.ipynb_
      - Codes to create figure 4 (A-C) and supplementary figure S4
  - _F4_ptreseq_splicing_reanalysis.ipynb_
      - Codes to create figure 4D and supplementary figure S5
  - _F5_cryptic_splicing_mpra.ipynb_
      - Codes to create figure 5 and supplementary figure S6 and S7
     
## Supporting data files        
- **data**
    - Tape station data comparing DNA and RNA sequencing library. Related to Fig. 1B.
        - _2023-08-04 - 12-17-54-D1000_Electropherogram.csv_
    - Nucleotide sequences at boundaries of internal deletions in HeLa RNA sample. Related to Fig. 1D.
        - _HELA-1_gap_sequences.txt.gz_
    - Match between 9-nt barcode sequences and reporter name
        - _barcode_indexed.txt_
    - Single reporter validation with westernblot and RT-PCR. Related to Fig. 2 and Fig. S2
        - _single_reporter_validation_quantification.xlsx_
    - Maximum Entropy score of annotated human splice sites (Ensembl hg38 assembly). Related to Fig. 2
        - _human_5UTR_5ss_MAXENT.txt_
        - _human_5UTR_3ss_MAXENT.txt_
        - _human_CDS_5ss_MAXENT.txt_
        - _human_CDS_3ss_MAXENT.txt_
        - _human_3UTR_5ss_MAXENT.txt_
        - _human_3UTR_3ss_MAXENT.txt_
    - Measurements from original PTRE-seq publication
        - _sup2_readscount_plasmid.xlsx_ - read count for input DNA plasmid library
        - _sup3_readscount_HeLa_total.xlsx_ - read count for total RNA
        - _sup4_readscount_HeLa_polysome.xlsx_ - read count for polysome-associated RNA
    - Expression measurements from 3'UTR MPRA. Related to Fig. 5, Fig. S6, and Fig. S7.
        - _griesemer_expression_measurements.xlsx_
        - _siegel_expression_measurements_jurkat.csv_
        - _zhao_expression_measurements.xls_
    - SpliceAI prediction for 3'UTR MPRAs. Related to Fig. 5, Fig. S6, and Fig. S7.
        - _griesemer_spliceai_prediction.csv_
        - _siegel_spliceai_prediction.csv_
        - _zhao_spliceai_prediction.csv_
- **ptreseq_raw_count**
    - Raw count for RNA samples from different cell lines and biological replicates. Counts have been categorized using pipeline described in Fig. S3.
- **ptreseq_splicing_quantification**
    - Quantification of splicing efficiency using raw count


## Changelog
07/30/24
  - Intial commit: add analysis codes and data files to generate all figures for publication
