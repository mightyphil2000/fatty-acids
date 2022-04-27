#===============================================
# README for cis-QTLs for the v8 GTEx release
#===============================================

This document describes the output format and column headers
for cis-eQTL and cis-sQTL files from the v8 GTEx data release.

#-----------------------------------------------
# eGenes
#-----------------------------------------------

  File extension: *.egenes.txt.gz, *.sgenes.txt.gz
  Column headers:

    gene_id:                  GENCODE/Ensembl gene ID
    gene_name:                GENCODE gene name
    gene_chr:                 chromosome (gene)
    gene_start:               gene start position (in base pairs; 1-based coordinates)
    gene_end:                 gene end position (in base pairs; 1-based coordinates)
    strand:                   genomic strand
    num_var:                  number of variants in cis-window
    beta_shape1:              1st shape parameter of the fitted Beta distribution: B(shape1, shape2)
    beta_shape2:              2nd shape parameter of the fitted Beta distribution: B(shape1, shape2)
    true_df:                  Effective degrees of freedom the Beta distribution approximation
    variant_id:               variant ID in the format {chr}_{pos_first_ref_base}_{ref_seq}_{alt_seq}_b38
    tss_distance:             distance between variant and transcription start site (TSS). Positive when variant is downstream of the TSS, negative otherwise
    chr:                      chromosome (for the variant, same as gene_chr for cis-eQTLs)
    variant_pos:              position of the first reference base of the variant
    ref:                      reference sequence of the variant
    alt:                      alternate sequence of the variant
    num_alt_per_site:         number of alternative alleles observed at this site
    rs_id_dbSNP151_GRCh38p7:  dbSNP151 rsID
    minor_allele_samples:     number of samples carrying the minor allele
    minor_allele_count:       total number of minor alleles across individuals
    maf:                      minor allele frequency observed in the set of donors for a given tissue
    ref_factor:               '1', when the minor allele is the alt base, '-1' when the minor allele is the reference base
    pval_nominal:             nominal p-value associated with the most significant variant for this gene
    slope:                    regression slope
    slope_se:                 standard error of the regression slope
    pval_perm:                permutation p-value
    pval_beta:                beta-approximated permutation p-value
    qval:                     Storey q-value derived from pval_beta
    pval_nominal_threshold:   nominal p-value threshold for calling a variant-gene pair significant for the gene
    log2_aFC:                 allelic fold change (aFC) in log2 scale (see Mohammadi et al., Genome Res. 2017)
    log2_aFC_lower:           lower bound of the 95% confidence interval of log2(aFC)
    log2_aFC_upper:           upper bound of the 95% confidence interval of log2(aFC)

#-----------------------------------------------
# Significant variant-gene pairs
#-----------------------------------------------

  File extensions: *.signif_variant_gene_pairs.txt.gz, *.sqtl_signifpairs.txt.gz
  Column headers:

    variant_id:               variant ID in the format {chr}_{pos_first_ref_base}_{ref_seq}_{alt_seq}_b38
    gene_id:                  GENCODE/Ensembl gene ID
    tss_distance:             distance between variant and transcription start site. Positive when variant is downstream of the TSS, negative otherwise
    ma_samples:               number of samples carrying the minor allele
    ma_count:                 total number of minor alleles across individuals
    maf:                      minor allele frequency observed in the set of donors for a given tissue
    pval_nominal:             nominal p-value
    slope:                    regression slope
    slope_se:                 standard error of the regression slope
    pval_nominal_threshold:   nominal p-value threshold for calling a variant-gene pair significant for the gene
    min_pval_nominal:         smallest nominal p-value for the gene
    pval_beta:                beta-approximated permutation p-value for the gene


#-----------------------------------------------
# All variant-gene pairs tested
#-----------------------------------------------

  File extension: *.allpairs.txt.gz, *.sqtl_allpairs.txt.gz
  Column headers:

    gene_id:                  GENCODE/Ensembl gene ID
    variant_id:               variant ID in the format {chr}_{pos_first_ref_base}_{ref_seq}_{alt_seq}_b38
    tss_distance:             distance between variant and transcription start site. Positive when variant is downstream of the TSS, negative otherwise
    ma_samples:               number of samples carrying the minor allele
    ma_count:                 total number of minor alleles across individuals
    pval_nominal:             nominal p-value
    slope:                    regression slope
    slope_se:                 standard error of the regression slope
