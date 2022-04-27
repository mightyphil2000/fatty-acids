cd /projects/MRC-IEU/users/ph14916/eqtl_bbj

gunzip *.gz

tar -xvf eQTL_B_cells.tar 
tar -xvf eQTL_CD4+T_cells.tar  
tar -xvf eQTL_CD8+T_cells.tar  
tar -xvf eQTL_Monocytes.tar  
tar -xvf eQTL_NK_cells.tar  
tar -xvf eQTL_Peripheral_blood.tar

cd /projects/MRC-IEU/users/ph14916/eqtl_bbj/CD8+T_cells
mkdir fattyacids 

mv chr11_cis_eqtl_mapping_nofilt_nomulti_with_alleles.txt.gz fattyacids/
mv chr6_cis_eqtl_mapping_nofilt_nomulti_with_alleles.txt.gz fattyacids/
mv chr10_cis_eqtl_mapping_nofilt_nomulti_with_alleles.txt.gz fattyacids/
mv chr2_cis_eqtl_mapping_nofilt_nomulti_with_alleles.txt.gz fattyacids/
mv chr16_cis_eqtl_mapping_nofilt_nomulti_with_alleles.txt.gz fattyacids/
mv chr20_cis_eqtl_mapping_nofilt_nomulti_with_alleles.txt.gz fattyacids/