# first run extract_regions_eqtlgen.sh
cd /projects/MRC-IEU/users/ph14916/eQTLGen
head -1 eqtlgen_ELOVL2.txt > head.txt
grep -v -h cis.trans eqtlgen_*.txt | cat head.txt - > eqtlgen_data.txt

cd ~/fatty-acids/colocalisation/data/
scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/eQTLGen/eqtlgen_data.txt .
 
# head eqtlgen_data.txt
# wc eqtlgen_data.txt

# way to do it in one line:
# head -1 eqtlgen_ELOVL2.txt >> grep -v -h cis.trans eqtlgen_*.txt > eqtlgen_data.txt