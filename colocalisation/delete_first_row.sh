
# cat /projects/MRC-IEU/users/ph14916/gtex2/all/GTex.header.txt Esophagus_Gastroesophageal_Junction.allpairs.txt.gz_ELOVL2?.txt  > temp.txt

sed '1d' Esophagus_Mucosa.allpairs.txt.gz_FADS?.txt > temp1 
cat /projects/MRC-IEU/users/ph14916/gtex2/all/GTex.header.txt temp1 > Esophagus_Mucosa.allpairs.txt.gz_FADS.txt

sed '1d' Esophagus_Muscularis.allpairs.txt.gz_FADS?.txt > temp2
cat /projects/MRC-IEU/users/ph14916/gtex2/all/GTex.header.txt temp2 > Esophagus_Muscularis.allpairs.txt.gz_FADS.txt

sed '1d' Esophagus_Gastroesophageal_Junction.allpairs.txt.gz_FADS?.txt > temp3
cat /projects/MRC-IEU/users/ph14916/gtex2/all/GTex.header.txt temp3 > Esophagus_Gastroesophageal_Junction.allpairs.txt.gz_FADS.txt

sed '1d' Lung.allpairs.txt.gz_FADS?.txt > temp4
cat /projects/MRC-IEU/users/ph14916/gtex2/all/GTex.header.txt temp4 > Lung.allpairs.txt.gz_FADS.txt 

sed '1d' Skin_Not_Sun_Exposed_Suprapubic.allpairs.txt.gz_FADS?.txt > temp5
cat /projects/MRC-IEU/users/ph14916/gtex2/all/GTex.header.txt temp5 > Skin_Not_Sun_Exposed_Suprapubic.allpairs.txt.gz_FADS.txt

sed '1d' Skin_Sun_Exposed_Lower_leg.allpairs.txt.gz_FADS?.txt > temp6
cat /projects/MRC-IEU/users/ph14916/gtex2/all/GTex.header.txt temp6 > Skin_Sun_Exposed_Lower_leg.allpairs.txt.gz_FADS.txt


rm temp*

wc Esophagus_Mucosa.allpairs.txt.gz_FADS?.txt 
wc Esophagus_Mucosa.allpairs.txt.gz_FADS.txt
wc Esophagus_Muscularis.allpairs.txt.gz_FADS?.txt 
wc Esophagus_Muscularis.allpairs.txt.gz_FADS.txt
wc Esophagus_Gastroesophageal_Junction.allpairs.txt.gz_FADS?.txt 
wc Esophagus_Gastroesophageal_Junction.allpairs.txt.gz_FADS.txt



