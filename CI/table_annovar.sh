#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail


# these are the script lines used to download annotations to the humandb folder
# annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
# annotate_variation.pl -buildver hg19 -downdb cytoBand humandb/
# annotate_variation.pl -buildver hg19 -downdb genomicSuperDups humandb/ 
# annotate_variation.pl -buildver hg19 -downdb -webfrom annovar esp6500siv2_all humandb/
# annotate_variation.pl -buildver hg19 -downdb -webfrom annovar 1000g2014oct humandb/
# annotate_variation.pl -buildver hg19 -downdb -webfrom annovar snp138 humandb/ 
# annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ljb26_all humandb/

# usage example
# nohup ./table_annovar.sh NB11.snp.Somatic.pass

# alternatie parallel usage example
# ls *.pass | nohup parallel ./table_annovar.sh {}


#tools
TABLE_ANNOVAR="/home/rmgzshd/annovar/table_annovar.pl"
#annotations
DB="/home/rmgzshd/annovar/humandb"
#inputs
FILE=`echo $1`
PREFIX=`echo "${FILE%%.*}"`

export TABLE_ANNOVAR; export DB; export FILE; export PREFIX


paste <(cut -f1,2 $FILE) <(cut -f2,3,4 $FILE) > $PREFIX.FILE

$TABLE_ANNOVAR $PREFIX.FILE \
$DB \
-buildver hg19 \
-out $PREFIX \
-remove \
-protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,snp138,ljb26_all \
-operation g,r,r,f,f,f,f,f,f,f \
-nastring NA

rm $PREFIX.FILE
