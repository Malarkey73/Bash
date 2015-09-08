-bash-4.1$ cat table_annovar.sh
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
# nohup ./table_annovar.sh S1.anno

# alternatie parallel usage example
# ls *.pass | parallel table_annovar.sh {}


#tools
TABLE_ANNOVAR="/home/rmgzshd/annovar/table_annovar.pl"
#annotations
DB="/home/rmgzshd/annovar/humandb"
#inputs
FILE=`echo $1`
PREFIX=`echo "${FILE%%.*}"`

export TABLE_ANNOVAR; export DB; export FILE; export PREFIX


#paste <(cut -f1,2 $FILE) <(cut -f2,4,5 $FILE) > $PREFIX.FILE

$TABLE_ANNOVAR $FILE $DB \
-buildver hg19 \
-out $PREFIX \
-protocol refGene,cytoBand,genomicSuperDups,1000g2014oct_all,snp138,cosmic70 \
-operation g,r,r,f,f,f \
-nastring .

