### use EDirect NCBI service to replace SNPID to SNP name for OAgwas file.
### this script only output snpOA.txt. which is quicker to obtain mapped refSNP id.

#!/bin/bash
echo "start  at "$(date)""

### remove header
#sed '1d' KP.Format.GO.FILTER.GW.AllOA.FULL.09052019.txt > temp.txt

### find how many rows in original OA gwas file. run outside of this script
rowN=$(cat temp.txt | wc -l)

for ((k=1; k<=(($rowN/1000)); k++))
do
  ### lines < 10,0000 1s/per snp; <100,0000 5s/per snp. Seperate files into small blocks to speed up, here we choose 1000 lines as a block. 
  sed -n $((1000*(k-1)+1)),$((1000*k))p temp.txt > OApre2.txt 

  ### extract Chromosome number and base position for each entry 
  for row in {1..1000}
  do
  thisCHR=$(awk -v thisRow=$row '{if(NR==thisRow) print $8}' OApre2.txt)
  thisPOS=$(awk -v thisRow=$row '{if(NR==thisRow) print $9}' OApre2.txt)

  ### API for snp database in NCBI
  searchResult=$(esearch -db snp -query '('$thisCHR'[Chromosome]) AND '$thisPOS'[Base Position Previous]' | efetch -format uid)

  if [[ -z "$searchResult" ]]  ###if not search found 
  then
  NewsnpOA="None"
  else
  NewsnpOA=${searchResult%%$'\n'*} ### if multiple returned searches
  fi

  ### concatenate retrieved refSNP for each row
  snpOA+="$row:${NewsnpOA} "  ### row number here is only used for error checking.

  echo "row $row finished"
  done

### write the output refSNP into file 
echo $snpOA | tr " " "," >> snpOA.txt

unset snpOA  ### snpOA variable should be limited to each small block

echo "$k round end"
done

### write the output refSNP into seprate lines and add "s" symbol
### append the refSNP column to original GWAS file. The final file OAGWAS.txt will be input for ldsc analysis.
cat snpOA.txt | tr " " ","  | sed -e 's/\,/\n/g' | sed 's/^.*:/rs/' > snpOAall.txt 
paste -d'\t' snpOAall.txt temp.txt > OAGWAS1.txt
   
### original header SNPID	EffectAllele	AlternateAllele	EffectAlleleFrequency	EffectSize.Beta	Pvalue	SampleSize	Chromosome	Position	EffectSize.OR
### append to OAGWAS2.txt
awk '{print $1,$3,$4,$5,$6,$7,$8,$9,$10}' OAGWAS1.txt >> OAGWAS2.txt 
 
### add colnames 
echo -e "SNP\tA1\tA2\tREF\tBeta\tP\tN\tCHR\tPOS" | cat -  OAGWAS2.txt > OAGWAS3.txt

### further check search results, delete records with failure search
awk '{if(NF!=9) print $0}' OAGWAS3.txt > OAGWAS4.txt
