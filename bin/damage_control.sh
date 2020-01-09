#/bin/bash

test=$(ls ../fastq/ | awk -F '[-_]' '{print $1"-"$2}' | uniq)

IFS=' ' read -r -a array <<< $test
printf '%s\n' "${array[@]}"

#ls ../fastq/ | awk -F "[_-]" '{print $1"-"$2}' | uniq | parallel --jobs 1 du -h ../fastq/{.}*




#ls ../fastq/ | awk -F "[_-]" '{print $1"-"$2}' | uniq | parallel --jobs 11 bash ../bin/damage_pe_script.sh {.} Artemisia /usr/share/data/ncbi/genomes/plastids/plastid_all.fna
#ls ../fastq/ | awk -F "[_-]" '{print $1"-"$2}' | uniq | parallel --jobs 11 bash ../bin/damage_pe_script.sh {.} Lupinus /usr/share/data/ncbi/genomes/plastids/plastid_all.fna
 
#ls ../fastq/ | awk -F "[_-]" '{print $1"-"$2}' | uniq | parallel --jobs 11 bash ../bin/damage_pe_script.sh {.} Pinus /usr/share/data/ncbi/genomes/plastids/plastid_all.fna
#ls ../fastq/ | awk -F "[_-]" '{print $1"-"$2}' | uniq | parallel --jobs 11 bash ../bin/damage_pe_script.sh {.} Malus /usr/share/data/ncbi/genomes/plastids/plastid_all.fna
#ls ../fastq/ | awk -F "[_-]" '{print $1"-"$2}' | uniq | parallel --jobs 11 bash ../bin/damage_pe_script.sh {.} Juniperus /usr/share/data/ncbi/genomes/plastids/plastid_all.fna

#ls ../fastq/ | awk -F "[_-]" '{print $1"-"$2}' | uniq | parallel --jobs 11 bash ../bin/damage_pe_script.sh {.} "Homo " /usr/share/data/ncbi/genomes/mito/mito_all.fna
#ls ../fastq/ | awk -F "[_-]" '{print $1"-"$2}' | uniq | parallel --jobs 11 bash ../bin/damage_pe_script.sh {.} Neotoma /usr/share/data/ncbi/genomes/mito/mito_all.fna
 
#ls ../fastq/ | awk -F "[_-]" '{print $1"-"$2}' | uniq | parallel --jobs 11 bash ../bin/damage_pe_script.sh {.} Plastid /usr/share/data/ncbi/genomes/plastids/plastid_all.fna
#ls ../fastq/ | awk -F "[_-]" '{print $1"-"$2}' | uniq | parallel --jobs 11 bash ../bin/damage_pe_script.sh {.} Mito /usr/share/data/ncbi/genomes/mito/mito_all.fna
 
ls ../results/ | uniq | parallel --jobs 11 bash ../bin/damage_pmdtools.sh {.} Plastid /usr/share/data/ncbi/genomes/plastids/plastid_all.fna

exit 1


#printf '%s\n' "${array[@]}" | parallel --jobs 11 bash ../bin/damage_pe_script.sh ../fastq/{.}* 