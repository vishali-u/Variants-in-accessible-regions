#!/bin/bash

trait='sz'
folder_name='BCB330/gwas_dar/sz/sz_intersect'

for ((i=1;i<=20585;i+=5));
do 
  echo "#!/bin/bash" > "/nethome/kcni/vumaiyalan/${folder_name}/${trait}_${i}.sh"
  echo "Rscript /nethome/kcni/vumaiyalan/BCB330/gwas_dar/sz/intersect.R ${i}" >> "/nethome/kcni/vumaiyalan/${folder_name}/${trait}_${i}.sh"
  echo "/nethome/kcni/vumaiyalan/${folder_name}/${trait}_${i}.sh" >> "/nethome/kcni/vumaiyalan/${folder_name}/script_list.txt" 
  chmod +x /nethome/kcni/vumaiyalan/${folder_name}/${trait}_${i}.sh
  echo "File ${trait}_${i}.sh generated"
done
