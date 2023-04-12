#!/bin/bash

trait='sz'
folder_name='shell_scripts'
file_path='/nethome/kcni/vumaiyalan/BCB330/gwas_dar/sz_trevino'

for ((i=1;i<=20585;i+=5));
do 
  echo "#!/bin/bash" > "${file_path}/${folder_name}/${trait}_${i}.sh"
  echo "Rscript ${file_path}/intersect.R ${i}" >> "${file_path}/${folder_name}/${trait}_${i}.sh"
  echo "${file_path}/${folder_name}/${trait}_${i}.sh" >> "${file_path}/${folder_name}/script_list.txt" 
  chmod +x ${file_path}/${folder_name}/${trait}_${i}.sh
  echo "File ${trait}_${i}.sh generated"
done
