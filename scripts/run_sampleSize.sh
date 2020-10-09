#!/usr/bin/bash

### shell 脚本实现细胞数分析流程

#sed -i "s/sub20w/sub20w/g" get_sample.R num_PC_8_level.R num_PC_allCell.R num_clusters.R 

cd /home/chenhy/mouse_brain/sub_20w_rep/scripts/

sed -i "s/5678/1234/g" 04_cluster_singler_anno.R
rawseed=1234
newseed=(3456 6789 8901 9023 5678)
n=0
while [ $n -lt 5 ]
do 
	echo "change pre seed $rawseed with new seed ${newseed[n]}"
    sed -i "s/$rawseed/${newseed[n]}/g" 04_cluster_singler_anno.R 02_num_PC_8_level.R 01_get_sample.R
    rawseed=${newseed[n]} 
 #    Rscript 01_get_sample.R
#     wait
 #    Rscript 02_num_PC_8_level.R
	# wait
	Rscript 04_cluster_singler_anno.R
	wait
    let n+=1
    echo $n
done


# Rscript /home/chenhy/mouse_brain/scripts/num_PC_8_level.R
# wait

# Rscript /home/chenhy/mouse_brain/scripts/num_PC_allCell.R
# wait

# Rscript /home/chenhy/mouse_brain/scripts/num_clusters.R
# wait


