AP_DIR=$1
AP_N=$2

for i in `ls -1v ${AP_DIR}/V_t_*`; do sed -n "${AP_N}p" $i | awk -F ',' '{print $5}'  ; done > ${3}.txt
