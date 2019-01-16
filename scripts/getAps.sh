AP_DIR=$1
AP_PREFIX=$2
AP_LINE=$3
AP_OUT=$4

for i in `ls -1v ${AP_DIR}/${AP_PREFIX}*`; do sed -n "${AP_LINE}p" $i | awk -F ',' '{print $7}'  ; done > ${AP_OUT}.txt
