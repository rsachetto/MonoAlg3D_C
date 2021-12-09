#!/bin/bash

# Create File List
#echo file videos/shocker_rafa_sebastian_LV_iter_0.mp4 >  mylist.txt
#echo file videos/shocker_rafa_sebastian_LV_iter_1.mp4 >> mylist.txt
#echo file videos/shocker_rafa_sebastian_LV_iter_2.mp4 >> mylist.txt
#echo file videos/shocker_rafa_sebastian_LV_iter_3.mp4 >> mylist.txt

# Concatenate Files
ffmpeg -f concat -i mylist.txt -c copy output.mp4

# Clean Files
rm mylist.txt
