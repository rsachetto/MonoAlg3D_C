# =======================================================================================
# Bash script to convert a set of frames stored in .png to video animation in .mp4
# Author: Lucas Berg
# =======================================================================================
#!/bin/bash

# Variables
FILENAME="frames/frame"
FRAME_RATE="5"
START_FRAME="10"
END_FRAME="40"
OUTPUT_VIDEO_FILENAME="videos/oxford_DTI003_only_root_nodes_simulation_comparison_with_ecg_leads"
RESOLUTION="1408x738"

# Execute the converting command using FFMPEG
ffmpeg -r ${FRAME_RATE} -f image2 -s ${RESOLUTION} -start_number ${START_FRAME} -i ${FILENAME}.%04d.png -vframes ${END_FRAME} -vcodec libx264 -crf 25  -pix_fmt yuv420p ${OUTPUT_VIDEO_FILENAME}.mp4
