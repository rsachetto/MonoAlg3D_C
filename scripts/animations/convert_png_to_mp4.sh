# =======================================================================================
# Bash script to convert a set of frames stored in .png to video animation in .mp4
# Author: Lucas Berg
# =======================================================================================
#!/bin/bash

# Variables
FILENAME="frames/frame"
FRAME_RATE="10"
END_FRAME="200"
OUTPUT_VIDEO_FILENAME="video/ddm_3d"
RESOLUTION="1020x720"

# Execute the converting command using FFMPEG
ffmpeg -r $FRAME_RATE -f image2 -s $RESOLUTION -start_number 1 -i $FILENAME.%04d.png -vframes $END_FRAME -vcodec libx264 -crf 25  -pix_fmt yuv420p $OUTPUT_VIDEO_FILENAME.mp4
