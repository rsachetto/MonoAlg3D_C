from  paraview import simple
import os
from sys import argv, exit
import glob
import re

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-s", "--step", type=int, help="Animate a step fraction of the files")
parser.add_argument("-v", "--video", help="save animation video",  action="store_true")
parser.add_argument("-i", "--images", help="save animation images",  action="store_true")
parser.add_argument("-a", "--azimuth", type=float, help="a degress horizontal rotation of the scene.", default=0.0)
parser.add_argument("-e", "--elevation", type=float, help="y degress vertical rotation of the scene.", default=0.0)
parser.add_argument("source_dir", type=str, help="Directory containing the vtu files")

args = parser.parse_args()

video_file = None

if args.video or args.images:
    video_file = 'video/video.avi'
    os.system('mkdir video')

    if args.video:
        print "Video will be saved in ./video dir"
    if args.images:
        print "Images will be saved in ./images dir"
        os.system('mkdir images')


files = glob.glob(args.source_dir+"/*.vtk")
files = natural_sort(files)

jump = args.step

count = 0

files = files[::jump] #some_list[start:stop:step]

reader = simple.OpenDataFile(files)
simple.Show(reader)
dp = simple.GetDisplayProperties(reader)
dp.Representation = 'Surface'

simple.GetActiveView().GetRenderWindow().SetSize(800, 800)
dp.LookupTable = simple.MakeBlueToRedLT(-86.2, 40.0)
dp.ColorArrayName = 'Scalars_'

camera = simple.GetActiveCamera()
camera.Elevation(args.elevation)
camera.Azimuth(args.azimuth)

simple.Render()

#simple.AnimateReader(reader, filename=video_file)
simple.SaveAnimation(video_file)

if args.images:
    os.system('ffmpeg -i ' + video_file + ' ./images/output_%04d.png')
    if not args.video:
        os.system('rm -fr ./video')

