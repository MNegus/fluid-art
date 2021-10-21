#!/bin/bash

ffmpeg -f image2 -i frames/frame_%d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p droplet.mp4

ffmpeg -i droplet.mp4 -vf palettegen palette.png
ffmpeg -i droplet.mp4 -i palette.png -lavfi paletteuse droplet.gif
