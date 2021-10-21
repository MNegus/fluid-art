# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from PIL import Image

# First image
sourceImg = Image.open("edited_1-20.png")

count = 1
w = 1128
x0 = 719
y0 = 246

rows = 4
cols = 5

for i in range(rows):
    for e in range(cols):
        x = x0 + e * w
        y = y0 + i * w
        
        crop = sourceImg.crop((x, y, x + w, y + w))
        crop.save("frames/frame_"+str(count)+".png")
        
        count += 1
        
# Second image
sourceImg = Image.open("edited_21-40.png")

w = 1127
x0 = 705
y0 = 258

rows = 4
cols = 5

for i in range(rows):
    for e in range(cols):
        x = x0 + e * w
        y = y0 + i * w
        
        crop = sourceImg.crop((x, y, x + w, y + w))
        crop.save("frames/frame_"+str(count)+".png")
        
        count += 1