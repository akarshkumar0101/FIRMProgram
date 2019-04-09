# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 20:46:59 2018

@author: akars
"""
import numpy as math
import cv2



def ultimateEfficiency(vg, v0, linewidth):
    try:
        a = linewidth
        answer = 1
        answer *= vg
        answer /= (math.pi + 2 * math.arctan(v0 / a))
        answer /= (math.power(a, 2) + math.power(v0, 2))
        mainfunc = 0
    
        mainfunc += (math.pi * v0)
        mainfunc += (a * math.log(1 - 2 * v0 / vg + (math.power(a, 2) + math.power(v0, 2)) / math.power(vg, 2)))
        mainfunc += (2 * v0 * math.arctan((v0 - vg) / a))
    
        answer *= mainfunc
        
        return answer
    except:
        return 0

width = 2000
height = 2000
img_gray = math.zeros((width,height),dtype=math.uint8);


J_TO_eV = 6.241509343260179E18
H_IN_J_S = 6.62607004081E-34
H_IN_eV_S = H_IN_J_S * J_TO_eV


xs,ys = math.mgrid[:width,:height];

v0ev = ys*8/height
vgev = xs*8/width

v0 = v0ev/H_IN_eV_S
vg = vgev/H_IN_eV_S

efficiencies = ultimateEfficiency(vg,v0,math.power(10,9))

img_gray = math.uint8(efficiencies*255);        

img_color = cv2.applyColorMap(img_gray, cv2.COLORMAP_JET)




topspace = 200
leftspace = 200
bordercolor = 255
img_final = math.zeros((width+leftspace,height+topspace,3),dtype=math.uint8);
img_final[:,:,:] = bordercolor

img_final[leftspace:width+leftspace,topspace:height+topspace,:] = img_color

"""
lasers = open("../lasers.txt").read().split("\n")
semiconductors = open("../semiconductors.txt")

lasers = []
freqs = []
for line in open("../semiconductors.txt"):
    try:
        name,freq = line.split(" ")
        freq = float(freq)
        lasers.append(line)
        freqs.append(freq)
    except:
        pass

for laser, freq in zip(lasers,freqs):
    print(laser)
    y = int(freq*height/8);
    print(freq)
    cv2.putText(img_final,laser, (50,y+200), cv2.FONT_HERSHEY_SIMPLEX, .3, (0,0,0))    
"""

cv2.imwrite("../python_heatmapgray.png",img_gray)
cv2.imwrite("../python_heatmapcolor.png",img_color)
cv2.imwrite("../python_heatmapfinal.png",img_final)



