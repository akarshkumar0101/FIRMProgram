# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 20:46:59 2018

@author: akars
"""
import numpy as np
import cv2

img_gray = cv2.imread(input(),cv2.IMREAD_GRAYSCALE)
img_gray[0][0]=255
img_color = cv2.applyColorMap(img_gray, cv2.COLORMAP_JET)

cv2.imshow("img",img_gray)
cv2.waitKey(0)
cv2.destroyAllWindows()

width,height= img_color.shape[:2]


topspace = 50
leftspace = 50
bordercolor = 255
img_final = np.zeros((width+2*leftspace,height+2*topspace,3),dtype=np.uint8);
img_final[:,:,:] = bordercolor

img_final[leftspace:width+leftspace,topspace:height+topspace,:] = img_color


#cv2.imwrite("../testcolor1.png",img_color)
cv2.imwrite(input(),img_final)



