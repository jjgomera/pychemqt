import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

# Manually get point for a chart as image
# Load first the image as plot, get the x,y coord from status bar and next run
# to get the translated values

###############################################################################
# Input parameters
###############################################################################
fname = '/home/jjgomera/Imágenes/Captura de pantalla de 2023-05-22 16-26-23.png'

# x axis coords
x0 = 41
x1 = 609
xt0 = 0
xt1 = 7e3

# y axis coords
y0 = 32
y1 =448
yt0 = 100+273
yt1 = -100+273

# Point to transform
x= np.r_[88, 106, 113, 113, 125, 118, 146, 156, 162, 186]
y= np.r_[291, 264, 243, 227, 205, 187, 126, 88, 61, 35]
###############################################################################


img = np.asarray(Image.open(fname))
imgplot = plt.imshow(img)
plt.show()

bx = (xt1-xt0)/(x1-x0)
ax = xt0-bx*x0
def xt(x):
    return ax + bx*x

by = (yt1-yt0)/(y1-y0)
ay = yt0-by*y0
def yt(y):
    return ay + by*y

print(xt(x))
print(yt(y))
