"""
Python functions and structures for COSMOS programs

"""
import numpy as np
from astropy.io import fits
from pyds9 import *
import os.path
import sys
from ctypes import *

class obsdata(Structure):
    _fields_ = [("instrument",c_char*10),
                ("gr_order",c_int),
                ("gr_angle",c_float),
                ("temp",c_float),
                ("alignrot",c_float),
                ("grating",c_char*20),
                ("dewoff",c_char*80),
                ("camera",c_char*6),
                ("mode",c_char*7),
                ("mask",c_char*80),
                ("nshuffle",c_int),
                ("ranod",c_float),
                ("decnod",c_float),
                ("dewar",c_char*15),
                ("distor",c_char*80)]

#
#display8
#
def dispn(imagefile,display,hd,plane,frame,nchip):

    lane="frame frameno " + frame
    display.set(lane)
    lane="scale mode zscale"
    display.set(lane)
    if(nchip==1):
        filename=imagefile + '.fits'
        hdu=fits.open(filename)
        omage=hdu[hd].data
        dims=omage.ndim
        if(dims>2):
            image=omage[plane,:,:]
        else:
            image=omage
        binfac=hdu[0].header["BINNING"]
        xb,yb=binfac.split('x')
        xbin=int(xb)
        ybin=int(yb)
    else:
        for nchip in xrange(1,9,1):
            filename=imagefile + 'c' +str(nchip) + '.fits'
            if not(os.path.isfile(filename)):
                print "\nfile", filename, "cannot be found"
                sys.exit()
            print "reading file %s" % filename
            hdu=fits.open(filename)
            if nchip==1:
                binfac=hdu[hd].header["BINNING"]
                xb,yb=binfac.split('x')
                xbin=int(xb)
                ybin=int(yb)
                if xbin != ybin:
                    print "display8 does not support unequal x and y binning\n"
                    sys.exit(1)
                if xbin==1:
                    xoff = [0,2048,4096,6144,2048,0,6144,4096]
                    image=np.zeros((8192,8192),dtype=np.float32)
                    big=4096
                    small=2048
                elif xbin==2:
                    xoff=[0,1024,2048,3072,1024,0,3072,2048]
                    image=np.zeros((4096,4096),dtype=np.float32)
                    big=2048
                    small=1024
                elif xbin==4:
                    xoff[0,512,1024,1536,512,0,1536,1024]
                    big=1024
                    small=512
                else:
                    print "display8 only supports binning of 1x1, 2x2, or 4x4\n"
                    sys.exit(1)

            omage=hdu[0].data
            dims=omage.ndim
            if(dims>2):
                chap=omage[0,:,:]
            else:
                chap=omage
            yoff=big
            chip=chap[plane:big,0:small]
            if nchip > 4:
                chip  =  np.fliplr(chip)
                yoff=0
            else:
                chip = np.flipud(chip)
            image[yoff:yoff+big, xoff[nchip-1]:xoff[nchip-1]+small] = chip[0:big,0:small]
        #put binning info in pixels 0,0 and 0,1

    image[0,0]=float(xbin)
    image[0,1]=float(ybin)
    display.set_np2arr(image)

#    binn=np.zeros((2,1),dtype=np.int)
#    binn[0,0]=xbin
#    binn[1,0]=ybin
#    display.set_np2arr(binn)

#
#die
#
def die(message):

    print message, "\n"
    sys.exit(0)
