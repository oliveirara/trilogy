import os
import sys
import string


import astropy.io.fits as pyfits
import numpy as np


from glob import glob
from os.path import exists
from os.path import join
from PIL import Image
from PIL import ImageDraw
from scipy.optimize import golden


defaultvalues = {
    "indir": "",
    "outname": None,
    "outdir": "",
    "noiselum": 0.15,
    "noiselums": {},
    "satpercent": 0.001,
    "colorsatfac": 1,
    "thumbnail": None,
    "samplesize": 1000,
    "sampledx": 0,
    "sampledy": 0,
    "stampsize": 1000,
    "testfirst": 1,
    "show": 0,
    "showstamps": 0,
    "showwith": "open",
    "deletetests": 0,
    "deletefilters": 1,
    "scaling": None,
    "maxstampsize": 6000,
    "legend": 0,
    "invert": 0,
    "combine": "average",
    "noise": None,
    "saturate": None,
    "bscale": 1,
    "bzero": 0,
    "correctbias": 0,
    "noisesig": 1,
    "noisesig0": 2,
}

# Luminance vector
# All pretty similar; yellow galaxy glow extended a bit more in NTSC
# rw, gw, bw = 0.299,  0.587,  0.114   # NTSC (also used by PIL in "convert")
# rw, gw, bw = 0.3086, 0.6094, 0.0820  # linear
# D65: red boosted, blue muted a bit, I like it
rw, gw, bw = 0.212671, 0.715160, 0.072169


imfilt = ""  # Initialize


# Some functions
################


def rms(x):
    return np.sqrt(np.mean(x ** 2))


class meanstd_robust:
    # Generates robust statistics using a sigma clipping
    # algorithm. It is controlled by the parameters n_sigma
    # and n, the number of iterations
    # ADAPTED from Txitxo's stat_robust
    def __init__(self, x, n_sigma=3, n=5, sortedalready=False):
        self.x = x
        self.n_sigma = n_sigma
        self.n = n
        self.sortedalready = sortedalready

    def run(self):
        ihi = nx = len(self.x)
        ilo = 0

        if not self.sortedalready:
            print("sorting...")
            self.xsort = np.sort(self.x)
        else:
            self.xsort = self.x

        for i in range(self.n):
            xs = self.xsort[ilo:ihi]
            imed = int((ilo + ihi) / 2)
            aver = xs[imed]
            std1 = np.std(xs)
            std1 = rms(xs - aver)
            lo = aver - self.n_sigma * std1
            hi = aver + self.n_sigma * std1
            ilo = np.searchsorted(self.xsort, lo)
            ihi = np.searchsorted(self.xsort, hi, side="right")
            nnx = ihi - ilo
            if nnx == nx:
                break
            else:
                nx = nnx

        xrem = xs[ilo:ihi]
        self.mean = np.mean(xrem)
        self.std = rms(xrem - self.mean)


def strend(str, phr):
    return str[-len(phr) :] == phr


def decapfile(name, ext=""):
    """Remove extension from filename if present.
    If ext left blank, then any extension will be removed"""
    if ext:
        if ext[0] != ".":
            ext = "." + ext
        n = len(ext)
        if name[-n:] == ext:
            name = name[:-n]
    else:
        if strend(name, ".gz"):
            name = name[:-3]
        if strend(name, ".fz"):
            name = name[:-3]
        i = name.rfind(".")
        if i > -1:
            name = name[:i]
    return name


def str2num(str, rf=0):
    """Converts a string to a number (int or float) if possible
    Also returns format if rf=1"""
    try:
        num = int(str)
        format = "d"
    except:
        try:
            num = float(str)
            format = "f"
        except:
            if not str.strip():
                num = None
                format = ""
            else:
                words = str.split()
                if len(words) > 1:
                    num = map(str2num, tuple(words))
                    format = "l"
                else:
                    num = str
                    format = "s"
    if rf:
        return (num, format)
    else:
        return num


def clip2(m, m_min=None, m_max=None):
    if m_min == None:
        m_min = min(m)
    if m_max == None:
        m_max = max(m)
    return np.clip(m, m_min, m_max)


def striskey(str):
    """Is str an option like -C or -ker
    (It's not if it's -2 or -.9)"""
    iskey = 0
    if str:
        if str[0] == "-":
            iskey = 1
            if len(str) > 1:
                iskey = str[1] not in [
                    "0",
                    "1",
                    "2",
                    "3",
                    "4",
                    "5",
                    "6",
                    "7",
                    "8",
                    "9",
                    ".",
                ]
    return iskey


def params_cl(converttonumbers=True):
    """Returns parameters from command line ('cl') as dictionary:
    Keys are options beginning with '-'
    Values are whatever follows keys: either nothing (''), a value, or a list of values
    all values are converted to int / float when appropriate"""
    list = sys.argv[:]
    i = 0
    dict = {}
    key = ""
    list.append("")  # extra element so we come back and assign the last value
    while i < len(list):
        if striskey(list[i]) or not list[i]:  # (or last value)
            if key:  # assign values to old key
                if value:
                    if len(value) == 1:  # list of 1 element
                        value = value[0]  # just element
                dict[key] = value
            if list[i]:
                key = list[i][1:]  # remove leading '-'
                value = None
                dict[key] = value  # in case there is no value!
        else:  # value (or haven't gotten to keys)
            if key:  # (have gotten to keys)
                if value:
                    if converttonumbers:
                        value.append(str2num(list[i]))
                    else:
                        value = value + " " + list[i]
                else:
                    if converttonumbers:
                        value = [str2num(list[i])]
                    else:
                        value = list[i]
        i += 1

    return dict


# TRILOGY-specific tools
########################


def determinescaling(data, unsatpercent, noisesig=1, correctbias=True, noisesig0=2):
    """Determines data values (x0,x1,x2) which will be scaled to (0,noiselum,1)"""
    # Robust mean & standard deviation
    datasorted = np.sort(data.flat)
    datasorted[np.isnan(datasorted)] = 0  # set all nan values to zero
    if datasorted[0] == datasorted[-1]:
        levels = 0, 1, 100  # whatever
    else:
        s = meanstd_robust(datasorted, sortedalready=True)
        s.run()
        m = s.mean
        r = s.std
        if correctbias:
            x0 = m - noisesig0 * r
        else:
            x0 = 0
        x1 = m + noisesig * r
        x2 = setlevels(datasorted, np.array([unsatpercent]), sortedalready=True)[0]
        levels = x0, x1, x2
    return levels


def setlevels(data, pp, stripneg=False, sortedalready=False):
    if sortedalready:
        vs = data
    else:
        print("sorting...")
        vs = np.sort(data.flat)
    if stripneg:  # Get rid of negative values altogether!
        i = np.searchsorted(vs, 0)
        vs = vs[i + 1 :]
    else:  # Clip negative values to zero
        vs = clip2(vs, 0, None)
    ii = np.array(pp) * len(vs)
    ii = ii.astype(int)
    ii = np.clip(ii, 0, len(vs) - 1)
    levels = vs.take(ii)
    return levels


def da(k):
    a1 = k * (x1 - x0) + 1
    a2 = k * (x2 - x0) + 1
    a1n = a1 ** n
    a1n = abs(a1n)  # Don't want the solutions where a1 & a2 are both negative!
    da1 = a1n - a2
    k = abs(k)
    if k == 0:
        return da(1e-10)
    else:
        da1 = da1 / k  # To avoid solution k = 0!
    return abs(da1)


def imscale2(data, levels, y1):
    # x0, x1, x2  yield  0, y1, 1,  respectivelly
    global n, x0, x1, x2  # So that golden can use them
    # Normalize?  No.  Unless the data is all ~1e-40 or something...
    x0, x1, x2 = levels
    if y1 == 0.5:
        k = (x2 - 2 * x1 + x0) / float(x1 - x0) ** 2
    else:
        n = 1 / y1
        k = abs(golden(da))
    r1 = np.log10(k * (x2 - x0) + 1)
    v = np.ravel(data)
    v = clip2(v, 0, None)
    d = k * (v - x0) + 1
    d = clip2(d, 1e-30, None)
    z = np.log10(d) / r1
    z = np.clip(z, 0, 1)
    z.shape = data.shape
    z = z * 255
    z = z.astype(np.uint8)
    return z


def satK2m(K):
    m00 = rw * (1 - K) + K
    m01 = gw * (1 - K)
    m02 = bw * (1 - K)

    m10 = rw * (1 - K)
    m11 = gw * (1 - K) + K
    m12 = bw * (1 - K)

    m20 = rw * (1 - K)
    m21 = gw * (1 - K)
    m22 = bw * (1 - K) + K

    m = np.array([[m00, m01, m02], [m10, m11, m12], [m20, m21, m22]])
    return m


def adjsat(RGB, K):
    """Adjust the color saturation of an image.  K > 1 boosts it."""
    m = satK2m(K)
    three, nx, ny = RGB.shape
    RGB.shape = three, nx * ny
    RGB = np.dot(m, RGB)
    RGB.shape = three, nx, ny
    return RGB


def RGB2im(RGB):
    """r, g, b = data  (3, ny, nx)
    Converts to an Image"""
    data = RGB
    data = np.transpose(data, (1, 2, 0))  # (3, ny, nx) -> (ny, nx, 3)
    data = np.clip(data, 0, 255)
    data = data.astype(np.uint8)
    three = data.shape[-1]  # 3 if RGB, 1 if L
    if three == 3:
        im = Image.fromarray(data)
    elif three == 1:
        im = Image.fromarray(data[:, :, 0], "L")
    else:
        print(
            "Data shape not understood: expect last number to be 3 for RGB, 1 for L",
            data.shape,
        )
        raise Exception  # Raise generic exception and exit
    im = im.transpose(Image.FLIP_TOP_BOTTOM)
    return im


def RGBscale2im(RGB, levdict, noiselums, colorsatfac, mode="RGB", invlum=0):
    three, nx, ny = RGB.shape  # if 'L', then three = 1 !
    if nx * ny > 2000 * 2000:
        print("Warning: You should probably feed smaller stamps into RGBscale2im.")
        print("This may take a while...")
    scaled = np.zeros(RGB.shape, float)
    for i in range(three):
        channel = mode[i]  # 'RGB' or 'L'
        levels = levdict[channel]
        noiselum = noiselums[channel]
        scaled[i] = imscale2(RGB[i], levels, noiselum)
    if (colorsatfac != 1) and (mode == "RGB"):
        scaled = adjsat(scaled, colorsatfac)
    if invlum:
        scaled = 255 - scaled
    im = RGB2im(scaled)
    return im


def loadfile(filename, dir="", silent=0, keepnewlines=0):
    infile = join(dir, filename)
    if not silent:
        print("Loading ", infile, "...\n")
    fin = open(infile, "r")
    sin = fin.readlines()
    fin.close()
    if not keepnewlines:
        for i in range(len(sin)):
            sin[i] = sin[i][:-1]
    return sin


def loaddict(filename, dir="", silent=0):
    lines = loadfile(filename, dir, silent)
    dict = {}
    for line in lines:
        if line[0] != "#":
            words = line.split()
            key = str2num(words[0])
            val = ""  # if nothing there
            valstr = " ".join(words[1:])
            valtuple = False
            valarray = True
            if valstr[0] in "[(" and valstr[-1] in "])":  # LIST / TUPLE!
                valtuple = valstr[0] == "("
                valstr = valstr[1:-1].replace(",", "")
                words[1:] = valstr.split()
            if len(words) == 2:
                val = str2num(words[1])
            elif len(words) > 2:
                val = []
                for word in words[1:]:
                    val.append(str2num(word))
                if valtuple:
                    val = tuple(val)
                if valarray:
                    val = np.array(val)
            dict[key] = val
    return dict


# Apply offsets
###############

offsets = {}


def offsetarray(data, offset):
    new = np.zeros(data.shape)
    dx, dy = offset
    if dy >= 0:
        if dx >= 0:
            new[dy:, dx:] = data[:-dy, :-dx]
        else:
            new[dy:, :-dx] = data[:-dy, dx:]
    else:
        if dx >= 0:
            new[:-dy, dx:] = data[dy:, :-dx]
        else:
            new[:-dy, :-dx] = data[dy:, dx:]
    return new


# for channel in offsets.keys():
#     dataRGB[channel] = offsetarray(dataRGB[channel], offsets[channel])


def extractfilter(header):
    """Extracts filter from a FITS header"""
    filt = header.get("FILTER", "")
    if filt == "":
        filt1 = header.get("FILTER1", "")
        if filt1 != "":
            if type(filt1) == str:
                if filt1[:5] == "CLEAR":
                    filt2 = header.get("FILTER2", "")
                    filt = filt2
                else:
                    filt = filt1
    return filt


def loadfitsimagedata(image, indir="", silent=1, bscale=1, bzero=0):
    global imfilt, imfilts
    if image[-1] == "]":
        iext = int(image[-2])
        image = image[:-3]  # Remove [0]
    else:
        iext = 0

    image0 = decapfile(image)

    image = join(indir, image)
    hdulist = pyfits.open(image, memmap=1)
    hdu = hdulist[iext]
    data = hdu.data
    if image0 in imfilts.keys():
        imfilt = imfilts[image0]
    else:
        imfilt = extractfilter(hdu.header)
        imfilts[image0] = imfilt
    if not silent:
        print(image + "[%d]" % iext, data.shape, imfilt)
    return data


imfilts = {}


def processimagename(image):
    global imfilts
    if image[-1] == ")":
        i = image.find("(")
        filt = image[i + 1 : -1]
        image = image[:i]
        imfilts[image] = filt
    if image[-1] == "]":
        ext = image[-3:]
        image = image[:-3]
    else:
        ext = ""
    if (
        not strend(image, "fits")
        and not strend(image, "fits.gz")
        and not strend(image, "fits.fz")
    ):
        image += ".fits"
    image = image + ext
    return image


def datascale(data, bscale, bzero):
    if (bscale != 1) or (bzero != 0):
        return bscale * data + bzero
    else:
        return data


class Trilogy:
    def __init__(self, infile=None, images=None, imagesorder="BGR", **inparams):
        self.nx = None  # image size
        self.ny = None  # image size
        self.xlo = 0
        self.ylo = 0
        self.xhi = None
        self.yhi = None
        self.imagesRGB = {"R": [], "G": [], "B": [], "L": []}  # File names
        self.inkeys = []
        self.mode = "L"  # reset below if color
        self.weightext = None  # No weighting unless weight images are declared
        # Can use either:
        # weightext drz wht
        # weightext drz -> wht

        print("From input file", infile, ":")
        self.infile = infile
        if infile:
            self.loadinputs()

        self.images = images
        if images:
            self.setimages()

        self.inparams = inparams
        if inparams:
            self.setinparams()

        self.setdefaults()

    def setinparams(self):
        print("From input parameters:")
        bekeys = "invert ".split()
        inkeys = self.inparams.keys()
        for key in inkeys:
            if key in bekeys:
                continue
            val = self.inparams[key]
            cmd = "self.%s = val" % key
            print(key, "=", val)
            exec(cmd)
            self.inkeys.append(key)
        for key in bekeys:
            val = key in inkeys
            cmd = "self.%s = val" % key
            print(key, "=", val)
            exec(cmd)
            self.inkeys.append(key)  # added later

    def setdefaults(self):
        print("Default:")
        for key in defaultvalues.keys():
            if key not in self.inkeys:
                val = defaultvalues[key]
                cmd = "self.%s = val" % key
                exec(cmd)
                print(key, "=", val, "(default)")

    def setimages(self, images=None):
        images = images or self.images
        if images != None:
            if type(images) == str:  # Single image
                images = processimagename(images)
                self.imagesRGB["L"] = [images]
                self.mode = "L"
            elif type(images[0]) == str:  # List of images
                images = map(processimagename, images)
                self.imagesRGB["L"] = images
                self.mode = "L"
            else:  # List of 3 lists of images, one for each channel
                self.mode = "RGB"
                for i in range(3):
                    channel = imagesorder[i]
                    channelimages = map(processimagename, images[i])
                    self.imagesRGB[channel] = channelimages

    def setnoiselums(self):
        for channel in self.mode:
            self.noiselums[channel] = self.noiselum

    def loadinputs(self):
        """Load R,G,B filenames and options"""
        f = open(self.infile)
        prevline = ""
        channel = "L"  # if no channel declared, then it's grayscale!
        self.noiselums = {}
        for line in f:
            if line[0] == "#":
                continue

            word = line.strip()
            if len(word):
                words = word.split()
                if len(words) == 1:  # Channel or image name
                    if (word in "RGB") and (prevline == ""):
                        channel = word
                        self.mode = "RGB"
                    else:
                        image = word
                        image = processimagename(image)
                        self.imagesRGB[channel].append(image)
                        print(channel, image)
                else:  # parameter and value(s)
                    key = words[0]
                    val = str2num("".join(words[1:]))
                    if key == "weightimages":
                        print(words)
                        if len(words[1:]) == 2:  # drz wht
                            self.imext, self.weightext = words
                        else:  # drz -> wht
                            self.imext = words[1]
                            self.weightext = words[3]
                    elif key == "noiselums":
                        if "," in val:
                            val = val.split(",")
                            val = map(float, val)
                        for i, channel in enumerate(self.mode[::-1]):
                            self.noiselums[channel] = val[i]
                    else:
                        cmd = "self.%s = val" % key
                        # print(cmd)
                        exec(cmd)
                    print(key, "=", val)
                    self.inkeys.append(key)
            prevline = word

        if self.noiselums == {}:
            if "noiselum" in self.inkeys:
                for channel in self.mode:
                    self.noiselums[channel] = self.noiselum
                    self.inkeys.append("noiselums")
        f.close()

    def setoutfile(self, outname=None):
        self.outname = outname or self.outname or self.infile
        if self.outname == None:
            self.outname = self.images
            if self.outname:
                self.outname = os.path.basename(self.outname)
                self.outname = decapfile(self.outname)
        if (len(self.outname) > 4) and (self.outname[-4] == "."):
            # Has extension
            self.outfile = self.outname  # Use whatever extension they picked
            self.outname = self.outname[:-4]  # Remove extension
        else:  # Just root
            self.outfile = self.outname + ".png"

    def outfilterfile(self):
        outfile = self.outname + "_filters.txt"
        outfile = join(self.outdir, outfile)
        return outfile

    def loadimagesize(self):
        global filters
        print(
            "Loading image data.",
        )
        print("If multiple filters per channel, adding data.")
        filters = {"B": [], "G": [], "R": [], "L": []}
        fout = open(self.outfilterfile(), "w")

        for channel in self.mode[::-1]:
            if channel in "RGB":
                ichannel = "RGB".index(channel)
            else:  # L
                ichannel = 0
            outline = channel + " = "
            for iimage, image in enumerate(self.imagesRGB[channel]):
                print(channel)
                sgn = 1
                if image[0] == "-":
                    sgn = -1
                    image = image[1:]
                if iimage:
                    sgnsym = ".+-"[sgn]
                    outline += " %s " % sgnsym
                elif sgn == -1:
                    outline += "- "
                data = loadfitsimagedata(image, self.indir, silent=0)
                filt = imfilt
                if filt == "":
                    filt = imfilts.get(image, "")
                filters[channel].append(filt)
                outline += filt
                ny, nx = data.shape
                if self.ny == None:
                    self.ny = ny
                    self.nx = nx
                    self.yc = ny / 2
                    self.xc = nx / 2
                else:
                    if (self.ny != ny) or (self.nx != nx):
                        print(
                            "Input FAIL.  Your images are not all the same size as (%d,%d)."
                            % (self.ny, self.nx)
                        )
                        for channel in self.mode[::-1]:  # 'BGR'
                            for image in self.imagesRGB[channel]:
                                data = loadfitsimagedata(image, self.indir, silent=0)
                        raise  # Raise Exception (error) and quit

            fout.write(outline + "\n")
            print(outline)
            print()
        fout.close()

        fout = open(self.outfilterfile(), "w")
        for channel in "BGR":
            filtstr = " + ".join(filters[channel])
            fout.write("%s = %s\n" % (channel, filtstr))
        fout.close()

        # Allow the user to just create part of the image
        # xlo  1000
        # xhi  5000
        # xhi -1000
        # if hi is negative, then trims from that edge (as 1000:-1000)
        if self.xhi == None:
            self.xhi = self.nx
        elif self.xhi < 0:
            self.xhi = self.nx + self.xhi

        if self.yhi == None:
            self.yhi = self.ny
        elif self.yhi < 0:
            self.yhi = self.ny + self.yhi

    def loadstamps(self, limits, silent=1):
        ylo, yhi, xlo, xhi = limits

        ylo = np.clip(ylo, 0, self.ny)
        yhi = np.clip(yhi, 0, self.ny)
        xlo = np.clip(xlo, 0, self.nx)
        xhi = np.clip(xhi, 0, self.nx)

        ny = yhi - ylo
        nx = xhi - xlo

        three = len(self.mode)
        stampRGB = np.zeros((three, ny, nx), float)
        weighting = self.weightext != None
        if weighting:
            weightstampRGB = np.zeros((three, ny, nx), float)

        for ichannel, channel in enumerate(self.mode):
            for image in self.imagesRGB[channel]:
                if not silent:
                    print(
                        channel,
                    )
                sgn = 1
                if image[0] == "-":
                    sgn = -1
                    image = image[1:]

                data = loadfitsimagedata(image, self.indir, silent=silent)
                stamp = data[ylo:yhi, xlo:xhi]
                stamp = datascale(stamp, self.bscale, self.bzero)

                # weight image?
                if weighting:
                    weightimage = image.replace(self.imext, self.weightext)
                    weightfile = join(self.indir, weightimage)
                    if exists(weightfile):
                        weight = loadfitsimagedata(
                            weightimage, self.indir, silent=silent
                        )
                        weightstamp = weight[ylo:yhi, xlo:xhi]
                        # FLAG IMAGE!!  EITHER 1 or 0
                        weightstamp = np.greater(weightstamp, 0)
                        weightstampRGB[ichannel] = (
                            weightstampRGB[ichannel] + sgn * weightstamp
                        )
                        stamp = stamp * weightstamp
                    else:
                        print(weightfile, "DOES NOT EXIST")

                stampRGB[ichannel] = stampRGB[ichannel] + sgn * stamp

        if weighting:
            for ichannel, channel in enumerate(self.mode):
                stampRGB[ichannel] = np.where(
                    weightstampRGB[ichannel],
                    stampRGB[ichannel] / weightstampRGB[ichannel],
                    0,
                )
        elif self.combine == "average":
            for ichannel, channel in enumerate(self.mode):
                stampRGB[ichannel] = stampRGB[ichannel] / len(self.imagesRGB[channel])

        return stampRGB

    def determinescalings(self):
        """Determine data scalings
        will sample a (samplesize x samplesize) region of the (centered) core
        make color image of the core as a test if desired"""
        self.testimages = []
        redo = True
        while redo:  # Until user is happy with test image of core
            dx = dy = self.samplesize
            print()
            unsatpercent = 1 - 0.01 * self.satpercent
            self.levdict = {}

            if dx * dy == 0:
                print(
                    "By setting samplesize = 0, you have asked to sample the entire image to determine the scalings."
                )
                print(
                    "(Note this will be clipped to a maximum of %dx%d.)"
                    % (self.maxstampsize, self.maxstampsize)
                )
                dx = dy = self.maxstampsize  # Maximum size possible

            ylo = np.clip(self.yc - dy / 2 + self.sampledy, 0, self.ny)
            yhi = np.clip(self.yc + dy / 2 + self.sampledy, 0, self.ny)
            xlo = np.clip(self.xc - dx / 2 + self.sampledx, 0, self.nx)
            xhi = np.clip(self.xc + dx / 2 + self.sampledx, 0, self.nx)
            dy = yhi - ylo
            dx = xhi - xlo
            print(
                "Determining image scaling based on %dx%d core sample" % (dx, dy),
            )
            if self.sampledx or self.sampledy:
                print(
                    "offset by (%d,%d)" % (self.sampledx, self.sampledy),
                )

            print("...")

            limits = int(ylo), int(yhi), int(xlo), int(xhi)
            stampRGB = self.loadstamps(limits)
            for ichannel, channel in enumerate(self.mode):
                self.levdict[channel] = determinescaling(
                    stampRGB[ichannel],
                    unsatpercent,
                    noisesig=self.noisesig,
                    correctbias=self.correctbias,
                    noisesig0=self.noisesig0,
                )
                print(
                    channel,
                )
                print(" %f  %f  %f" % self.levdict[channel])

            redo = False
            if self.testfirst:
                im = RGBscale2im(
                    stampRGB,
                    self.levdict,
                    self.noiselums,
                    self.colorsatfac,
                    self.mode,
                    self.invert,
                )

                outfile = "%s_test_%g_%g_%g.png" % (
                    self.outname,
                    self.satpercent,
                    self.noiselum,
                    self.colorsatfac,
                )
                outfile = join(self.outdir, outfile)
                self.testimages.append(outfile)

                print("Creating test image", outfile)
                im.save(outfile)

                # NOTE I use "open" instead of im.show()
                # because the latter converts the image to a jpg for display
                # which degrades it slightly
                if self.show:
                    self.showimage(outfile, Image)

                print("Like what you see?")
                print("If so, press <Enter> a few times")

                print("Otherwise, enter new values:")

                line = "  noise yields brightness: %g? " % self.noiselum
                inp = input(line)
                if inp.strip() != "":
                    self.noiselum = float(inp)
                    for channel in self.mode:
                        self.noiselums[channel] = self.noiselum
                    redo = True

                line = "  %% of pixels that saturate: %g? " % self.satpercent
                inp = input(line)
                if inp.strip() != "":
                    self.satpercent = float(inp)
                    redo = True

                if self.mode == "RGB":
                    line = "  color saturation factor: %g? " % self.colorsatfac
                    inp = input(line)
                    if inp.strip() != "":
                        self.colorsatfac = float(inp)
                        redo = True

                line = "  Sample size: %d? " % self.samplesize
                inp = input(line)
                if inp.strip() != "":
                    self.samplesize = int(inp)
                    redo = True

                line = "  Sample offset x: %d? " % self.sampledx
                inp = input(line)
                if inp.strip() != "":
                    self.sampledx = int(inp)
                    redo = True

                line = "  Sample offset y: %d? " % self.sampledy
                inp = input(line)
                if inp.strip() != "":
                    self.sampledy = int(inp)
                    redo = True

    def determinescalings2(self):
        """Determine data scalings
        will sample a (samplesize x samplesize) region of the (centered) core
        make color image of the core as a test if desired"""
        self.testimages = []
        redo = True
        while redo:  # Until user is happy with test image of core
            dx = dy = self.samplesize
            print()
            unsatpercent = 1 - 0.01 * self.satpercent
            self.levdict = {}

            if dx * dy == 0:
                print(
                    "By setting samplesize = 0, you have asked to sample the entire image to determine the scalings."
                )
                print(
                    "(Note this will be clipped to a maximum of %dx%d.)"
                    % (self.maxstampsize, self.maxstampsize)
                )
                dx = dy = self.maxstampsize  # Maximum size possible

            ylo = np.clip(self.yc - dy / 2 + self.sampledy, 0, self.ny)
            yhi = np.clip(self.yc + dy / 2 + self.sampledy, 0, self.ny)
            xlo = np.clip(self.xc - dx / 2 + self.sampledx, 0, self.nx)
            xhi = np.clip(self.xc + dx / 2 + self.sampledx, 0, self.nx)

            dy = yhi - ylo
            dx = xhi - xlo
            print(
                "Determining image scaling based on %dx%d core sample" % (dx, dy),
            )

            if self.sampledx or self.sampledy:
                print(
                    "offset by (%d,%d)" % (self.sampledx, self.sampledy),
                )

            print("...")

            limits = ylo, yhi, xlo, xhi
            stampRGB = self.loadstamps(limits)
            for channel in self.mode:
                self.levdict[channel] = 0, self.noise, self.saturate

            redo = False
            if self.testfirst:
                im = RGBscale2im(
                    stampRGB,
                    self.levdict,
                    self.noiselums,
                    self.colorsatfac,
                    self.mode,
                    self.invert,
                )

                outfile = "%s_test_%g_%g_%g.png" % (
                    self.outname,
                    self.noiselum,
                    self.noise,
                    self.saturate,
                )
                outfile = join(self.outdir, outfile)
                self.testimages.append(outfile)

                print("Creating test image", outfile)
                im.save(outfile)

                # Note: I use "open" instead of im.show()
                # because the latter converts the image to a jpg for display
                # which degrades it slightly
                if self.show:
                    self.showimage(outfile, Image)

                print("Like what you see?")
                print("If so, press <Enter> a few times")

                print("Otherwise, enter new values:")

                line = "  noise yields brightness: %g? " % self.noiselum
                inp = input(line)
                if inp.strip() != "":
                    self.noiselum = float(inp)
                    for channel in self.mode:
                        self.noiselums[channel] = self.noiselum
                    redo = True

                line = "  noise image input data value: %g? " % self.noise
                inp = input(line)
                if inp.strip() != "":
                    self.noise = float(inp)
                    redo = True

                line = "  saturation level image input data value: %g? " % self.saturate
                inp = input(line)
                if inp.strip() != "":
                    self.saturate = float(inp)
                    redo = True

                if self.mode == "RGB":
                    line = "  color saturation factor: %g? " % self.colorsatfac
                    inp = input(line)
                    if inp.strip() != "":
                        self.colorsatfac = float(inp)
                        redo = True

                line = "  Sample size: %d? " % self.samplesize
                inp = input(line)
                if inp.strip() != "":
                    self.samplesize = int(inp)
                    redo = True

                line = "  Sample offset x: %d? " % self.sampledx
                inp = input(line)
                if inp.strip() != "":
                    self.sampledx = int(inp)
                    redo = True

                line = "  Sample offset y: %d? " % self.sampledy
                inp = input(line)
                if inp.strip() != "":
                    self.sampledy = int(inp)
                    redo = True

    def makecolorimage(self):
        """Make color image (in sections)"""
        if (self.stampsize == self.samplesize == 0) and self.testfirst:
            # Already did the full image!
            print("Full size image already made.")
            imfile = self.testimages[-1]
            outfile = join(self.outdir, self.outfile)
            if self.deletetests:
                print("Renaming to", outfile)
                os.rename(imfile, outfile)
            else:
                print("Copying to", outfile)
                os.copy(imfile, outfile)
            imfull = Image.open(outfile)
            return imfull

        # Clean up: Delete test images
        if self.deletetests:
            for testimage in self.testimages:
                if exists(testimage):
                    os.remove(testimage)

        dx = dy = self.stampsize
        if dx * dy == 0:
            dx = dy = self.maxstampsize

        imfull = Image.new(self.mode, (self.nx, self.ny))

        print()
        if self.mode == "RGB":
            print("Making full color image, one stamp (section) at a time...")
        elif self.mode == "L":
            print("Making full grayscale image, one stamp (section) at a time...")
        for yo in range(self.ylo, self.yhi, dy):
            dy1 = min([dy, self.yhi - yo])
            for xo in range(self.xlo, self.xhi, dx):
                dx1 = min([dx, self.xhi - xo])
                print("%5d, %5d  /  (%d x %d)" % (xo, yo, self.nx, self.ny))
                limits = yo, yo + dy, xo, xo + dx
                stamps = self.loadstamps(limits)
                im = RGBscale2im(
                    stamps,
                    self.levdict,
                    self.noiselums,
                    self.colorsatfac,
                    self.mode,
                    self.invert,
                )
                if self.show and self.showstamps:
                    im.show()

                imfull.paste(im, (xo, self.ny - yo - dy1, xo + dx1, self.ny - yo))

        outfile = join(self.outdir, self.outfile)
        if self.legend:
            self.addlegend(im=imfull)
        else:
            print("Saving", outfile, "...")
            imfull.save(outfile)

        if self.show:
            self.showimage(outfile, Image)

        return imfull

    def makethumbnail1(self, outroot, width, fmt="jpg"):
        nx = width
        ny = 1000 * (nx / float(nx))
        im = open(outroot + ".png")
        im2 = im.resize((nx, ny))
        im2.save(self.outname + "_%d.%s" % (width, fmt))
        return im2

    def makethumbnail(self):
        if self.thumbnail not in [None, "None"]:
            outname = self.thumbnail
            if outname[-4] == ".":
                outname = outname[:-4]
                fmt = outname[-3:]
            width = int(outname)
            self.makethumbnail1(self.outname, width, fmt)

    def showsample(self, outfile):
        dx = dy = self.samplesize
        if dx * dy == 0:
            print(
                "By setting samplesize = 0, you have asked to sample the entire image to determine the scalings."
            )
            print(
                "(Note this will be clipped to a maximum of %dx%d.)"
                % (self.maxstampsize, self.maxstampsize)
            )
            dx = dy = self.maxstampsize  # Maximum size possible

        ylo = np.clip(self.yc - dy / 2 + self.sampledy, 0, self.ny)
        yhi = np.clip(self.yc + dy / 2 + self.sampledy, 0, self.ny)
        xlo = np.clip(self.xc - dx / 2 + self.sampledx, 0, self.nx)
        xhi = np.clip(self.xc + dx / 2 + self.sampledx, 0, self.nx)

        dy = yhi - ylo
        dx = xhi - xlo
        print(
            "Showing %dx%d core sample" % (dx, dy),
        )
        if self.sampledx or self.sampledy:
            print(
                "offset by (%d,%d)" % (self.sampledx, self.sampledy),
            )

        print("...")

        limits = ylo, yhi, xlo, xhi
        stampRGB = self.loadstamps(limits)
        im = RGBscale2im(
            stampRGB,
            self.levdict,
            self.noiselums,
            self.colorsatfac,
            self.mode,
            self.invert,
        )

        outfile = join(self.outdir, outfile)

        print("Creating test image", outfile)
        im.save(outfile)

        # Note: I use "open" instead of im.show()
        # because the latter converts the image to a jpg for display
        # which degrades it slightly
        if self.show:
            self.showimage(outfile, Image)

        print("Like what you see?")
        inp = input()

    def addlegend(self, outfile=None, im=None):
        if im == None:
            outfile1 = join(self.outdir, self.outfile)
            print("Adding legend to", outfile1, "...")
            im = Image.open(outfile1)
        else:
            print("Adding legend...")

        nx, ny = im.size
        draw = ImageDraw.Draw(im)

        x = 20
        y0 = 20
        dy = 15

        txt = loadfile(self.outfilterfile(), silent=1)

        if self.mode == "L":
            white = 255
            line = txt[0][4:]  # get rid of leading "L = "
            draw.text((x, y0), line, fill=white)
        else:
            blue = tuple(255 * np.array([0, 0.5, 1]))
            green = tuple(255 * np.array([0, 1, 0]))
            red = tuple(255 * np.array([1, 0, 0]))

            colors = blue, green, red
            colors = np.array(colors).astype(int)

            for i, line in enumerate(txt):
                y = y0 + dy * i
                ichannel = "BGR".index(line[0])
                color = tuple(colors[ichannel])
                draw.text((x, y), line, fill=color)

        if outfile == None:
            outfile = join(self.outdir, self.outfile)

        print("Saving", outfile, "...")
        im.save(outfile)

    def showimage(self, outfile, Image):
        cmd = self.showwith
        if (not cmd) or (cmd.upper() == "PIL"):
            Image.open(outfile).show()
        else:
            try:
                os.system(cmd + " " + outfile)
            except:
                # In case "open" doesn't work on their system (not a Mac)
                # Although it may not work but not raise an error either!
                # Should do better error handling here
                Image.open(outfile).show()

    def run(self):
        self.setimages()  # not needed from command line
        self.setoutfile()  # adds .png if necessary to outname
        self.loadimagesize()

        if "justaddlegend" in self.inkeys:
            self.addlegend()
            quit()

        if self.noiselums == {}:
            self.setnoiselums()
        if self.scaling == None:
            if self.noise and self.saturate:
                self.determinescalings2()
            else:
                self.determinescalings()
        else:
            print("Loading scaling saved in", self.scaling)
            self.levdict = loaddict(self.scaling)
            scaleroot = os.path.basename(self.scaling)[:-4]
            if self.testfirst:
                self.showsample(self.outname + "_" + scaleroot + ".png")
        print("Scalings:")
        for channel in self.mode:
            print(channel, self.levdict[channel])
        self.makecolorimage()
        self.makethumbnail()

        # Clean up: Delete filter files
        if self.deletefilters:
           if exists(self.outfilterfile()):
               os.remove(self.outfilterfile())



if __name__ == "__main__":
    if len(sys.argv) == 1:
        infile = "trilogy.in"
        images = None
        Trilogy(infile, images=images, **params_cl()).run()
    else:  # > 1
        input1 = sys.argv[1]
        if ("*" in input1) or ("?" in input1):
            indir = params_cl().get("indir", "")
            input1 = join(indir, input1)
            images = glob(input1)
            for image in images:
                Trilogy(images=image, **params_cl()).run()
        else:
            images = None
            if (
                strend(input1, ".fits")
                or strend(input1, ".fits.gz")
                or strend(input1, ".fits.fz")
            ):
                images = input1
                infile = None
            else:
                infile = input1

            Trilogy(infile, images=images, **params_cl()).run()