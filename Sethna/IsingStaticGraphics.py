import os, random
import Image, ImageDraw
import tempfile

# Graphics for the Ising model

# -----------------------------------------------------------------------

# image display

def Display(image_file='tmpf.jpg'):
    """Display(image_file) attempts to display the specified image_file on
    the screen, using the Preview application on Mac OS X, the ImageMagick
    display utility on other posix platforms (e.g., Linux)
    and the Microsoft mspaint utility on Windows platforms."""
    os_name = os.name
    if os_name == 'nt': # Windows
        try:
            os.system('mspaint %s &' % image_file)
        except:
            raise OSError, "Cannot display %s with Windows mspaint" % \
                  image_file
    else:
      os_uname = os.uname()
      if os_uname[0] == 'Darwin':  # Mac OS X, assume no X server running
        try:
            os.system('open /Applications/Preview.app %s &' % image_file)
        except:
            raise OSError, "Cannot display %s with Preview application" % \
                  image_file

      elif os_name == 'posix':  # Linux, Unix, etc.
        try:
            os.system('display %s &' % image_file)
        except:
            raise OSError, "Cannot display %s with ImageMagick display. ImageMagick display requires a running X server." % \
                  image_file
      else:
        raise OSError, "no known display function for OS %s" % os_name


# -----------------------------------------------------------------------


def DrawIsingLattice(lattice, scale=0, imsize=800, 
                           imfile=None):
    """DrawIsingLattice(lattice) will draw an image file of the 2D LxL
    square--lattice Ising model, and then will display the result.
    By default, the image file will be stored in a uniquely named png file
    in /tmp, although the image file name can be supplied optionally with
    the imfile argument."""
    # Set up image file
    if imfile is None:
        imfile = tempfile.mktemp()  # make unique filename in /tmp
        imfile += "_square_network_sites.png"
    L = len(lattice) 
    if (scale==0):
        scale = max(1,int(imsize/L)) # Size of squares for each node
    # Background white (in case some nodes missing)
    white = (255,255,255)
    im = Image.new('RGB', (scale*L, scale*L), white)
    if (scale>1):
        draw = ImageDraw.Draw(im)
    black = (0,0,0) # starting color
    # Draw up spins black
    for i in range(L):
        for j in range(L):
            if lattice[i,j]==1:
                if (scale==1):
		    im.putpixel((i,j), black)
                else:
	            x = i*scale
		    y = j*scale
		    draw.rectangle(((x,y),(x+scale,y+scale)), fill=black)
    im.save(imfile)
    Display(imfile)
    return im
