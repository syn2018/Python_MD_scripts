"""
PlotAtoms
"""

import numpy, math, time
import Tkinter 
import Image, ImageTk, ImageDraw

DefaultImageSize = (400,400)

root = Tkinter.Tk()
root.withdraw()

class PlotAtoms (Tkinter.Label):

    def __init__(self, scale, size=DefaultImageSize, mode='P'):
        top = Tkinter.Toplevel()
        top.title('PlotAtoms')
        Tkinter.Label.__init__(self, top)
        self.scale = scale
	self.size = size
	self.mode = mode
        self.canvas = Tkinter.Canvas(top, 
				     width=self.size[0], height=self.size[1])
	self.im = Image.new('RGB', self.size)
        #self.displayIm = self.im.resize(self.size)
        self.blank = self.im.copy()
        self.draw = ImageDraw.Draw(self.im)
        #self.tkimage = ImageTk.PhotoImage(self.displayIm)
        self.tkimage = ImageTk.PhotoImage(self.im)
	self.canvas.create_image(0, 0, anchor=Tkinter.NW, image=self.tkimage)
        self.canvas.pack()
        self.counter = 0

    def reset(self):
        self.im.paste(self.blank)

    def setTitle(self, title):
        self.master.title(title)

    def display(self, positions, colors=None, sizes=None,
                pause=None, clear=False, saveFrames=False):
        if pause is not None:
            time.sleep(pause)
        if clear:
            self.reset()
        windowSize = self.size[0]
        color = (255,255,255)
        for id, pos in enumerate(positions):
            x,y = pos
            dotsize = sizes[id]
            self.draw.ellipse( ((x-dotsize/2, y-dotsize/2), 
                                (x+dotsize/2, y+dotsize/2)), fill=colors[id] )
        #self.displayIm = self.im.resize(self.size)
        #self.tkimage.paste(self.displayIm)
        self.tkimage.paste(self.im)
        if saveFrames:
            self.im.save('tmp/frame%06d.gif' % self.counter)
        self.counter += 1
        self.canvas.update()

	    
def test(num_nodes = 1000):
    import numpy
    dc = PlotAtoms((400,400))
    #colortab = [(200, 20, 20), (20, 200, 20), (20, 20, 200)]
    colortab = ['blue', 'red', 'yellow']
    pos = 400.*numpy.random.random((num_nodes, 2))    
    sizes = [2*((i%3)+1) for i in range(num_nodes)]
    #colors = [(i%255, i%255, i%255) for i in range(num_nodes)]
    colors = [colortab[i%3] for i in range(num_nodes)]
    for i in range(100):
        pos = pos + (10.*numpy.random.random((num_nodes, 2))-5.)
        dc.display(pos, colors, sizes, pause=None, clear=True, saveFrames=False)
    return dc



    
