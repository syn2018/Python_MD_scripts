"""
DynamicNetwork
"""

import numpy, math
import Tkinter 
import Image, ImageTk, ImageDraw

DefaultImageSize = (400,400)

root = Tkinter.Tk()
root.withdraw()

class DynamicNetwork (Tkinter.Label):

    def __init__(self, size=DefaultImageSize, mode='P', 
                 zmin=0.0, zmax=1.0):
        top = Tkinter.Toplevel()
        top.title('DynamicNetwork')
        Tkinter.Label.__init__(self, top)
	self.size = size
	self.mode = mode
	self.zmin = zmin
	self.zmax = zmax
	self.zrange = zmax-zmin
        self.canvas = Tkinter.Canvas(top, 
				     width=self.size[0], height=self.size[1])
	self.im = Image.new('RGB', self.size)
        self.blank = self.im.copy()
        self.draw = ImageDraw.Draw(self.im)
	if self.mode == '1':
	    self.tkimage = \
			 ImageTk.BitmapImage(self.im,
					     foreground="white")
	else:
	    self.tkimage = \
			 ImageTk.PhotoImage(self.im)
	self.canvas.create_image(0, 0, anchor=Tkinter.NW, image=self.tkimage)
        self.canvas.pack()

    def reset(self):
        self.im.paste(self.blank)

    def setTitle(self, title):
        self.master.title(title)

    def displayFromLists(self, all_nodes, edgelist, nodec={}, edgec={},
                         clear=False):
        if clear:
            self.reset()
        windowMargin = 0.1
        dotsize=4
        windowSize = self.size[0]
        radius = (1.-2*windowMargin)*windowSize/2.
        color = (255,255,255)
        center = windowSize/2.
        L = len(all_nodes)
        # Create a dictionary that maps nodes to their positions
        # around the circle
        nodePosition = {}
        for index, node in enumerate(all_nodes):
            theta = 2.*numpy.pi*float(index)/L
            x = radius * numpy.cos(theta) + center
            y = -radius * numpy.sin(theta) + center
            nodePosition[node] = (x, y)
            nc = nodec.get(node, color)
            self.draw.ellipse( ((x-dotsize/2, y-dotsize/2),
                                (x+dotsize/2, y+dotsize/2)), fill=nc )

        # Draw the lines between the nodes
        for node, neighbor in edgelist:
            # We want to draw bonds only once, even though two bonds connect
            #   node and neighbor.  We can test "if neighbor > node" to
            #   implement this.  For any pair of objects, Python will
            #   consistently define this operation (e.g., if x1<x2 is True,
            #   then x1>x2 is False).  For arbitrary node IDs, we may not
            #   know what it means for ID1 to be greater than ID2, but it
            #   does not matter for the purpose here.  For Python classes,
            #   one define the special __cmp__ (compare) method that
            #   indicates how to compare (>,<,==) two instances of the
            #   class.
            ec = edgec.get((node, neighbor), color)
            self.draw.line((nodePosition[node],
                            nodePosition[neighbor]),
                           fill=ec)
        self.tkimage.paste(self.im)
        self.canvas.update()

	    
