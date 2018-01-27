import visual as v
import math
v.scene.autoscale = 1
v.scene.center = (1, 1, 0)
v.scene.height=700
v.scene.width=500
v.scene.forward=(0,1,-1)
v.scene.range=3.75

eps = 1.0e-08

PROTEIN_FACTOR = 250.
RNA_FACTOR = 10.
PROMOTER_FACTOR = 1.0
SCALE =  0.031622  # SCALE * sqrt(RNA_FACTOR) = 0.1 for display

class GeneDisplay:

    def __init__(self, index, name, colors, mRNA, protein, P0, P1, P2):
        self.index = index
        self.name = name
        self.pos = (index, 0, 0)
        self.color = colors[index]
        self.mRNA_amt = mRNA
        self.protein_amt = protein
        self.P0_amt = P0
        self.P1_amt = P1
        self.P2_amt = P2
        self.f = v.frame()
        j = (index+1)%3
        self.P0 = v.cylinder(display=v.scene,
                             frame=self.f,
                             pos=(-0.1,-1,0), axis=(0,0,1),
                             radius=SCALE*math.sqrt(PROMOTER_FACTOR),
                             length=self.P0_amt*math.sqrt(PROMOTER_FACTOR),
                             color=colors[j])
        self.P1 = v.cylinder(display=v.scene,frame=self.f,
                            pos=(0.0,-1,0), axis=(0,0,1),
                            radius=SCALE*math.sqrt(PROMOTER_FACTOR),
                            length=self.P1_amt*math.sqrt(PROMOTER_FACTOR),
                            color=colors[j])
        self.P2 = v.cylinder(display=v.scene,frame=self.f,
                            pos=(0.1,-1,0), axis=(0,0,1),
                            radius=SCALE*math.sqrt(PROMOTER_FACTOR),
                            length=self.P2_amt*math.sqrt(PROMOTER_FACTOR),
                            color=colors[j])
        self.m = v.cylinder(display=v.scene,frame=self.f, pos=(0,0,0), axis=(0,0,1),
                       radius = SCALE*math.sqrt(RNA_FACTOR),
                       length=self.mRNA_amt/RNA_FACTOR,
                       color=self.color)
        self.p = v.cylinder(display=v.scene,frame=self.f, pos=(0,1,0), axis=(0,0,1),
                       radius = SCALE*math.sqrt(PROTEIN_FACTOR),
                       length=self.protein_amt/PROTEIN_FACTOR,
                       color=self.color)
        self.lab = v.label(frame=self.f, text=self.name, pos=(0,-1.5,0),
                           height=100, color=self.color,
                           box=0, font='helvetica')
        self.f.pos = self.pos
        self.f.axis = (1, 0, 0)


class RepressilatorAnimator:

    def __init__(self, trajectory=None):
        for obj in v.scene.objects:
            obj.visible = 0
        self.trajectory = trajectory
        self.colors = [v.color.red, v.color.yellow, v.color.green]
        self.names = ['lacI', 'tetR', 'cI']
        self.displays = []
        for i in range(3):
            gd = GeneDisplay(i, self.names[i],
                             self.colors, eps, eps, eps, eps, eps)
            self.displays.append(gd)
            v.scene.autocenter = 0
            v.scene.autoscale = 0
        if trajectory is not None:
            self.AnimateTrajectory(self.trajectory)
            
    def AnimateTrajectory(self, trajectory):
        for point in trajectory:
            self.mRNAs = point[0:3]
            self.proteins = point[3:6]
            self.P0 = point[6:9]
            self.P1 = point[9:12]
            self.P2 = point[12:15]
            self.Draw()

    def Draw(self):
        v.rate(100)
        for i in range(3):
            self.displays[i].m.length = (self.mRNAs[i]+eps)/RNA_FACTOR
            self.displays[i].p.length = (self.proteins[i]+eps)/PROTEIN_FACTOR
            self.displays[i].P0.length = (self.P0[i] + eps)/\
                                         PROMOTER_FACTOR
            self.displays[i].P1.length = (self.P1[i] + eps)/\
                                         PROMOTER_FACTOR
            self.displays[i].P2.length = (self.P2[i] + eps)/\
                                         PROMOTER_FACTOR


def ClearScene():
    for obj in v.scene.objects:
        v.scene.objects.remove(obj)
        
