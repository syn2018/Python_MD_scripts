import DynamicLattice
import RandomArray
    
def test(shape=(100,100)):
    dl = DynamicLattice.DynamicLattice(shape)
    for n in range(20):
        a = RandomArray.randint(0, 2, shape)
        dl.display(a)
    

def test2(shape=(100,100)):
    dl = DynamicLattice.DynamicLattice(shape)
    a = RandomArray.randint(0, 2, shape)
    dl.display(a)
    for i in range(shape[0]/2):
        for j in range(shape[0]/2):
            a[i,j] = 0
            dl.display(a, (i,j))
            
