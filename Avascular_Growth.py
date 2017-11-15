# Final version of the Tumoral Growth Model (without genetic network)
# authors = 'Charles GAYDON - Helio WANG'
# Last edited : 06/2016

# USE :
# Python 2.7

import random
import numpy as np
import copy

PROLIFERATING = 0
NON_TUMOR = -1
QUIESCENT = 1
NECROTIC = 2

diffusion_constant = [0.0594, 0.00152, 0.002124,
                      10 ** (-6), 10 ** -6]
co = 1
metabolic_rate = {PROLIFERATING: np.array([-108, -162, 240, 1, 0])/co,
                  QUIESCENT: np.array([-50, -80, 110, 0.5, 1])/co,
                  NECROTIC: np.array([0, 0, 0, 0, 2])/co}

n_molecules = 3 # Those are Oxygen, Glucose and Wastes


def laplacian(Z, dx):
    Ztop = Z[0:-2, 1:-1]
    Zleft = Z[1:-1, 0:-2]
    Zbottom = Z[2:, 1:-1]
    Zright = Z[1:-1, 2:]
    Zcenter = Z[1:-1, 1:-1]
    return (Ztop + Zleft + Zbottom + Zright - 4 * Zcenter) / dx**2


class Site:
    def __init__(self, i, j, cell_id, cell_type):
        self.i = i
        self.j = j
        self.cell_id = cell_id
        self.cell_type = cell_type
        if self.cell_type == NON_TUMOR:
            self.concentration = [0.28, 5.5, 0.0]  # , 1.0, 0.0]
        else:
            self.concentration = [0.0, 0.0, 0.0]  # , 0.0, 0.0]

    def invaded_by(self, k):
        self.cell_id = k
        self.cell_type = PROLIFERATING


class Cell:
    def __init__(self, sites):
        self.cell_type = PROLIFERATING
        self.sites = sites

    def add_site(self, site):
        self.sites.add(site)

    def update_state(self, grid, radius):
        concentration = [0] * n_molecules
        grid2 = copy.copy(grid)
        rad = 0
        n = len(grid)
        for i, j in self.sites:
            if np.abs(i-n/2) > rad:
                rad = np.abs(i-n/2)
            if np.abs(j-n/2) > rad:
                rad = np.abs(j-n/2)
            site = grid[i][j]
            for u in [-1,1,0] :
                for v in [-1,1,0]:
                    vois = grid[(i+u)%size][(j+v)%size]
                    for k in xrange(n_molecules):
                        concentration[k] += vois.concentration[k]
            for k in xrange(n_molecules):
                concentration[k] /= float(9)
            site = grid2[i][j]
            if site.cell_type == QUIESCENT and concentration[0] < 0.02 and concentration[1] < 0.06 and concentration[2] > 8:
                p = 1
                if random.uniform(0, 1) < p:
                    grid2[i][j].cell_type = NECROTIC
                site.cell_type = NECROTIC
            elif site.cell_type == PROLIFERATING and concentration[0] < 0.02 and concentration[1] < 0.06 and concentration[2] > 8 and rad < radius - 2:
                p = 1
                if random.uniform(0,1)<p:
                    grid2[i][j].cell_type = QUIESCENT
                    site.cell_type = QUIESCENT
        grid = copy.copy(grid2)


class Tumor:
    def __init__(self, size):
        self.n = size
        self.grid = []
        for i in xrange(size):
            self.grid.append([])
            for j in xrange(size):
                self.grid[i].append(Site(i, j, -1, NON_TUMOR))
        self.grid[size/2][size/2] = Site(0, 0, 0, PROLIFERATING)
        self.cell = {0: Cell({(size/2, size/2)})}
        self.n_cells = 1
        self.radius = 1

    def growth(self):
		# Aggressivenes of proliferating cells
        p = 0.3
        # Parameter of celerity of the proliferation
        vit = 0.55
        rad = self.radius
        hits = 0
        Lu = [-1, 0, 1] 
        Lv = [-1, 0, 1]
        
        for k in xrange(self.n_cells):
            sites = set()
            SSS = list(self.cell[k].sites)
            random.shuffle(SSS)
            lim = 1.0 + vit*float(len(SSS))
            for i, j in SSS[0:int(lim)]:
                if self.grid[i][j].cell_type == PROLIFERATING:
                    random.shuffle(Lu)
                    ok = True
                    for u in Lu:
                        random.shuffle(Lv)
                        for v in Lv:
                            x = (i + u) % self.n
                            y = (j + v) % self.n
                            
                            if self.grid[x][y].cell_type == NON_TUMOR and ok:
                                if random.uniform(0, 1) < p:
                                    self.grid[x][y].invaded_by(self.n_cells+hits)
                                    sites.add((x, y))
                                    ok = False # a cell can reproduce only once per round
                                if x - self.n/2 > rad:
                                    rad = x - self.n/2
                                if y - self.n/2 > rad:
                                    rad = y - self.n/2
            if len(sites) > 0:
                self.cell[self.n_cells+hits] = Cell(sites)
                hits += 1
        self.n_cells += hits
        self.radius = rad

    def numerous(self):
		nb = [0,0,0,0]
		for i in xrange(size):
			for j in xrange(size):
				nb[self.grid[i][j].cell_type]+=1
		for i in [0,1,2,3] :
			nb[i] = float(nb[i])/(float(size*size))
		return nb

    def update_state(self):
        for k in self.cell:
            self.cell[k].update_state(self.grid, self.radius)
	
    def reaction_diffusion(self):
        oxy = np.zeros((self.n, self.n))
        glu = np.zeros((self.n, self.n))
        lac = np.zeros((self.n, self.n))
        molecules = [oxy, glu, lac]
        boundary_conditions = (0.28, 5.5, 0.0)
        for _ in xrange(n):
            for i in xrange(self.n):
                for j in xrange(self.n):
                    for k in xrange(n_molecules):
                        molecules[k][i, j] = self.grid[i][j].concentration[k]

            delta_oxy = laplacian(oxy, dx)
            delta_glu = laplacian(glu, dx)
            delta_lac = laplacian(lac, dx)
            deltas = [delta_oxy, delta_glu, delta_lac]

            for i in xrange(1, self.n-1):
                for j in xrange(1, self.n-1):
                    cell_type = self.grid[i][j].cell_type
                    if cell_type == NON_TUMOR:
                        oxy[i][j], glu[i][j], lac[i][j] = boundary_conditions
                    else:
                        for k in xrange(n_molecules):
                            molecules[k][i, j] += dt * (diffusion_constant[k] * deltas[k][i-1, j-1]
                                                        + metabolic_rate[cell_type][k])

            for i in xrange(1, self.n-1):
                for j in xrange(1, self.n-1):
                    for k in xrange(n_molecules):
                        self.grid[i][j].concentration[k] = max(molecules[k][i, j], 0)

            for i in xrange(self.n):
                self.grid[i][0].concentration = list(boundary_conditions)
                self.grid[i][self.n-1].concentration = list(boundary_conditions)
            for j in xrange(self.n):
                self.grid[0][j].concentration = list(boundary_conditions)
                self.grid[self.n-1][j].concentration = list(boundary_conditions)

    def divide(self):
        pass

########################################################
## PArameters of the Simulation
########################################################

size = 15
dx = 1./size
T = 0.75
dt = .9 * dx**2/2
n = int(T/dt)
Tc = 3
turns = 5

name = 'p'+str(0.9)+'size'+str(size)+'co'+str(co)+'vit'+str(0.05)+'Tc'+str(Tc)+'.txt'
print("Name of the file where results shall be printed is : " + name)

## Simulation

tum = Tumor(size)
oxy = []
nb = []
nb.append(tum.numerous())
for t in range(turns):
    print(""+str(t)+"/" + str(turns))
    oxy.append([])
    nb.append(tum.numerous())
    print("[Cells prop] ="+ str(map(lambda x : round(x*100,1),nb[-1])))
    for i in xrange(size):
        for j in xrange(size):
            # oxy[t].append(tum.grid[i][j].concentration[0])
            oxy[t].append(tum.grid[i][j].cell_type)
    oxy[t] = ' '.join(map(str, oxy[t]))

    for i in xrange(Tc):
        tum.growth()
    tum.reaction_diffusion()
    tum.update_state()

print("End of simulation.")
f = open('History/'+name, 'w')
f.write('\n'.join(oxy))
f.close()

infile = 'History/Concentrations_'+name
with file(infile, 'w') as outfile: 
	np.savetxt(outfile,nb,fmt='%-7.2f')




