from vpython import *
import numpy as np
import matplotlib.pyplot as plt
#Web VPython 3.2

# Hard-sphere gas.

# Bruce Sherwood

win = 500

Natoms = 500  # change this to have more or fewer atoms

# Typical values
L = 1 # container is a cube L on a side
gray = color.gray(0.7) # color of edges of container
mass = 4E-3/6E23 # helium mass
Ratom = 0.03 # wildly exaggerated size of helium atom
k = 1.4E-23 # Boltzmann constant
T = 300 # around room temperature
dt = 1E-5




deltav = 100 # binning for v histogram

def checkCollisions():
    hitlist = []
    r2 = 2*Ratom
    r2 *= r2
    for i in range(Natoms):
        ai = apos[i]
        for j in range(i) :
            aj = apos[j]
            dr = ai - aj
            if mag2(dr) < r2: hitlist.append([i,j])
    return hitlist

R0 =  0.001

Radis = []
Col_sec_teo = []
T_avg_teo =[]
Col_sec_sim = []
T_avg_sim =[]
for i in range(10):
    Radis.append(R0+i*R0*0.5)
for i in range(20):
    Radis.append(6*R0*+i*R0*2) 

for j in range(len(Radis)):

    N_colisions = 0
    T_total = 0
    Ratom = Radis[j]
    v_avg=np.sqrt(8*k*T/(mass*np.pi))
    d = L/2+Ratom
    def col_avg(R,D):
        return np.sqrt(2)*(R)**2*np.pi*v_avg*Natoms/(D**3)
    Col_sec_teo.append(col_avg(Ratom,d))
    T_avg_teo.append((col_avg(Ratom,d))**(-1))


 #   Atoms = []
    p = []
    apos = []
    pavg = sqrt(2*mass*1.5*k*T) # average kinetic energy p**2/(2mass) = (3/2)kT
        


    dv = 10

    accum = []
    for i in range(int(3000/deltav)): accum.append([deltav*(i+.5),0])


    d = L/2+Ratom
    r = 0.005
 
    for i in range(Natoms):
        x = L*random()-L/2
        y = L*random()-L/2
        z = L*random()-L/2
 #       if i == 0:
 #           Atoms.append(sphere(pos=vector(x,y,z), radius=Ratom, color=color.cyan, make_trail=True, retain=100, trail_radius=0.3*Ratom))
 #       else: Atoms.append(sphere(pos=vector(x,y,z), radius=Ratom, color=gray))
        apos.append(vec(x,y,z))
        theta = pi*random()
        phi = 2*pi*random()
        px = pavg*sin(theta)*cos(phi)
        py = pavg*sin(theta)*sin(phi)
        pz = pavg*cos(theta)
        p.append(vector(px,py,pz))

    while T_total<0.005:
        rate(300)
        T_total +=dt
        # Accumulate and average histogram snapshots


        # Update all positions
#        for i in range(Natoms): Atoms[i].pos = apos[i] = apos[i] + (p[i]/mass)*dt
        
        # Check for collisions
        hitlist = checkCollisions()
        N_colisions += len(hitlist)
        # If any collisions took place, update momenta of the two atoms
        for ij in hitlist:
            i = ij[0]
            j = ij[1]
            ptot = p[i]+p[j]
            posi = apos[i]
            posj = apos[j]
            vi = p[i]/mass
            vj = p[j]/mass
            vrel = vj-vi
            a = vrel.mag2
            if a == 0: continue;  # exactly same velocities
            rrel = posi-posj
            if rrel.mag > Ratom: continue # one atom went all the way through another
        
            # theta is the angle between vrel and rrel:
            dx = dot(rrel, vrel.hat)       # rrel.mag*cos(theta)
            dy = cross(rrel, vrel.hat).mag # rrel.mag*sin(theta)
            # alpha is the angle of the triangle composed of rrel, path of atom j, and a line
            #   from the center of atom i to the center of atom j where atome j hits atom i:
            alpha = asin(dy/(2*Ratom)) 
            d = (2*Ratom)*cos(alpha)-dx # distance traveled into the atom from first contact
            deltat = d/vrel.mag         # time spent moving from first contact to position inside atom
            
            posi = posi-vi*deltat # back up to contact configuration
            posj = posj-vj*deltat
            mtot = 2*mass
            pcmi = p[i]-ptot*mass/mtot # transform momenta to cm frame
            pcmj = p[j]-ptot*mass/mtot
            rrel = norm(rrel)
            pcmi = pcmi-2*pcmi.dot(rrel)*rrel # bounce in cm frame
            pcmj = pcmj-2*pcmj.dot(rrel)*rrel
            p[i] = pcmi+ptot*mass/mtot # transform momenta back to lab frame
            p[j] = pcmj+ptot*mass/mtot
            apos[i] = posi+(p[i]/mass)*deltat # move forward deltat in time
            apos[j] = posj+(p[j]/mass)*deltat

        for i in range(Natoms):
            loc = apos[i]
            if abs(loc.x) > L/2:
                if loc.x < 0: p[i].x =  abs(p[i].x)
                else: p[i].x =  -abs(p[i].x)
            
            if abs(loc.y) > L/2:
                if loc.y < 0: p[i].y = abs(p[i].y)
                else: p[i].y =  -abs(p[i].y)
            
            if abs(loc.z) > L/2:
                if loc.z < 0: p[i].z =  abs(p[i].z)
                else: p[i].z =  -abs(p[i].z)
    if N_colisions>0:
        Col_sec_sim.append((N_colisions)/T_total/Natoms)  
        T_avg_sim.append(((N_colisions)/T_total/Natoms)**(-1))
    else:
        Col_sec_sim.append(0)  
        T_avg_sim.append(0.07) 
    print(len(Col_sec_sim))
plt.plot(Radis,Col_sec_sim)
plt.plot(Radis,Col_sec_teo)
plt.show()
plt.plot(Radis,T_avg_sim)
plt.plot(Radis,T_avg_teo)
plt.show()