
from vpython import *
import random

# Hard-sphere gas.

# Bruce Sherwood



# Hard-sphere gas.

# Bruce Sherwood

win = 500



Natoms = 500  # change this to have more or fewer atoms
graph2 = graph(title='Pressió en funció del temps', xtitle='Temps', ytitle='Pressió')
# Typical values
L = 1 # container is a cube L on a side
gray = color.gray(0.7) # color of edges of container
mass = 4E-3/6E23 # helium mass
Ratom = 0.03 # wildly exaggerated size of helium atom
k = 1.4E-23 # Boltzmann constant
T = 300 # around room temperature
dt = 5*1E-5


R0 =  0.001
Radis = []

for i in range(10):
    Radis.append(R0+i*0.001*10)


    
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

nhisto = 0 # number of histogram snapshots to average
Volum=[]

for i in range(len(Radis)):
    P_mean=[]
    L_mobil = L/2
    Volum = []
    Volum2 = []
    P_ext = 0.3*1E-17
    P_teo = P_ext
    P_int = []

    Ratom = Radis[i]
   
    animation = canvas( width=1, height=1, align='left')
    animation.range = L
    animation.title = 'Gas'
    s = """  T
      
    """
    animation.caption = s
    
    d = L/2+Ratom
    r = 0.005
    boxbottom = curve(color=gray, radius=r)
    boxbottom.append([vector(-d,-d,-d), vector(-d,-d,d), vector(d,-d,d), vector(d,-d,-d), vector(-d,-d,-d)])
    boxtop = curve(color=gray, radius=r)
    boxtop.append([vector(-d,d,-d), vector(-d,d,d), vector(d,d,d), vector(d,d,-d), vector(-d,d,-d)])
    vert1 = curve(color=gray, radius=r)
    vert2 = curve(color=gray, radius=r)
    vert3 = curve(color=gray, radius=r)
    vert4 = curve(color=gray, radius=r)
    vert1.append([vector(-d,-d,-d), vector(-d,d,-d)])
    vert2.append([vector(-d,-d,d), vector(-d,d,d)])
    vert3.append([vector(d,-d,d), vector(d,d,d)])
    vert4.append([vector(d,-d,-d), vector(d,d,-d)])
    
    Atoms = []
    p = []
    apos = []
    pavg = sqrt(2*mass*1.5*k*T) # average kinetic energy p**2/(2mass) = (3/2)kT
        
    for i in range(Natoms):
        x = L*random.random()-L/2
        y = L*random.random()-L/2
        z = L*random.random()-L/2
        if i == 0:
            Atoms.append(sphere(pos=vector(x,y,z), radius=Ratom, color=color.cyan, make_trail=True, retain=100, trail_radius=0.3*Ratom))
        else: Atoms.append(sphere(pos=vector(x,y,z), radius=Ratom, color=gray))
        apos.append(vec(x,y,z))
        theta = pi*random.random()
        phi = 2*pi*random.random()
        px = pavg*sin(theta)*cos(phi)
        py = pavg*sin(theta)*sin(phi)
        pz = pavg*cos(theta)
        p.append(vector(px,py,pz))
    
    deltav = 100 # binning for v histogram
    
    def barx(v):
        return int(v/deltav) # index into bars array
    
    nhisto = int(4500/deltav)
    histo = []
    for i in range(nhisto): histo.append(0.0)
    histo[barx(pavg/mass)] = Natoms
    
    gg = graph( width=1, height=0.4, xmax=3000, align='left',
        xtitle='speed, m/s', ytitle='Number of atoms', ymax=Natoms*deltav/1000)
    
    theory = gcurve( color=color.blue, width=2 )
    dv = 10
    for v in range(0,3001+dv,dv):  # theoretical prediction
        theory.plot( v, (deltav/dv)*Natoms*4*pi*((mass/(2*pi*k*T))**1.5) *exp(-0.5*mass*(v**2)/(k*T))*(v**2)*dv )
    
    accum = []
    for i in range(int(3000/deltav)): accum.append([deltav*(i+.5),0])
    vdist = gvbars(color=color.red, delta=deltav )
    
    def interchange(v1, v2):  # remove from v1 bar, add to v2 bar
        barx1 = barx(v1)
        barx2 = barx(v2)
        if barx1 == barx2:  return
        if barx1 >= len(histo) or barx2 >= len(histo): return
        histo[barx1] -= 1
        histo[barx2] += 1
        

        
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
    
    nhisto = 0 # number of histogram snapshots to average
    T_total=0
    Temps2= []
    Temps=[]
    while T_total<0.1:
        rate(300)
        I = 0
        T_total += dt
        # Accumulate and average histogram snapshots
        for i in range(len(accum)): accum[i][1] = (nhisto*accum[i][1] + histo[i])/(nhisto+1)
        if nhisto % 10 == 0:
            vdist.data = accum
        nhisto += 1
    
        # Update all positions
        for i in range(Natoms): Atoms[i].pos = apos[i] = apos[i] + (p[i]/mass)*dt
        
        # Check for collisions
        hitlist = checkCollisions()
    
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
            interchange(vi.mag, p[i].mag/mass)
            interchange(vj.mag, p[j].mag/mass)
        
        #This checks for colisions against the walls
   #     if len(P_int)>0 and P_ext>P_int[len(P_int)-1]:
#            L_mobil= L_mobil-dx
#        else len(P_int)>0 and P_ext<P_int[len(P_int)-1]:
#            L_mobil= L_mobil+dx
        
        
        for i in range(Natoms):
            loc = apos[i]
            if loc.x > L_mobil:#Aquest terme s'aegeix per afegir el canvi de moment a causa de la col·lisió
                p[i].x =  -abs(p[i].x)
                I += 2* abs(p[i].x)
            
            if loc.x < -L/2:
                p[i].x =  abs(p[i].x)
                I += 2* abs(p[i].x) #Aquest terme s'aegeix per afegir el canvi de moment a causa de la col·lisió            
            
            if abs(loc.y) > L/2:
                if loc.y < 0: 
                    p[i].y = abs(p[i].y)
                    I += 2* abs(p[i].y)
                else: 
                    p[i].y =  -abs(p[i].y)
                    I += 2* abs(p[i].y)
            
            if abs(loc.z) > L/2:
                if loc.z < 0: 
                    p[i].z =  abs(p[i].z)
                    I += 2* abs(p[i].z)
                else: 
                    p[i].z =  -abs(p[i].z)
                    I += 2* abs(p[i].z)
        
        A = (2*L**2+4*(L_mobil+L/2)*L)
        P_int.append(I/(A*dt))
    
        if len(P_int)==50:
            Temps2.append(T_total)
            sum=0
            for i in range(10):
                sum+=P_int[i]
            
            P_mean.append(sum/10)
            P_int=[]
        
        dx2=0.0005*L
        if len(P_int)>0 and P_ext>P_int[len(P_int)-1]:
            L_mobil= L_mobil-dx2
        else:
            L_mobil= L_mobil+dx2
            
        Temps.append(T_total)
        Volum.append((L_mobil+1/2))
        #print(Volum)
    f1 = gcurve(graph=graph2, color=vec(random.random(),random.random(),random.random()),label = f"R={round(Ratom, 4)}")
    for x, y in zip(Temps2, P_mean):
        f1.plot(x, y)
    #plt.plot(Temps,P_int, label = f"R={round(Ratom,4)}")

    Volum2.append(Volum[len(Volum)-1])
    print(Volum[len(Volum)-1])

P_teo_list = []
for i in range(len(Temps)):
    P_teo_list.append(P_teo)



#Create a curve to plot points
f3 = gcurve(graph=graph2,color=color.red)


#Plot the points
for x, y in zip(Temps, P_teo_list):
    f3.plot(x, y)
    
graph3 = graph(title='Volum final en funció del radi', xtitle='Radi (m)', ytitle='Volum (m^3)')
f4 = gcurve(graph=graph3,color=color.red)
for x, y in zip(Radis, Volum2):
    f4.plot(x, y)


#plt.plot(Temps,P_teo_list,label="Pressió teórica")

#plt.grid(True)
#plt.legend(fontsize="small")
#plt.title("Pressió en funció del temps per diferents radis")
#plt.show()

