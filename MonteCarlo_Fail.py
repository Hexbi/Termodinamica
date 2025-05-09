import numpy as np
from random import randrange
from vpython import *

#System variables
temperature = 1
particleNumber = 100 
dimensions = 3 
cubeSide = 1
timePerFrame = 1E-5
n_steps = 10000 
gray = color.gray(0.7) # color of edges of container

#Particle Variables
particleRadius = 0.03 # wildly exaggerated size of helium atom
mass = 4E-3/6E23 # helium mass
characteristicTemperature = 33.3 # K
sigma = 1

#Constants
k = 1#1.4E-23 # Boltzmann constant
pi2 = pi * 2

# Derived variables
avgKineticEnergy = sqrt(2*mass*1.5*k*temperature) # average kinetic energy p**2/(2mass) = (3/2)kT
beta = 1/(k*temperature)
avgVelocity = sqrt(8*k*temperature/(pi * mass))
d = cubeSide/2+particleRadius
volume = d**3
distanceDiferential = 0.1*cubeSide#avgVelocity * timePerFrame
epsilon = 1#k * characteristicTemperature

#Control variables
atoms = np.empty(particleNumber, dtype=sphere)
totalTime = 0
loops = 0
valid = True
acceptedMovements = 0
avgEnergy = 0
avgEnergy2 = 0

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

min_dist = 0.9 * sigma
def check_overlap(pos, atoms, particleRadius):
    for atom in atoms:
        if atom is not None:
            distance = mag(pos - atom.pos)
            if distance < 2 * particleRadius:  # Si la distancia entre las partículas es menor que la suma de los radios, hay solapamiento
                return True
        else:
            continue
    return False

if dimensions == 1:
    for i in range(particleNumber):
        while True:
            pos = cubeSide * vector(random() - 0.5, 0, 0)
            if not check_overlap(pos, atoms, particleRadius):
                atoms[i] = sphere(pos=pos, radius=particleRadius, color=gray)
                break
elif dimensions == 2:
    for i in range(particleNumber):
        while True:
            pos = cubeSide * vector(random() - 0.5, random() - 0.5, 0)
            if not check_overlap(pos, atoms, particleRadius):
                atoms[i] = sphere(pos=pos, radius=particleRadius, color=gray)
                break
elif dimensions == 3:
    for i in range(particleNumber):
        while True:
            pos = cubeSide * vector(random() - 0.5, random() - 0.5, random() - 0.5)
            if not check_overlap(pos, atoms, particleRadius):
                atoms[i] = sphere(pos=pos, radius=particleRadius, color=gray)
                break
else:
    valid = False
    print("ValueError: The variable 'dimensions' can only be 1, 2 or 3.")

def RandomDisplacement(dim):
    if dim == 1:
        return distanceDiferential * vector(random.choice([-1, 1]), 0, 0)
    if dim == 2:
        polarAngle = np.random.uniform(0, 2 * np.pi)
        return distanceDiferential * vector(np.cos(polarAngle), np.sin(polarAngle), 0)
    if dim == 3:
        polarAngle =  np.random.uniform(0, 2 * np.pi)
        azimuthalAngle =  np.random.uniform(0, np.pi)
        return distanceDiferential * vector(np.cos(polarAngle) * np.sin(azimuthalAngle), np.sin(polarAngle) * np.sin(azimuthalAngle), np.cos(azimuthalAngle))
    return 0

def EnergyWithoutParticleI(gas, index):
    u = 0
    for i in range(particleNumber):
        if i != index:
            for j in range(i + 1, particleNumber):
                if j != index and i != j:
                    a = (sigma/mag(gas[i].pos - gas[j].pos))**6
                    u += a * (a - 1)
    u *= 4
    return u

def EnergyOnlyParticleI(gas, index):
    u = 0
    for i in range(particleNumber):
        if i != index:
            a = (sigma/mag(gas[i].pos - gas[index].pos))**6
            u += a * (a - 1)
    u *= 4
    return u

def TotalEnergy(gas):
    u = 0.0
    for i in range(particleNumber):
            for j in range(i + 1, particleNumber):
                a = (2.5e-10/mag(gas[i].pos - gas[j].pos))**6
                u += a * (a - 1)
    u *= 4 
    return u

u_i = TotalEnergy(atoms)
u_f = 0
u = TotalEnergy(atoms)
savePos = np.empty(dimensions, dtype=np.float32)
halfUsefullSpace = ((cubeSide*0.5)-particleRadius)
if valid: 
    for n in range(n_steps):

        particle = randrange(particleNumber)  # Choose random particle
        prevPos = vector(atoms[particle].pos.x, atoms[particle].pos.y, atoms[particle].pos.z) # Save initial position

        # Generates and ads random displacement
        dr = RandomDisplacement(dimensions) 
        atoms[particle].pos += dr
        
        # Checks if it's inside the box
        if ((abs(atoms[particle].pos.x) > halfUsefullSpace) or 
            (abs(atoms[particle].pos.y) > halfUsefullSpace) or 
            (abs(atoms[particle].pos.z) > halfUsefullSpace)):
            atoms[particle].pos = vector(prevPos.x, prevPos.y, prevPos.z) #discard movement
            # If not inside the box, revert displacement and try again
            continue

        # Calculates energy before and after the displacement
        atoms[particle].pos += dr
        u_f = TotalEnergy(atoms)


        du = u_f - u_i
        if du <= 0:
            acceptedMovements += 1
            u_i = u_f
        else:
            rnd = np.random.rand()
            expo = np.exp(-du * beta)
            if rnd < expo and n > n_steps * 0.5:
                acceptedMovements += 1
                u_i = u_f
            else:
                atoms[particle].pos = vector(prevPos.x, prevPos.y, prevPos.z) #discard movement
            

        
        # with open("energia_total.txt", "a") as archivo:
        #     archivo.write(f"{u_i},")

        loops +=1
        avgEnergy += u
        avgEnergy2 += u * u
        e_avg = avgEnergy / loops
        e2_avg = avgEnergy2 / loops
        cv = (e2_avg - e_avg ** 2) / (k*temperature)**2
        print('Cv Teoric = ', (dimensions/2)*particleNumber)
        print('Cv Empiric = ', cv)