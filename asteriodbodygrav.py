# -*- coding: utf-8 -*-

import numpy as np

PI = np.pi

dataFile =  open('gravDataFull.dat', 'w')

t = 0.0
tHalf = 0.0
dt = 0.005 #years

kinEnergy1=0.0
kinEnergy2=0.0
potEnergy = 0.0

class asteroid:
    def __init__(self, info, mass):
        self.now = np.zeros(2)
        self.nextstep = np.zeros(2)
        self.prev = np.zeros(2) #previous step (i-1)
        self.vel = np.zeros(2) #velocities
        
        self.now[0] = info[0] #x position
        self.now[1] = info[1] #y position
        self.vel[0] = info[2] #x velocity
        self.vel[1] = info[3] #y velocity
        self.distanceJ = 0
        self.distanceS =0
        self.speed = 0
        self.kinetic = 0
        self.potential = 0
        self.momentum = 0
        self.AngMom = 0
        self.maxDist = info[0]
        self.minDist = info[0]
        self.mass = mass
        
def initVel(a):
    vel = np.sqrt((4*pow(PI,2))/(a))
    return vel

#dataFile.write('Time X1 Y2 VX1 VY1 X2 Y2 VX2 VY2')
#sun is at the origin
def findDistances(body):
    body.distanceJ = pow( pow(body.now[0]-Jupiter.now[0], 2) + pow(body.now[1] - Jupiter.now[1], 2), 0.5)
    body.distanceS = pow( pow(body.now[0]-Sun[0], 2) + pow(body.now[1] - Sun[1], 2), 0.5)
    if body.distanceS >= body.maxDist:
        body.maxDist = body.distanceS
    if body.distanceS <= body.minDist:
        body.minDist = body.distanceS
    
    
def findSpeed(body):
    body.speed = pow( pow(body.vel[0],2) + pow(body.vel[1],2), 0.5)
    
def findKinetic(body):
    body.kinetic = body.mass*pow(body.speed,2)/2
    
def findKineticJ(body):
    body.kinetic = body.mass*pow(body.speed,2)/2
    
def findPotentialJ(body):
    body.potential = -4*pow(PI,2)*body.mass/body.distanceS
    
def findPotential(body):
    body.potential = -4*pow(PI,2)*body.mass/body.distanceS - 4*pow(PI,2)*body.mass/body.distanceJ
    
def findMomentum(body):
    body.momentum = body.mass*body.speed
    
def findAngMom(body):
    body.AngMom = body.mass*body.speed*body.distanceS
    
def writeBodyToFile(body):
    #body.now.tofile(dataFile, ' ')
    #dataFile.write(' ')
    #body.vel.tofile(dataFile, ' ')
    #dataFile.write(' ')
    dataFile.write(str(body.distanceS))
    dataFile.write(' ')
    #dataFile.write(str(body.kinetic))
    #dataFile.write(' ')
    #dataFile.write(str(body.potential))
    #dataFile.write(' ')
    #dataFile.write(str(body.momentum))
    #dataFile.write(' ')
    #dataFile.write(str(body.AngMom))
    #dataFile.write(' ')
    dataFile.write(str(body.maxDist))
    dataFile.write(' ')
    dataFile.write(str(body.minDist))
    dataFile.write(' ')
    
def verlet(body):
    body.nextstep[0] = 2*body.now[0] - body.prev[0] + (-4*pow(PI,2)*(body.now[0]- Sun[0]) /pow(body.distanceS, 3) - 4*pow(PI,2)*massRatio*(body.now[0]-Jupiter.now[0])/pow(body.distanceJ,3))*pow(dt,2)
    body.nextstep[1] = 2*body.now[1] - body.prev[1] + (-4*pow(PI,2)*(body.now[1]-Sun[1]) /pow(body.distanceS,3) - 4*pow(PI,2)*massRatio*(body.now[1]-Jupiter.now[1])/pow(body.distanceJ,3))*pow(dt,2)
    body.vel[0] = (body.nextstep[0] - body.prev[0])/(2*dt)
    body.vel[1] = (body.nextstep[1] - body.prev[1])/(2*dt)
    
    
def verletJ(body):
    body.nextstep[0] = 2*body.now[0] - body.prev[0] + (-4*pow(PI,2)*(body.now[0] - Sun[0])/pow(body.distanceS, 3))*pow(dt,2)
    body.nextstep[1] = 2*body.now[1] - body.prev[1] + (-4*pow(PI,2)*(body.now[1]- Sun[1])/pow(body.distanceS, 3))*pow(dt,2)
    body.vel[0] = (body.nextstep[0] - body.prev[0])/(2*dt)
    body.vel[1] = (body.nextstep[1] - body.prev[1])/(2*dt)

def updatebody(body):
    body.prev[0] = body.now[0]
    body.prev[1] = body.now[1]
    body.now[0] = body.nextstep[0]
    body.now[1] = body.nextstep[1]
    
def Euler1(body):
    findDistances(body)
    body.vel[0] = body.vel[0] + (-4*pow(PI,2)*(body.now[0]- Sun[0]) /pow(body.distanceS, 3) - 4*pow(PI,2)*massRatio*(body.now[0]-Jupiter.now[0])/pow(body.distanceJ,3))*dt
    body.vel[1] = body.vel[1] + (-4*pow(PI,2)*(body.now[1]- Sun[1]) /pow(body.distanceS, 3) - 4*pow(PI,2)*massRatio*(body.now[1]-Jupiter.now[1])/pow(body.distanceJ,3))*dt
    
    body.nextstep[0] = body.now[0] + body.vel[0]*dt
    body.nextstep[1] = body.now[1] + body.vel[1]*dt
    
    updatebody(body)
    
def Euler1J(body):
    findDistances(body)
    body.vel[0] = body.vel[0] + (-4*pow(PI,2)*(body.now[0]- Sun[0]) /pow(body.distanceS, 3))*dt
    body.vel[1] = body.vel[1] + (-4*pow(PI,2)*(body.now[1]- Sun[1]) /pow(body.distanceS, 3))*dt
    
    body.nextstep[0] = body.now[0] + body.vel[0]*dt
    body.nextstep[1] = body.now[1] + body.vel[1]*dt
    
    updatebody(body)
    
    
#begin the model:
#scale all masses by the sun mass; energies are then kg*AU^2/yr^2 per sun mass
sunMass = 1.989*pow(10, 30)
x = np.linspace(3.0, 3.6, 40)
body1 = asteroid([x[0], 0, 0, initVel(x[0])], pow(10, 18)/sunMass)
body2 = asteroid([x[1], 0, 0, initVel(x[1])], pow(10, 18)/sunMass)
body3 = asteroid([x[2], 0, 0, initVel(x[2])], pow(10, 18)/sunMass)
body4 = asteroid([x[3], 0, 0, initVel(x[3])], pow(10, 18)/sunMass)
body5 = asteroid([x[4], 0, 0, initVel(x[4])], pow(10, 18)/sunMass)
body6 = asteroid([x[5], 0, 0, initVel(x[5])], pow(10, 18)/sunMass)
body7 = asteroid([x[6], 0, 0, initVel(x[6])], pow(10, 18)/sunMass)
body8 = asteroid([x[7], 0, 0, initVel(x[7])], pow(10, 18)/sunMass)
body9 = asteroid([x[8], 0, 0, initVel(x[8])], pow(10, 18)/sunMass)
body10 = asteroid([x[9], 0, 0, initVel(x[9])], pow(10, 18)/sunMass)

body11 = asteroid([x[10], 0, 0, initVel(x[10])], pow(10, 18)/sunMass)
body12 = asteroid([x[11], 0, 0, initVel(x[11])], pow(10, 18)/sunMass)
body13 = asteroid([x[12], 0, 0, initVel(x[12])], pow(10, 18)/sunMass)
body14 = asteroid([x[13], 0, 0, initVel(x[13])], pow(10, 18)/sunMass)
body15 = asteroid([x[14], 0, 0, initVel(x[14])], pow(10, 18)/sunMass)
body16 = asteroid([x[15], 0, 0, initVel(x[15])], pow(10, 18)/sunMass)
body17 = asteroid([x[16], 0, 0, initVel(x[16])], pow(10, 18)/sunMass)
body18 = asteroid([x[17], 0, 0, initVel(x[17])], pow(10, 18)/sunMass)
body19 = asteroid([x[18], 0, 0, initVel(x[18])], pow(10, 18)/sunMass)
body20 = asteroid([x[19], 0, 0, initVel(x[19])], pow(10, 18)/sunMass)

body21 = asteroid([x[20], 0, 0, initVel(x[20])], pow(10, 18)/sunMass)
body22 = asteroid([x[21], 0, 0, initVel(x[21])], pow(10, 18)/sunMass)
body23 = asteroid([x[22], 0, 0, initVel(x[22])], pow(10, 18)/sunMass)
body24 = asteroid([x[23], 0, 0, initVel(x[23])], pow(10, 18)/sunMass)
body25 = asteroid([x[24], 0, 0, initVel(x[24])], pow(10, 18)/sunMass)
body26 = asteroid([x[25], 0, 0, initVel(x[25])], pow(10, 18)/sunMass)
body27 = asteroid([x[26], 0, 0, initVel(x[26])], pow(10, 18)/sunMass)
body28 = asteroid([x[27], 0, 0, initVel(x[27])], pow(10, 18)/sunMass)
body29 = asteroid([x[28], 0, 0, initVel(x[28])], pow(10, 18)/sunMass)
body30 = asteroid([x[29], 0, 0, initVel(x[29])], pow(10, 18)/sunMass)

body31 = asteroid([x[30], 0, 0, initVel(x[30])], pow(10, 18)/sunMass)
body32 = asteroid([x[31], 0, 0, initVel(x[31])], pow(10, 18)/sunMass)
body33 = asteroid([x[32], 0, 0, initVel(x[32])], pow(10, 18)/sunMass)
body34 = asteroid([x[33], 0, 0, initVel(x[33])], pow(10, 18)/sunMass)
body35 = asteroid([x[34], 0, 0, initVel(x[34])], pow(10, 18)/sunMass)
body36 = asteroid([x[35], 0, 0, initVel(x[35])], pow(10, 18)/sunMass)
body37 = asteroid([x[36], 0, 0, initVel(x[36])], pow(10, 18)/sunMass)
body38 = asteroid([x[37], 0, 0, initVel(x[37])], pow(10, 18)/sunMass)
body39 = asteroid([x[38], 0, 0, initVel(x[38])], pow(10, 18)/sunMass)
body40 = asteroid([x[39], 0, 0, initVel(x[39])], pow(10, 18)/sunMass)
#x, y, vx, vy
Jupiter = asteroid([5.2, 0, 0, initVel(5.2)], 1.9*pow(10, 27)/sunMass)
massRatio = (1.9*pow(10, 27))/(1.989*pow(10, 30))  #mJ/mS
Sun = [0, 0, 0, 0, 1]

#run Euler first for some small time steps
Euler1J(Jupiter)
Euler1(body1)
Euler1(body2)
Euler1(body3)
Euler1(body4)
Euler1(body5)
Euler1(body6)
Euler1(body7)
Euler1(body8)
Euler1(body9)
Euler1(body10)
                   
Euler1(body11)
Euler1(body12)
Euler1(body13)
Euler1(body14)
Euler1(body15)
Euler1(body16)
Euler1(body17)
Euler1(body18)
Euler1(body19)
Euler1(body20)
                   
Euler1(body21)
Euler1(body22)
Euler1(body23)
Euler1(body24)
Euler1(body25)
Euler1(body26)
Euler1(body27)
Euler1(body28)
Euler1(body29)
Euler1(body30)

Euler1(body31)
Euler1(body32)
Euler1(body33)
Euler1(body34)
Euler1(body35)
Euler1(body36)
Euler1(body37)
Euler1(body38)
Euler1(body39)
Euler1(body40)
t = t+dt

for ii in range(1, 25000):
    #  Calculate current values of separation, speeds, and energies.
    findDistances(Jupiter)
    findDistances(body1)
    findDistances(body2)
    findDistances(body3)
    findDistances(body4)
    findDistances(body5)
    findDistances(body6)
    findDistances(body7)
    findDistances(body8)
    findDistances(body9)
    findDistances(body10)
                   
    findDistances(body11)
    findDistances(body12)
    findDistances(body13)
    findDistances(body14)
    findDistances(body15)
    findDistances(body16)
    findDistances(body17)
    findDistances(body18)
    findDistances(body19)
    findDistances(body20)
                   
    findDistances(body21)
    findDistances(body22)
    findDistances(body23)
    findDistances(body24)
    findDistances(body25)
    findDistances(body26)
    findDistances(body27)
    findDistances(body28)
    findDistances(body29)
    findDistances(body30)
    
    findDistances(body31)
    findDistances(body32)
    findDistances(body33)
    findDistances(body34)
    findDistances(body35)
    findDistances(body36)
    findDistances(body37)
    findDistances(body38)
    findDistances(body39)
    findDistances(body40)
    
    findSpeed(Jupiter)
    findSpeed(body1)
    findSpeed(body2)
    findSpeed(body3)
    findSpeed(body4)
    findSpeed(body5)
    findSpeed(body6)
    findSpeed(body7)
    findSpeed(body8)
    findSpeed(body9)
    findSpeed(body10)
                   
    findSpeed(body11)
    findSpeed(body12)
    findSpeed(body13)
    findSpeed(body14)
    findSpeed(body15)
    findSpeed(body16)
    findSpeed(body17)
    findSpeed(body18)
    findSpeed(body19)
    findSpeed(body20)
                   
    findSpeed(body21)
    findSpeed(body22)
    findSpeed(body23)
    findSpeed(body24)
    findSpeed(body25)
    findSpeed(body26)
    findSpeed(body27)
    findSpeed(body28)
    findSpeed(body29)
    findSpeed(body30)
    
    findSpeed(body31)
    findSpeed(body32)
    findSpeed(body33)
    findSpeed(body34)
    findSpeed(body35)
    findSpeed(body36)
    findSpeed(body37)
    findSpeed(body38)
    findSpeed(body39)
    findSpeed(body40)

    # write current data to file; if at a time step of 20*dt
    if ii%20 == 0 :
        dataFile.write(str(t))
        dataFile.write(' ')
        writeBodyToFile(Jupiter)
        dataFile.write(' ')
        writeBodyToFile(body1)
        dataFile.write(' ')
        writeBodyToFile(body2)
        dataFile.write(' ')
        writeBodyToFile(body3)
        dataFile.write(' ')
        writeBodyToFile(body4)
        dataFile.write(' ')
        writeBodyToFile(body5)
        dataFile.write(' ')
        writeBodyToFile(body6)
        dataFile.write(' ')
        writeBodyToFile(body7)
        dataFile.write(' ')
        writeBodyToFile(body8)
        dataFile.write(' ')
        writeBodyToFile(body9)
        dataFile.write(' ')
        writeBodyToFile(body10)
        dataFile.write(' ')
                   
        writeBodyToFile(body11)
        dataFile.write(' ')
        writeBodyToFile(body12)
        dataFile.write(' ')
        writeBodyToFile(body13)
        dataFile.write(' ')
        writeBodyToFile(body14)
        dataFile.write(' ')
        writeBodyToFile(body15)
        dataFile.write(' ')
        writeBodyToFile(body16)
        dataFile.write(' ')
        writeBodyToFile(body17)
        dataFile.write(' ')
        writeBodyToFile(body18)
        dataFile.write(' ')
        writeBodyToFile(body19)
        dataFile.write(' ')
        writeBodyToFile(body20)
        dataFile.write(' ')
                   
        writeBodyToFile(body21)
        dataFile.write(' ')
        writeBodyToFile(body22)
        dataFile.write(' ')
        writeBodyToFile(body23)
        dataFile.write(' ')
        writeBodyToFile(body24)
        dataFile.write(' ')
        writeBodyToFile(body25)
        dataFile.write(' ')
        writeBodyToFile(body26)
        dataFile.write(' ')
        writeBodyToFile(body27)
        dataFile.write(' ')
        writeBodyToFile(body28)
        dataFile.write(' ')
        writeBodyToFile(body29)
        dataFile.write(' ')
        writeBodyToFile(body30)
        
        writeBodyToFile(body31)
        dataFile.write(' ')
        writeBodyToFile(body32)
        dataFile.write(' ')
        writeBodyToFile(body33)
        dataFile.write(' ')
        writeBodyToFile(body34)
        dataFile.write(' ')
        writeBodyToFile(body35)
        dataFile.write(' ')
        writeBodyToFile(body36)
        dataFile.write(' ')
        writeBodyToFile(body37)
        dataFile.write(' ')
        writeBodyToFile(body38)
        dataFile.write(' ')
        writeBodyToFile(body39)
        dataFile.write(' ')
        writeBodyToFile(body40)
        dataFile.write('\n')

    #calculate new values for position and velocity
    verletJ(Jupiter)
    verlet(body1)
    verlet(body2)
    verlet(body3)
    verlet(body4)
    verlet(body5)
    verlet(body6)
    verlet(body7)
    verlet(body8)
    verlet(body9)
    verlet(body10)
                   
    verlet(body11)
    verlet(body12)
    verlet(body13)
    verlet(body14)
    verlet(body15)
    verlet(body16)
    verlet(body17)
    verlet(body18)
    verlet(body19)
    verlet(body20)
                   
    verlet(body21)
    verlet(body22)
    verlet(body23)
    verlet(body24)
    verlet(body25)
    verlet(body26)
    verlet(body27)
    verlet(body28)
    verlet(body29)
    verlet(body30)
    
    verlet(body31)
    verlet(body32)
    verlet(body33)
    verlet(body34)
    verlet(body35)
    verlet(body36)
    verlet(body37)
    verlet(body38)
    verlet(body39)
    verlet(body40)
    
    #update values; now -> prev, next -> now
    updatebody(Jupiter)
    updatebody(body1)
    updatebody(body2)
    updatebody(body3)
    updatebody(body4)
    updatebody(body5)
    updatebody(body6)
    updatebody(body7)
    updatebody(body8)
    updatebody(body9)
    updatebody(body10)
                   
    updatebody(body11)
    updatebody(body12)
    updatebody(body13)
    updatebody(body14)
    updatebody(body15)
    updatebody(body16)
    updatebody(body17)
    updatebody(body18)
    updatebody(body19)
    updatebody(body20)
                   
    updatebody(body21)
    updatebody(body22)
    updatebody(body23)
    updatebody(body24)
    updatebody(body25)
    updatebody(body26)
    updatebody(body27)
    updatebody(body28)
    updatebody(body29)
    updatebody(body30)
    
    updatebody(body31)
    updatebody(body32)
    updatebody(body33)
    updatebody(body34)
    updatebody(body35)
    updatebody(body36)
    updatebody(body37)
    updatebody(body38)
    updatebody(body39)
    updatebody(body40)

    t = t + dt
    

dataFile.close()
