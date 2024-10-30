#Equations
#[1] - Surface evaporation of water (https://www.engineeringtoolbox.com/evaporation-water-surface-d_690.html)
#[2] - Saturation Pressure water vapour (https://www.engineeringtoolbox.com/water-vapor-saturation-pressure-air-d_689.html)
#[3] - Air - Humidity Ratio (https://www.engineeringtoolbox.com/humidity-ratio-air-d_686.html) 
#[4] - Orifice flow - (https://www.omnicalculator.com/physics/orifice-flow)

import math
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import fsolve

def celciusToKelvin(pTemperature):
    return pTemperature + 273.15 #K

def kelvinToCelcius(pTemperature):
    return pTemperature - 273.15 #°C

T0 = celciusToKelvin(23) #K (roomTemperature)

hr = 0.05 #m (reflectorHeight)
wr = 0.50 #m (reflectorWidth)
Kr = 4 * hr / wr**2

def reflectorParabolaFormulaDerivative(x):
    #Kr * x**2
    return 2 * Kr * x

def lengthOfLineCalc(f, lowerBound, UpperBound):
    def integrand(x):
        return math.sqrt(f(x)**2 + 1)
    
    return quad(integrand, lowerBound, UpperBound)

ri = 0.04 #m (innerRadius)
ro = 0.05 #m (outerRadius)
Lc = 0.30 #m (cellLength)

rOut = 0.02 #m (outletRadius)
CSAOut = math.pi * rOut ** 2 #m²

rIn = 0.02 #m (inletRadius)
CSAIn = math.pi * rIn ** 2 #m²

h0 = 0.03 #m (originalWaterLevel)
tankH = 0.2 #m (tankHeight)

ρw = 1000 #kg m⁻³ (denistyWater)
m0 = ρw * Lc * ((math.acos(1 - h0/ri) * ri**2) - (ri - h0)*math.sqrt(ri**2 - (ri - h0)**2)) #kg (originalMassOfWater)

ls = 2 * lengthOfLineCalc(reflectorParabolaFormulaDerivative, ro, wr/2)[0] #(lengthOfCurveExposed)
Lr = 0.25 #m (reflectorLength)
As = ls * Lr #m² (areaOfParabolaExposed)
Cs = 1370 #W m⁻² (solarConstant)
P = As * Cs #W (solarPower)


cw = 4200 #J Kg⁻¹ K⁻¹ (specificHeatCapacity) (of water)
dt = 10e-2 #delta time

#P dt = mc dT
def calculate_dT(pMass): #K change in temperature
    return P * dt / (pMass * cw)


#gs = Θ A (x s - x) / 3600 [1]
def calculate_gs(pHeight, pTemperature): #kg s⁻¹ rate of mass loss due to evaporation
    pArea = Lc * ((math.acos(1 - pHeight/ri) * ri**2) - (ri - pHeight)*math.sqrt(ri**2 - (ri - pHeight)**2))
    pVelocity = 3.5 #m s⁻¹
    #Θ = 25 + 19v
    pΘ = 25 + 19 * pVelocity

    #pws = e^(77.3450 + 0.0057 T - 7235 / T) / T^8.2 [2]
    A = 77.3450
    B = 0.0057
    C = 7235
    D = 8.2

    pPws = math.exp(A + B * pTemperature - C / pTemperature) / (pTemperature ** D)

    #xs = 0.62198 pws / (pa - pws)  [3] 
    K = 0.62198 

    pXs = K * pPws / (101325 - pPws)

    R = 0
    return abs(pΘ * pArea * pXs * (1 - R) / 3600)

#dm = -gs dt
def calculate_dm(pMassLoss, pFlowRateOut, pFlowRateIn):
    return (-pMassLoss -pFlowRateOut + pFlowRateIn) * dt

#Q = Cd * A * sqrt(2*g*H) [4]
def calculate_QOut(pValveOpen, pHeightWaterAbove = 0, pHoleCSA = 0):
    Cd = 0.8 #() water discharge constant
    g = 9.81 #m s⁻² #acceleration due to graivty
    
    if pValveOpen:
        PUMP = 0.04
        return Cd * pHoleCSA * math.sqrt(2 * g * pHeightWaterAbove) + PUMP
    else:
        return 0

def calculate_QIn(pValveOpen, pHeightWaterAbove = 0, pHoleCSA = 0, pump = False):
    Cd = 0.8 #() water discharge constant
    g = 9.81 #m s⁻² #acceleration due to graivty
    
    PUMP = 0.04
    if pValveOpen:
        if pump:
            return Cd * pHoleCSA * math.sqrt(2 * g * pHeightWaterAbove) + PUMP
        return Cd * pHoleCSA * math.sqrt(2 * g * pHeightWaterAbove)
    else:
        return 0
    

#Calculate the Height of water
def calculate_h(pMass):
    def volumeEquation(h):
        return ρw * Lc * ((math.acos(1 - h/ri) * ri**2) - (ri - h)*math.sqrt(ri**2 - (ri - h)**2)) - pMass
    initialGuess = ri/2
    h_solution, = fsolve(volumeEquation, initialGuess)
    return h_solution

def updateVariables(pTime, pTemperature, pMass, pMassGained, pOutletValveOpen, pInletValveOpen, pInletPump = False):
    pHeight = calculate_h(pMass)
    pTime += dt
    pTemperature += calculate_dT(pMass)

    pWaterEvaporationRate =  calculate_gs(pHeight, pTemperature)
    pWaterOutRate = calculate_QOut(pOutletValveOpen, pHeight, CSAOut)
    pWaterInRate = calculate_QIn(pInletValveOpen, tankH, CSAIn, pInletPump)
    dm = calculate_dm(pWaterEvaporationRate, pWaterOutRate, pWaterInRate)

    pTemperature = (T0 * pWaterInRate * dt  + pTemperature * (M-dm))/((M-dm) + pWaterInRate * dt)

    pMass += dm
    pMassGained += pWaterEvaporationRate * dt
    return (pTime, pTemperature, pMass, pMassGained)

def updateLists(pL1, pL2 ,pL3 ,pL4, pElement1, pElement2, pElement3, pElement4):
    return (pL1 + [pElement1], pL2 + [pElement2], pL3 + [pElement3], pL4 + [pElement4])
    
time = 0
T = T0
baseMass = m0/4
drainMass = m0/10
M = baseMass
MG = 0

timeList = [time]
temperatureList = [kelvinToCelcius(T)]
massList = [M]
massGainedList = [MG]

for _ in range(4):
    print(T)
    
    while kelvinToCelcius(T) < 99:
        time, T, M, MG = updateVariables(time, T, M, MG, False, (M < m0))
        timeList, temperatureList, massList, massGainedList = updateLists(timeList, temperatureList, massList, massGainedList, time, kelvinToCelcius(T), M, MG)
    while M > drainMass:
        time, T, M, MG = updateVariables(time, T, M, MG, True, False)
        timeList, temperatureList, massList, massGainedList = updateLists(timeList, temperatureList, massList, massGainedList, time, kelvinToCelcius(T), M, MG)
    while M < baseMass:
        time, T, M, MG = updateVariables(time, T, M, MG, False, True, True)
        timeList, temperatureList, massList, massGainedList = updateLists(timeList, temperatureList, massList, massGainedList, time, kelvinToCelcius(T), M, MG)


fig, ax1 = plt.subplots()

ax1.plot(timeList, temperatureList, 'g-')
ax1.set_xlabel('time (s)')

ax2 = ax1.twinx()
ax2.plot(timeList, massList, 'b-')

fig.tight_layout()
fig.show()

fig2 = plt.figure(2)
plt.plot(timeList, massGainedList, color = "red")
fig2.show()
input()
