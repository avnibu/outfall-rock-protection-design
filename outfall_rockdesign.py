#Author: Francisco Chaves
#Creation Date: 06.01.2015
#Version 0.00

#Import Libraries
import math
from math import pi
from numpy import genfromtxt
import numpy
numpy.set_printoptions(suppress=True) #supress scientific notation for numbers.
import openpyxl

#Gobal Variables
global g
g = 9.81 #m/s
global rhoW
rhoW= 1030 #kg/m3
global rhoR
rhoR = 2650 #kg/m3

#Classes

class WavePeriod(object):
	"""Calculation relevant wave period parameters based on input wave period. \nargs = (t,opt)\nt (s): Known Wave period\n opt: 'tp', 'tm', 'tm10' whichever is known"""
	def __init__(self,t,opt):
		t = float(t)
		if opt == 'tp':
			self.tp = t
			self.tm = t/1.2
			self.tm10 = t/1.1
		elif opt == 'tm':
			self.tp = t*1.2
			self.tm = t
			self.tm10 = t*1.2/1.1
		elif opt == 'tm10':
			self.tp = t*1.1
			self.tm = t*1.1/1.2
			self.tm10 = t

class WaveLength(object):
	"""Calculation of offshore and nzearshore wavelength.
	tp (s) : peak wave period
	d (m) : local water depth
	"""
	def __init__(self,tp,d):
		self.L0 = g*pow(tp,2)/(2*pi)
		while True:
			self.L = self.L0 #inititate self.L
			self.LIteration = 0 #inititate self.LIteration
			while abs(self.L-self.LIteration)>0.0001:
				self.LIteration = self.L
				self.L = self.L0 * math.tanh(2*math.pi*d/self.LIteration)
			break

class WaveMotion(object):
	"""Calculation of vertical and horizontal wave velocity along a water column(Linear Wave Theory). Here, cos(theta) and sin (theta) are always 1.0, so the velocity and acceleration calculated is maximum
	Hs (m) : significant wave height
	tp (s) : peak wave period
	L (m) : wave length
	z (m) : maximum water depth (still water level - seabed level)
	d (m) : water depth along the water column at which the wave velocity is going to be calculated
	"""
	def __init__(self, Hs,tp,L,z,d):
		Hs = float(Hs)
		tp = float(tp)
		L = float(L)
		z = float(z)
		d = float(d)
		self.u = Hs/2 * (g*tp/L) * math.cosh(2*pi*(z+d)/L) / (math.cosh(2*pi*d/L))
		self.w = Hs/2 * (g*tp/L) * math.sinh(2*pi*(z+d)/L) / (math.cosh(2*pi*d/L))
		self.ax = (g*pi*Hs/L) * math.cosh(2*pi*(z+d)/L) / (math.cosh(2*pi*d/L))
		self.az = (g*pi*Hs/L) * math.sinh(2*pi*(z+d)/L) / (math.cosh(2*pi*d/L))

class WaveHeight(object):
	"""Calculates wave height for a given depth.
	Hs (m) : significant wave height
	d (m) : water depth
	slope : cot (alpha)  (1 in x)
	"""
	def __init__(self,Hs,d,slope):
		Hs = float(Hs)
		d = float(d)
		slope = float(slope)
		battjes = genfromtxt("battjes.csv",delimiter=',') #import table with normalized wave heights from Battjes&Groenendijk 2000, Wave height distribution on shallow foreshores
		if Hs/d >= 0.78:
			self.Hs = 0.78*d
		else:
			self.Htr = (0.35+5.8*1/slope)*d
			# Hrms equation .59 The Rock Manual (page 359)
			self.Hrms = (0.6725 + 0.2025*(Hs/d))*Hs
			# calculate the normalised Htr
			HtrNorm = self.Htr / self.Hrms
			#find nearest to self.Htr in column 1 of Battjes. Choose the value immediately next to it.
			index = int(HtrNorm / 0.05) + 1
			if index > 60:
				index = 60
			#extract the relevant wave heights from Battjes table.
			self.Hs = battjes[index,3] * self.Hrms
			self.H2Percent = battjes[index,5] * self.Hrms
			self.H1Percent = battjes[index,6] * self.Hrms
			self.Hmax = battjes[index,7] * self.Hrms

class Chezy(object):
	"""Calculation of Chezy coefficient based on the particle size of the rock protection
	h (m): water depth
	ks (m) : particle diameter
	"""
	def __init__(self,h,ks):
		h = float(h)
		ks = float(ks)
		if h / ks > 2:
			self.C = 18*math.log(1+12*h/ks)
		else:
			self.C = 18*math.log(12*h/ks)

class ShearStressCurrent(object):
	"""Calculation of bed shear stress induced by currents
	U (m/s): depth averaged current speed
	C (m^(1/2)/s) : Chezy coefficient
	"""
	def __init__(self,U,C):
		U = float(U)
		C = float(C)
		self.tawC=rhoW*g*pow(U,2)/(pow(C,2))

class ShearStressWaves(object):
	"""Calculation of bed shear stress induced by waves
	u0 (m) : peak orbital velocity
	T (s) : peak wave period
	ks (m) : particle diameter
	"""
	def __init__(self,u0,T,ks):
		u0 = float(u0)
		T = float(T)
		ks = float(ks)
		self.a0 = u0*T/(2*pi)
		if self.a0 > 0.636*ks:
			self.fw = 0.237*pow((self.a0/ks),-0.52)
		else:
			self.fw = 0.3
		self.tawW = 0.5 * rhoW * self.fw * pow(u0,2)

class ShearStress(object):
	"""Combines wave and current induced shear stress
	tawC (N/m2) : current induced shear stress
	tawW (N/m2) : wave induced shear stress
	"""
	def __init__(self,tawC,tawW):
		if tawC > 0.4 * tawW:
			self.tawCW = tawC + 0.5*tawW
		else:
			self.tawCW = tawC + tawW

class RockDesign(object):
	"""Calculation the required rock size for the rock protection
	tawCW (N/m2) : combined shear stress induced by waves and currents
	"""
	def __init__(self,tawCW):
		self.D = 0.20 #initialize rock diameter as 0.2m
		Diteration = 0.1 #initialize iteration variable
		while (self.D-Diteration)>0.00001:
			Diteration = self.D
			self.shields = tawCW/((rhoR-rhoW)*g*self.D)

# 1. Grabs input data from an excel file

inputFile = openpyxl.load_workbook('inputFile.xlsx').get_sheet_by_name('Sheet1')
inputHs = float(inputFile['B5'].value)
inputPeriodType = str(inputFile['B6'].value)
if inputPeriodType == "Tm-1,0":
	inputPeriodType = "tm10"
inputT = float(inputFile['B7'].value)
inputSWL = float(inputFile['B8'].value)
inputSlope = float(inputFile['B9'].value)
inputCurrentSpeed = float(inputFile['B10'].value)

# 2. Bathymetry - Extract bathymetric data from .csv file. The .csv file should be formatted as COL1:X COL2:Y COL3:Z

bathy = genfromtxt("bathy.csv",delimiter=',') #import data from bathy.csv using numpy genfromtxt tool

#generate fourth vector with distance between points

bathyDistance = numpy.empty(len(bathy)) #generate an empty array with the same number of elements in a column of bathy
bathyDistance = bathyDistance.reshape(len(bathyDistance),1) #transverse the previous array (could have been done in one line only. (Two lines used for better understanding)
bathyDepth = numpy.empty(len(bathy)).reshape(len(bathy),1) #generate an empty array with the same number of elements in a column of bathy already transversed

for i in range(2,len(bathy)): #calculate the distance between each point (which will be used for plotting later)
	bathyDistance[i] = ((bathy[i,0]-bathy[i-1,0])**2+(bathy[i,1]-bathy[i-1,1])**2)**0.5
for i in range(1,len(bathy)):
	bathyDepth[i]= inputSWL - bathy[i,2]
bathy = numpy.hstack([bathy,bathyDistance]) #appends bathyDistance as the last column of bathy.
bathy = numpy.hstack([bathy,bathyDepth]) #appends bathyDepth as the last column of bathy.

# 3. Calculation - Generates the results and compiles them in an array
#results array description CUMULATIVE_DISTANCE/WATER_DEPTH/WAVE_HEIGT

results = numpy.empty(len(bathy)).reshape(len(bathy),1) #generates an empty array with the same size as bathy

# 3.1 CUMULATIVE_DISTANCE - substitute the first column of results with distance between points
results[0] = 0 #set the value to zero (because it is text)
results[1] = 0 #set the value to zero (because it is the first point)
for i in range(2,len(bathy)):
	results[i] = results[i-1] + bathy[i,3]

# 3.2 WATER_DEPTH - append water depth extracted from bathy
results = numpy.hstack([results,bathy[:,4].reshape(len(bathy),1)])

# 3.3 WAVE_HEIGHT - calculate Wave Heights using Battjes
# first we create an empty vector with the same length as results.

resultsWaves = numpy.empty(len(results)).reshape(len(results),1)

for i in range(0,len(results)):
	resultsWaves[i] = WaveHeight(inputHs,results[i,1],inputSlope).Hs

#then we append resultsWaves to the results array

results = numpy.hstack([results,resultsWaves])

# 3.4 WAVE PERIOD calculation

## calculate wave periods (Tp, Tm, Tm-1,0) based on the wave period input

tp = WavePeriod(inputT,inputPeriodType.lower()).tp
tm = WavePeriod(inputT,inputPeriodType.lower()).tm
tm10 = WavePeriod(inputT,inputPeriodType.lower()).tm10

# 3.5 WAVE LENGTH calculation - Wave length varies with water depth so we need to calculate it for each point in results
## first we need to create an empty array with the same length as results

resultsWaveLength = numpy.empty(len(results)).reshape(len(results),1)

## then we calculate the wave length for each point

for i in range(1,len(results)):
	resultsWaveLength[i] = WaveLength(tp,results[i,1])
	print resultsWaveLength[i]

## and then we append it to the results array

#results = numpy.hstack([results,resultsWaveLength])

# 3.5 WAVE_VELOCITY_NEARBED
# first we create an empty array with the same length as results

#resultsWaveVelocityNearbed = numpy.empty(len(results)).reshape(len(results),1)

# calculate wave velocity nearbed

#resultsWaveVelocityNearbed[i] = WaveMotion(results[i,2],tp,)