#Author: Francisco Chaves
#Creation Date: 06.01.2015
#Version 0.00

#Import Libraries
import math
from numpy import genfromtxt

#Gobal Variables
global g
g = 9.81 #m/s
global pi
pi = math.pi
global rhoW
rhoW= 1030 #kg/m3
global rhoR
rhoR = 2650 #kg/m3

#Classes

class WavePeriod(object):
	"""Calculation relevant wave period parameters based on input wave period. \nargs = (t,opt)\nt (s): Known Wave period\n opt: 'tp', 'tm', 'tm10' whichever is known"""
	def __init__(self,t,opt):
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
		self.u = Hs/2 * (g*tp/L) * math.cosh(2*pi*(z+d)/L) / (math.cosh(2*pi*d/L))
		self.w = Hs/2 * (g*tp/L) * math.sinh(2*pi*(z+d)/L) / (math.cosh(2*pi*d/L))
		self.ax = (g*pi*Hs/L) * math.cosh(2*pi*(z+d)/L) / (math.cosh(2*pi*d/L))
		self.az = (g*pi*Hs/L) * math.sinh(2*pi*(z+d)/L) / (math.cosh(2*pi*d/L))

class WaveHeight(object):
	"""Calculates wave height for a given depth.
	Hs (m) : significant wave height
	Tp (m) : peak wave period
	d (m) : water depth
	slope : cot (alpha)  (1 in x)
	"""
	batjjes = genfromtxt('batjjes.cv',delimiter=',') #import table with normalized wave heights from batjjes&groenendijk 2000, Wave height distribution on shallow foreshores
	def __init__(self,Hs,Tp,d,slope):
			self.Htr = (0.35+5.8*1/slope)*d
			
			if Hs/d > 0.78:
				self.H = 0.78*d
			else:
				self.H = Hs

class Chezy(object):
	"""Calculation of Chezy coefficient based on the particle size of the rock protection
	h (m): water depth
	ks (m) : particle diameter
	"""
	def __init__(self,h,ks):
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
		self.tawC=rhoW*g*pow(U,2)/(pow(C,2))

class ShearStressWaves(object):
	"""Calculation of bed shear stress induced by waves
	u0 (m) : peak orbital velocity
	T (s) : peak wave period
	ks (m) : particle diameter
	"""
	def __init__(self,u0,T,ks):
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



