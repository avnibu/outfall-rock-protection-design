#Author: Francisco Chaves
#Creation Date: 06.01.2015
#Version 0.00

#Import Libraries
import math

#Gobal Variables
global g
g = 9.81 #m/s
global pi
pi = math.pi
global sinh
sinh = math.sinh
global cosh
cosh = math.cosh

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

