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

class WaveLength(WavePeriod):
	"""Calculation of offshore and nearshore wavelength. \nargs = (TP,D)\nTP (s) = peak wave period\n d (m) = water depth"""
	def __init__(self,t,opt,tp,d):
		super(WavePeriod,self).__init__(t,opt)
		self.L0 = g*pow(tp,2)/(2*pi)
		while True:
			self.L = self.L0 #inititate self.L
			self.LIteration = 0 #inititate self.LIteration
			while abs(self.L-self.LIteration)>0.0001:
				self.LIteration = self.L
				self.L = self.L0 * math.tanh(2*math.pi*d/self.LIteration)
			break

class WaveVelocity:
	"""Calculation of vertical and horizontal wave velocity along a water column(Linear Wave Theory)"""
	def __init__(self, Hs,tp,L):
		pass
