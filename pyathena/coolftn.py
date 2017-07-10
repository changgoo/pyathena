import astropy.constants as c
import astropy.units as u
import cPickle as p
import numpy as np

class coolftn(object):
	def __init__(self,fname='coolftn.p'):
		cdf=p.load(open(fname,'rb'))
		self.cool=np.array(cdf['cool'])
		self.heat=np.array(cdf['heat'])
		self.temp=np.array(cdf['temp'])
		self.T1=np.array(cdf['t1'])

		self.Tmin=self.T1.min()
		self.Tmax=self.T1.max()
		self.dT=np.log10(self.T1[1]/self.T1[0])
		self.nT=len(self.T1)

	def get_Tidx(self,T):
		if type(T)==np.ndarray:
			Tidx=np.log10(T/self.Tmin)/self.dT
			Tidx[np.where(T<self.Tmin)]=0
			Tidx[np.where(T>=self.Tmax)]=self.nT-2
			return Tidx.astype(int)
		else:
			if T < self.Tmin: return 0
			if T >= self.Tmax: return self.nT-2
			Tidx=np.log10(T/self.Tmin)/self.dT
			return int(Tidx)

	def get_temp(self,T1):
		T1idx=self.get_Tidx(T1)
		Ti=self.temp[T1idx]
		Tip1=self.temp[T1idx+1]
		T1i=self.T1[T1idx]
		T1ip1=self.T1[T1idx+1]
		T=Ti+(Tip1-Ti)*(T1-T1i)/(T1ip1-T1i)

		return T

	def get_cool(self,T1):
		T1idx=self.get_Tidx(T1)
		Li=self.cool[T1idx]
		Lip1=self.cool[T1idx+1]
		T1i=self.T1[T1idx]
		T1ip1=self.T1[T1idx+1]
		L=Li+(Lip1-Li)*(T1-T1i)/(T1ip1-T1i)

		return L
	
	def get_heat(self,T1):
		T1idx=self.get_Tidx(T1)
		Gi=self.heat[T1idx]
		Gip1=self.heat[T1idx+1]
		T1i=self.T1[T1idx]
		T1ip1=self.T1[T1idx+1]
		G=Gi+(Gip1-Gi)*(T1-T1i)/(T1ip1-T1i)

		return G
