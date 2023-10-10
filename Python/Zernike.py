import numpy as np

def zernike_polar(coefficient,r,u):
	"""
	Zernike Polynomials Caculation in polar coordinates
	inputs: coefficient: Zernike Polynomials Coefficient from input
				r: rho in polar coordinates
				u: theta in polar coordinates
	returns total value
	"""
	Z = [0]+coefficient
	Z1  =  Z[1]  * 1*(np.cos(u)**2+np.sin(u)**2)
	Z2  =  Z[2]  * 2*r*np.cos(u)
	Z3  =  Z[3]  * 2*r*np.sin(u)
	Z4  =  Z[4]  * np.sqrt(3)*(2*r**2-1)
	Z5  =  Z[5]  * np.sqrt(6)*r**2*np.sin(2*u)
	Z6  =  Z[6]  * np.sqrt(6)*r**2*np.cos(2*u)
	Z7  =  Z[7]  * np.sqrt(8)*(3*r**2-2)*r*np.sin(u)
	Z8  =  Z[8]  * np.sqrt(8)*(3*r**2-2)*r*np.cos(u)
	Z9  =  Z[9]  * np.sqrt(8)*r**3*np.sin(3*u)
	Z10 =  Z[10] * np.sqrt(8)*r**3*np.cos(3*u)
	Z11 =  Z[11] * np.sqrt(5)*(1-6*r**2+6*r**4)
	Z12 =  Z[12] * np.sqrt(10)*(4*r**2-3)*r**2*np.cos(2*u)
	Z13 =  Z[13] * np.sqrt(10)*(4*r**2-3)*r**2*np.sin(2*u)
	Z14 =  Z[14] * np.sqrt(10)*r**4*np.cos(4*u)
	Z15 =  Z[15] * np.sqrt(10)*r**4*np.sin(4*u)
	Z16 =  Z[16] * np.sqrt(12)*(10*r**4-12*r**2+3)*r*np.cos(u)
	Z17 =  Z[17] * np.sqrt(12)*(10*r**4-12*r**2+3)*r*np.sin(u)
	Z18 =  Z[18] * np.sqrt(12)*(5*r**2-4)*r**3*np.cos(3*u)
	Z19 =  Z[19] * np.sqrt(12)*(5*r**2-4)*r**3*np.sin(3*u)
	Z20 =  Z[20] * np.sqrt(12)*r**5*np.cos(5*u)
	Z21 =  Z[21] * np.sqrt(12)*r**5*np.sin(5*u)
	Z22 =  Z[22] * np.sqrt(7)*(20*r**6-30*r**4+12*r**2-1)
	Z23 =  Z[23] * np.sqrt(14)*(15*r**4-20*r**2+6)*r**2*np.sin(2*u)
	Z24 =  Z[24] * np.sqrt(14)*(15*r**4-20*r**2+6)*r**2*np.cos(2*u)
	Z25 =  Z[25] * np.sqrt(14)*(6*r**2-5)*r**4*np.sin(4*u)
	Z26 =  Z[26] * np.sqrt(14)*(6*r**2-5)*r**4*np.cos(4*u)
	Z27 =  Z[27] * np.sqrt(14)*r**6*np.sin(6*u)
	Z28 =  Z[28] * np.sqrt(14)*r**6*np.cos(6*u)
	Z29 =  Z[29] * 4*(35*r**6-60*r**4+30*r**2-4)*r*np.sin(u)
	Z30 =  Z[30] * 4*(35*r**6-60*r**4+30*r**2-4)*r*np.cos(u)
	Z31 =  Z[31] * 4*(21*r**4-30*r**2+10)*r**3*np.sin(3*u)
	Z32 =  Z[32] * 4*(21*r**4-30*r**2+10)*r**3*np.cos(3*u)
	Z33 =  Z[33] * 4*(7*r**2-6)*r**5*np.sin(5*u)
	Z34 =  Z[34] * 4*(7*r**2-6)*r**5*np.cos(5*u)
	Z35 =  Z[35] * 4*r**7*np.sin(7*u)
	Z36 =  Z[36] * 4*r**7*np.cos(7*u)
	Z37 =  Z[37] * 3*(70*r**8-140*r**6+90*r**4-20*r**2+1)


	Z = Z1 + Z2 +  Z3+  Z4+  Z5+  Z6+  Z7+  Z8+  Z9+ \
		Z10+ Z11+ Z12+ Z13+ Z14+ Z15+ Z16+ Z17+ Z18+ Z19+ \
		Z20+ Z21+ Z22+ Z23+ Z24+ Z25+ Z26+ Z27+ Z28+ Z29+ \
		Z30+ Z31+ Z32+ Z33+ Z34+ Z35+ Z36+ Z37
	return Z

def find_zernike(RAD,THET,VAL, n):
	"""
	find first n Zernike coefficients
	inputs: RAD=radial coordinates, 1d array
				THET=angular coordinates 1d array
				VAL=measured values 1d array
	return:	list of coefficients
	"""
	fitlist = []
	ll = len(RAD)
	for i in range(n):
		c = [0]*i+[1]+[0]*(37-i-1) # pick out desired term
		ZF = zernike_polar(c,RAD,THET)
		a = 1.3*sum(VAL*ZF)/ll
		fitlist.append(a)
	l1 = len(fitlist)
	fitlist = fitlist+[0]*(37-l1)
	return fitlist

def find_zernike_2d(RAD,THET,VAL, n):
	"""
	find first n Zernike coefficients
	inputs: RAD=radial coordinates, 1d array
				THET=angular coordinates 1d array
				VAL=measured values 1d array
	return:	list of coefficients
	"""
	fitlist = []
	ll = len(RAD)
	for i in range(n):
		c = [0]*i+[1]+[0]*(37-i-1) # pick out desired term
		ZF = zernike_polar(c,RAD,THET)
		a = 1.3*sum(sum(VAL*ZF))/ll/ll
		fitlist.append(a)
	l1 = len(fitlist)
	fitlist = fitlist+[0]*(37-l1)
	return fitlist
