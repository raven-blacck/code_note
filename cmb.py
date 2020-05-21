import numpy as np
import matplotlib.pyplot as plt
import camb

T_cmb=2.728*10**6 #uk
def cmb_powerspectrum(T_cmb=2.728*10**6, H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.06, omk=0, tau=0.06,As=2e-9, ns=0.965, r=0):
	pars = camb.CAMBparams()
	pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2, mnu=mnu, omk=omk, tau=tau)
	pars.InitPower.set_params(As=As, ns=ns, r=r)
	pars.set_for_lmax(10000, lens_potential_accuracy=0)
	results = camb.get_results(pars)
	powers =results.get_cmb_power_spectra(pars,CMB_unit="muK",raw_cl=True)
	l=np.arange(len(powers["total"][:,0]))
	coef=l*(l+1)/2/np.pi*T_cmb**2;coef[0]=1
	coef=T_cmb**2
	np.save("lensedcmbcl",powers["total"][:,0]/coef)
	np.save("unlensedcmbcl",powers["unlensed_scalar"][:,0]/coef)
	
def cl2map2d(cl,ang,grid):
	# unit of ang is rad, is the length of the square

	#get a array of cl with percise of 0.01
	def findcl(l):
		if l>=len(cl): return 0
		w2=l-int(l)
		w1=int(l)+1-l
		return w1*cl[int(l)]+w2*cl[int(l)]

	l_start=int(2*np.pi/ang)-5
	l_end=int(2*grid*np.pi/ang)+5
	l_array=np.arange(l_start,l_end,0.01)
	cl_array=np.zeros(len(l_array))
	for i in range(len(l_array)):
		cl_array[i]=findcl(l_array[i])

	#gen a noise map, pk=1 white noise
	map_noise=np.random.rand(grid,grid)
	coef1=2*np.pi*map_noise[0:grid:2,:]
	coef2=np.sqrt(-2*np.log(map_noise[1:grid:2,:]))
	map_noise[0:grid:2,:]=coef2*np.cos(coef1)
	map_noise[1:grid:2,:]=coef2*np.sin(coef1)

	map_noisek =np.fft.rfft2(map_noise)
	#coef of map_noisek from cl
	coefmap=np.zeros((grid,grid//2+1))
	for i in range(grid):
		if i>grid/2: ii=i-grid
		else: ii=i
		for j in range(grid//2+1):
			jj=j
			llen=np.sqrt(ii**2+jj**2)*2*np.pi/ang
			if llen==0:continue
			coefmap[i,j]=cl_array[int((llen-l_start)//0.01)]
	map_noisek*=np.sqrt(coefmap*(grid/ang)**2)
	map2d=np.fft.irfft2(map_noisek)
	return map2d

def cmb2d(ang,grid,T_cmb=2.728*10**6,Delta_T=2,sigma=7):
	cl=np.load("lensedcmbcl.npy")
	cmbmap=cl2map2d(cl,ang,grid)
	#add_noise	
	arcmin2rad=1.0/60/180*np.pi
	l=np.arange(10000)
	cln=(Delta_T*arcmin2rad/T_cmb)**2*np.exp(l*(l+1)*(sigma*arcmin2rad)**2/8/np.log(2))
	noisemap=cl2map2d(cln,ang,grid)
	return cmbmap+noisemap

	