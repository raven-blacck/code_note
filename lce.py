import numpy as np
import fftanalysis as fa
from astropy.cosmology import WMAP9 as cosmo


def postoden(x,grid):
	print("Position of particle to density field",end="\t")
	denhalo=np.histogramdd(x,bins=(grid,grid,grid))[0]
	norm=grid**3*1.0/len(x[:,0])
	denhalo=denhalo*norm
	print("end")
	return denhalo

def linear_continuous_eq_den(denhalo,z,Omega_m,boxlen):
	def lce_norm():
		E=Omega_m*(1+z)**3+(1-Omega_m)
		f=(Omega_m*(1+z)**3/E)**0.6
		H=cosmo.H(z).value
		k_unit=2*np.pi/boxlen
		return f*H/k_unit/(1+z)/300000
	print("Linear continuous equation, using den to get vel")
	t0=time.time()
	kmod,kx,ky,kz=fa.karray(grid,True)
	denk=np.fft.rfftn(denhalo)
	velxk=np.zeros(denk.shape,dtype=type(denk[0,0,0]))
	velyk=np.zeros(denk.shape,dtype=type(denk[0,0,0]))
	velzk=np.zeros(denk.shape,dtype=type(denk[0,0,0]))
	kmod[0,0,0]=1
	denk[0,0,0]=0+0j
	velxk+=kx*denk/kmod**2*(1j)
	velyk+=ky*denk/kmod**2*(1j)
	velzk+=kz*denk/kmod**2*(1j)
	momhalo=np.zeros((grid,grid,grid,3))
	momhalo2=np.zeros((grid,grid,grid,3))
	velhalo=np.zeros((grid,grid,grid,3))
	velhalo[:,:,:,0]=np.fft.irfftn(velxk)
	momhalo[:,:,:,0]=(denhalo-1)*velhalo[:,:,:,0]
	velhalo[:,:,:,1]=np.fft.irfftn(velyk)
	momhalo[:,:,:,1]=(denhalo-1)*velhalo[:,:,:,1]
	velhalo[:,:,:,2]=np.fft.irfftn(velzk)
	momhalo[:,:,:,2]=(denhalo-1)*velhalo[:,:,:,2]
	momhalo2[:,:,:,0]=(denhalo)*velhalo[:,:,:,0]
	momhalo2[:,:,:,1]=(denhalo)*velhalo[:,:,:,1]
	momhalo2[:,:,:,2]=(denhalo)*velhalo[:,:,:,2]
	
	print ("Momhalo over ----- ", np.round(time.time()-t0,3))
	return denhalo,velhalo*lce_norm(),momhalo*lce_norm(),momhalo2*lce_norm()	#vel is scale by light speed

def linear_continuous_eq_pos(x,grid,z,Omega_m,boxlen):
	den=postoden(x,grid)
	den,vel,mom,mom2=linear_continuous_eq_den(denhalo,z,Omega_m,boxlen)
	return den,vel,mom,mom2


