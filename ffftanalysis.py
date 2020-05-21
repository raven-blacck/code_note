import numpy as np

def karray(grid, BE=False):
	'''
	Return the mod of k in the fft order
	If BE is true, return kx,ky,kz
	'''
	kxx=np.zeros((grid,1,1))
	kyy=np.zeros((1,grid,1))
	kzz=np.zeros((1,1,grid//2+1))
	kzz[0,0,:]=np.arange(grid//2+1)
	for ii in [0,1]:
		kxx[ii*grid//2:(ii+1)*grid//2,0,0]=np.arange(grid//2)-ii*grid//2
		kyy[0,ii*grid//2:(ii+1)*grid//2,0]=np.arange(grid//2)-ii*grid//2
	kmod=np.sqrt(kxx**2+kyy**2+kzz**2)
	kx=kxx+0*kyy+0*kzz
	ky=0*kxx+kyy+0*kzz
	kz=0*kxx+0*kyy+kzz
	
	if BE: return kmod, kx,ky,kz
	return kmod

