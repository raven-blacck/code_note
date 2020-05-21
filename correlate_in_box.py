import numpy as np
import healpy as hp
from tqdm import tqdm
import sys 
import fftanalysis as fa
import time
import struct
from astropy.cosmology import WMAP9 as cosmo
zstep=[5000,4724,4338,4098,3982,3656,3356,3080,2826,2593,2448,2378,2181,2000,1942,1780,1631,1494,1328,1216]
vs=np.array([17400000.0, 15923808.3, 14020373.5, 12931188.4, 12430474.5, 11111965.1, 10011775.9, 9091734.5, 8318855.9, 7667665.8, 7288078.4, 7111356.2, 6634895.0, 6221216.6, 6092865.1, 5743409.9, 5431402.9, 5149821.4, 4811557.1, 4582649.9])/300000	#v/c
boxlen=1200.0
Omega_m=0.268
pm=2.77*10**11*Omega_m*boxlen**3/3072.**3	#particle mass
hmthup=10**16/pm
hmthdown=10**12/pm
datadir="/data/s4/chenzy/"
savedatadir="/data/s4/chenzy/halomom/"
grid=512
filename="/data/s1/simu/Jing6620/info/zlist.txt"
filedata=np.loadtxt(filename,skiprows=1)
step2z=[]
for s in zstep:
	step2z.append(filedata[np.where(filedata[:,0]==s),2])
step2z=np.array(step2z).flatten()

nsrange=[0，6,10,16]

def lce_norm(ns):
	E=Omega_m*(1+step2z[ns])**3+(1-Omega_m)
	f=(Omega_m*(1+step2z[ns])**3/E)**0.6
	H=cosmo.H(step2z[ns]).value
	k_unit=2*np.pi/boxlen
	return f*H/k_unit/(1+step2z[ns])



def ReadHalos(zstep):
	t0=time.time()
	name='/data/s1/simu/Jing6620/gcatdbb02.6620.'+str(zstep);print("Reading "+name)
	pos=open(name,'rb')
	pos.read(4)
	ng=struct.unpack('l',pos.read(8))[0]	#number of haloes
	print("The number of halo"+str(ng))
	nmin=pos.read(8);nmin=struct.unpack('l',nmin)[0]
	print('ng:',ng,'nmin:',nmin)
	pos.read(8)
	nrich=pos.read(8*ng);nrich=np.array(struct.unpack(str(ng)+'l',nrich))
	#print("nrich:",nrich)
	pos.read(8)
	halox=np.zeros((ng,3))
	haloxx=pos.read(4*ng);halox[:,0]=np.array(struct.unpack(str(ng)+'f',haloxx));pos.read(8)
	haloxy=pos.read(4*ng);halox[:,1]=np.array(struct.unpack(str(ng)+'f',haloxy));pos.read(8)
	haloxz=pos.read(4*ng);halox[:,2]=np.array(struct.unpack(str(ng)+'f',haloxz));pos.read(8)
	halov=np.zeros((ng,3))
	halovx=pos.read(4*ng);halov[:,0]=np.array(struct.unpack(str(ng)+'f',halovx));pos.read(8)
	halovy=pos.read(4*ng);halov[:,1]=np.array(struct.unpack(str(ng)+'f',halovy));pos.read(8)
	halovz=pos.read(4*ng);halov[:,2]=np.array(struct.unpack(str(ng)+'f',halovz));pos.read(8)
	#print(haloxx);print(haloxy);print(haloxz);
	print("Read data of "+str(zstep)+" costs "+str(np.round(time.time()-t0,3))+'s')
	pos.close()
	return halox,halov,nrich,ng
#halox 位置
#halov 速度
#nrich halo内dm粒子数
#ng halo总数


def halomom(x,ns):
	t0=time.time()
	denhalo=np.histogramdd(x,bins=(grid,grid,grid))[0]
	norm=grid**3*1.0/len(x[:,0])
	denhalo=denhalo*norm
	kmod,kx,ky,kz=fa.karray(grid,True)
	denk=np.fft.rfftn(denhalo)
	velxk=np.zeros(denk.shape,dtype=type(denk[0,0,0]))
	velyk=np.zeros(denk.shape,dtype=type(denk[0,0,0]))
	velzk=np.zeros(denk.shape,dtype=type(denk[0,0,0]))
	'''
	velxk[1:,1:,1:]=kx[1:,1:,1:]/kmod[1:,1:,1:]**2
	velxk[1:,1:,1:]*=denk[1:,1:,1:]*(1j)
	velyk[1:,1:,1:]=ky[1:,1:,1:]/kmod[1:,1:,1:]**2
	velyk[1:,1:,1:]*=denk[1:,1:,1:]*(1j)
	velzk[1:,1:,1:]=kz[1:,1:,1:]/kmod[1:,1:,1:]**2
	velzk[1:,1:,1:]*=denk[1:,1:,1:]*(1j)
	'''
	velxk[1:,1:,1:]+=kx[1:,1:,1:]*denk[1:,1:,1:]/kmod[1:,1:,1:]**2*(1j)
	velyk[1:,1:,1:]+=ky[1:,1:,1:]*denk[1:,1:,1:]/kmod[1:,1:,1:]**2*(1j)
	velzk[1:,1:,1:]+=kz[1:,1:,1:]*denk[1:,1:,1:]/kmod[1:,1:,1:]**2*(1j)
	momhalo=np.zeros((grid,grid,grid,3))
	momhalo2=np.zeros((grid,grid,grid,3))
	velhalo=np.zeros((grid,grid,grid,3))
	'''
	momhalo[:,:,:,0]=(denhalo-1)
	momhalo[:,:,:,0]*=np.fft.irfftn(velxk)
	momhalo[:,:,:,1]=(denhalo-1)
	momhalo[:,:,:,1]*=np.fft.irfftn(velyk)
	momhalo[:,:,:,2]=(denhalo-1)
	momhalo[:,:,:,2]*=np.fft.irfftn(velzk)
	'''
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
	return denhalo,velhalo*lce_norm(ns),momhalo*lce_norm(ns),momhalo2*lce_norm(ns)	#km/s

