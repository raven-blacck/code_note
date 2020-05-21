import numpy as np

def proj(f,grid,dirc,the):
	def proj1(f,the=np.pi/6):
		
		s=np.sin(the)
		c=np.cos(the)
		x=np.arange(grid).reshape(-1,1,1)
		y=np.arange(grid).reshape(1,-1,1)
		z=np.arange(grid).reshape(1,1,-1)

		rx=(c*x-s*y+0*z).astype(np.int)%grid
		ry=(s*x+c*y+0*z).astype(np.int)%grid
		rz=(0*x+0*y+1*z).astype(np.int)%grid

		fp=(c*f[rx,ry,rz,0]-s*f[rx,ry,rz,1])

		return np.sum(fp,axis=0)

	def proj2(f,the=np.pi/6):
		
		s=np.sin(the)
		c=np.cos(the)
		x=np.arange(grid).reshape(-1,1,1)
		y=np.arange(grid).reshape(1,-1,1)
		z=np.arange(grid).reshape(1,1,-1)

		rx=(1*x+0*y+0*z).astype(np.int)%grid
		ry=(0*x+c*y-s*z).astype(np.int)%grid
		rz=(0*x+s*y+c*z).astype(np.int)%grid

		fp=(c*f[rx,ry,rz,1]-s*f[rx,ry,rz,2])

		return np.sum(fp,axis=1)


	def proj3(f,the=np.pi/6):
		
		s=np.sin(the)
		c=np.cos(the)
		x=np.arange(grid).reshape(-1,1,1)
		y=np.arange(grid).reshape(1,-1,1)
		z=np.arange(grid).reshape(1,1,-1)

		rx=(c*x+0*y-s*z).astype(np.int)%grid
		ry=(0*x+1*y+0*z).astype(np.int)%grid
		rz=(s*x+0*y+c*z).astype(np.int)%grid

		fp=(s*f[rx,ry,rz,0]+c*f[rx,ry,rz,2])

		return np.sum(fp,axis=2)

	func=[proj1,proj2,proj3]
	result=func[dirc](f,the)
	return result
