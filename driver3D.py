import numpy as np
import math
import shapefunctions
import quadrature
# User defined inputs
h = 1/1
ngpt = 3
order = 2
dof = 3
dx = 0#1/2
dy = 0#1/2
dz = 0#1/2
E = 1
nu = 0.25
lam = E*nu/((1-2*nu)*(1+nu))
mu = E/(2*(1+nu)) 

# Allocate
order2 = order*order
order3 = order*order*order
harr = np.array([h,h,h])
dxarr = np.array([dx,dy,dz])
elecoord = np.zeros((order3,3))
ktotal = np.zeros((order3*dof,order3*dof))
kmat   = np.zeros((order3*dof,order3*dof))
x      = np.zeros((3,1))
eye    = np.identity(dof)
Tmat   = np.array([[1., 0., 0., 0., 0., 0., 0., 0., 0.],
                   [0., 0., 0., 1., 0., 0., 0., 0., 0.],
                   [0., 0., 0., 0., 0., 0., 1., 0., 0.],
                   [0., 1., 0., 0., 0., 0., 0., 0., 0.],
                   [0., 0., 0., 0., 1., 0., 0., 0., 0.],
                   [0., 0., 0., 0., 0., 0., 0., 1., 0.],
                   [0., 0., 1., 0., 0., 0., 0., 0., 0.],
                   [0., 0., 0., 0., 0., 1., 0., 0., 0.],
                   [0., 0., 0., 0., 0., 0., 0., 0., 1.]])

# Elements in a structured grid
for i in range(0,order):
  for j in range(0,order):
    for k in range(0,order):
      elecoord[i*order2+j*order+k][0] = dxarr[0] + k*harr[0]
      elecoord[i*order2+j*order+k][1] = dxarr[1] + j*harr[1]
      elecoord[i*order2+j*order+k][2] = dxarr[2] + i*harr[2]

# Get quadrature
wlist,zlist = quadrature.gauss3D(ngpt)

# Compute element stiffness matrix
for i in range(0,len(wlist)):
  zpair = zlist[i]
  N,DN = shapefunctions.shape3D(zpair[0],zpair[1],zpair[2],order)
  J = np.matmul(elecoord.T,DN)
  x = np.matmul(elecoord.T,N.T)
  invJ = np.linalg.inv(J)
  detJ = np.linalg.det(J)
  B = np.matmul(DN,invJ)
  vecB = np.reshape(B,(order3*dof,1))
  kmat = np.matmul(np.kron(B,eye),np.kron(B.T,eye))*mu + np.matmul(np.matmul(np.kron(B,eye),Tmat),np.kron(B.T,eye))*mu
  kmat = kmat + np.matmul(vecB,vecB.T)*lam
  ktotal = ktotal + wlist[i]*kmat*detJ
print(ktotal)