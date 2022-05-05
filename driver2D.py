import numpy as np
import math
import shapefunctions
import quadrature
# User defined inputs
h = 1/1
ngpt = 5
order = 2
dof = 2
dx = 0#1/2
dy = 0#1/2

E = 1.0
nu = 0.25
lam = E*nu/((1-2*nu)*(1+nu))
mu = E/(2*(1+nu)) 

# Allocate
order2 = order*order

harr = np.array([h,h])
dxarr = np.array([dx,dy])
elecoord = np.zeros((order2,2))
ktotal = np.zeros((order2*dof,order2*dof))
kmat   = np.zeros((order2*dof,order2*dof))
x      = np.zeros((2,1))
eye    = np.identity(dof)
Tmat   = np.array([[1,0,0,0],
                   [0,0,1,0],
                   [0,1,0,0],
                   [0,0,0,1]])






# Elements in a structured grid
for i in range(0,order):
  for j in range(0,order):
      elecoord[i*order+j][0] = dxarr[0] + j*harr[0]
      elecoord[i*order+j][1] = dxarr[1] + i*harr[1]



# Get quadrature
wlist,zlist = quadrature.gauss2D(ngpt)

# Compute element stiffness matrix
for i in range(0,len(wlist)):
  zpair = zlist[i]
  N,DN = shapefunctions.shape2D(zpair[0],zpair[1],order)
  J = np.matmul(elecoord.T,DN)
  x = np.matmul(elecoord.T,N.T)
  invJ = np.linalg.inv(J)
  detJ = np.linalg.det(J)
  B = np.matmul(DN,invJ)
  vecB = np.reshape(B,(order2*dof,1))
  kmat = np.matmul(np.kron(B,eye),np.kron(B.T,eye))*mu + np.matmul(np.matmul(np.kron(B,eye),Tmat),np.kron(B.T,eye))*mu
  kmat = kmat + np.matmul(vecB,vecB.T)*lam
  ktotal = ktotal + wlist[i]*kmat*detJ
print(ktotal)