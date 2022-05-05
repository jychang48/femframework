import numpy as np

def l2element(z):
  N = np.zeros((2))
  DN = np.zeros((2))
  N[0] = 0.5*(1-z)
  N[1] = 0.5*(1+z)
  DN[0] = -0.5
  DN[1] = 0.5
  return N,DN

def l3element(z):
  N = np.zeros((3))
  DN = np.zeros((3))
  N[0] = 0.5*z*(1-z)
  N[1] = (1-z*z)
  N[2] = 0.5*z*(1+z)
  DN[0] = 0.5*(2*z-1)
  DN[1] = -2*z
  DN[2] = 0.5*(2*z+1)
  return N,DN

def q4element(zx,zy):
  Nx,DNx = l2element(zx)
  Ny,DNy = l2element(zy)
  N = np.zeros(4)
  DN = np.zeros((4,2))
  for i in range(0,2):
    for j in range(0,2):
      N[i*2+j] = Nx[j]*Ny[i]
      DN[i*2+j][0] = DNx[j]*Ny[i]
      DN[i*2+j][1] = Nx[j]*DNy[i]
  return N,DN

def q9element(zx,zy):
  Nx,DNx = l3element(zx)
  Ny,DNy = l3element(zy)
  N = np.zeros(9)
  DN = np.zeros((9,2))
  for i in range(0,3):
    for j in range(0,3):
      N[i*3+j] = Nx[j]*Ny[i]
      DN[i*3+j][0] = DNx[j]*Ny[i]
      DN[i*3+j][1] = Nx[j]*DNy[i]
  return N,DN

def b8element(zx,zy,zz):
  Nx,DNx = l2element(zx)
  Ny,DNy = l2element(zy)
  Nz,DNz = l2element(zz)
  N = np.zeros(8)
  DN = np.zeros((8,3))
  for i in range(0,2):
    for j in range(0,2):
      for k in range(0,2):
        N[i*4+j*2+k] = Nx[k]*Ny[j]*Nz[i]
        DN[i*4+j*2+k][0] = DNx[k]*Ny[j]*Nz[i]
        DN[i*4+j*2+k][1] = Nx[k]*DNy[j]*Nz[i]
        DN[i*4+j*2+k][2] = Nx[k]*Ny[j]*DNz[i]
  return N,DN

def b27element(zx,zy,zz):
  Nx,DNx = l3element(zx)
  Ny,DNy = l3element(zy)
  Nz,DNz = l3element(zz)
  N = np.zeros(27)
  DN = np.zeros((27,3))
  for i in range(0,3):
    for j in range(0,3):
      for k in range(0,3):
        N[i*9+j*3+k] = Nx[k]*Ny[j]*Nz[i]
        DN[i*9+j*3+k][0] = DNx[k]*Ny[j]*Nz[i]
        DN[i*9+j*3+k][1] = Nx[k]*DNy[j]*Nz[i]
        DN[i*9+j*3+k][2] = Nx[k]*Ny[j]*DNz[i]
  return N,DN

def shape1D(z,order):
  if (order == 2):
    N,DN = l2element(z)
  elif (order == 3):
    N,DN = l3element(z)
  else:
    print("Invalid order, must be 2 or 3")
    exit(1)
  return N,DN

def shape2D(zx,zy,order):
  if (order == 2):
    N,DN = q4element(zx,zy)
  elif (order == 3):
    N,DN = q9element(zx,zy)
  else:
    print("Invalid order, must be 2 or 3")
    exit(1)
  return N,DN
  
def shape3D(zx,zy,zz,order):
  if (order == 2):
    N,DN = b8element(zx,zy,zz)
  elif (order == 3):
    N,DN = b27element(zx,zy,zz)
  else:
    print("Invalid order, must be 2 or 3")
    exit(1)
  return N,DN