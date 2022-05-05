import numpy as np
import math
def gauss2D(ngpt):
  if (ngpt == 1):
    zgpt = np.array([0.0])
    wgpt = np.array([2.0])
  elif (ngpt == 2):
    zgpt = np.array([-1/math.sqrt(3),1/math.sqrt(3)])
    wgpt = np.array([1.0,1.0])
  elif (ngpt == 3):
    zgpt = np.array([-0.744596669,0.0,0.744596669])
    wgpt = np.array([5.0/9.0,8.0/9.0,5.0/9.0])
  elif (ngpt == 4):
    zgpt = np.array([-0.861136312,-0.339981044, 0.339981044, 0.861136312])
    wgpt = np.array([0.3478548,0.6521452,0.6521452,0.3478548])
  elif (ngpt == 5):
    zgpt = np.array([-0.906179846,-0.538469310,0.0, 0.538469310,0.906179846])
    wgpt = np.array([0.2369269,0.4786287,0.5688889,0.4786287,0.2369269])
  elif (ngpt == 6):
    zgpt = np.array([-0.932469514,-0.661209386,-0.238619186, 0.238619186, 0.661209386, 0.932469514])
    wgpt = np.array([0.1713245,0.3607616,0.4679139,0.4679139,0.3607616,0.1713245])
  else:
    print("Invalid ngpt number, must be between 1 and 6")
    exit(1)

  wlist = np.zeros((ngpt*ngpt))
  zlist = np.zeros((ngpt*ngpt,2))
  for i in range(0,ngpt):
    for j in range(0,ngpt):
        wlist[i*ngpt+j] = wgpt[i]*wgpt[j]
        zlist[i*ngpt+j][0] = zgpt[j]
        zlist[i*ngpt+j][1] = zgpt[i]
  return wlist,zlist
  
def gauss3D(ngpt):
  if (ngpt == 1):
    zgpt = np.array([0.0])
    wgpt = np.array([2.0])
  elif (ngpt == 2):
    zgpt = np.array([-1/math.sqrt(3),1/math.sqrt(3)])
    wgpt = np.array([1.0,1.0])
  elif (ngpt == 3):
    zgpt = np.array([-0.744596669,0.0,0.744596669])
    wgpt = np.array([5.0/9.0,8.0/9.0,5.0/9.0])
  elif (ngpt == 4):
    zgpt = np.array([-0.861136312,-0.339981044, 0.339981044, 0.861136312])
    wgpt = np.array([0.3478548,0.6521452,0.6521452,0.3478548])
  elif (ngpt == 5):
    zgpt = np.array([-0.906179846,-0.538469310,0.0, 0.538469310,0.906179846])
    wgpt = np.array([0.2369269,0.4786287,0.5688889,0.4786287,0.2369269])
  elif (ngpt == 6):
    zgpt = np.array([-0.932469514,-0.661209386,-0.238619186, 0.238619186, 0.661209386, 0.932469514])
    wgpt = np.array([0.1713245,0.3607616,0.4679139,0.4679139,0.3607616,0.1713245])
  else:
    print("Invalid ngpt number, must be between 1 and 6")
    exit(1)

  wlist = np.zeros((ngpt*ngpt*ngpt))
  zlist = np.zeros((ngpt*ngpt*ngpt,3))
  for i in range(0,ngpt):
    for j in range(0,ngpt):
      for k in range(0,ngpt):
        wlist[i*ngpt*ngpt+j*ngpt+k] = wgpt[i]*wgpt[j]*wgpt[k]
        zlist[i*ngpt*ngpt+j*ngpt+k][0] = zgpt[k]
        zlist[i*ngpt*ngpt+j*ngpt+k][1] = zgpt[j]
        zlist[i*ngpt*ngpt+j*ngpt+k][2] = zgpt[i]
  return wlist,zlist