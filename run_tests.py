#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Program to test the modules to be used in inversions ofsubsurface
temperature profiles.

Francisco Jose Cuesta-Valero
2021-06-23
"""
import os
import numpy as np
import matplotlib as mpl
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon


# Flags
flags = {
  'svd'  : True,
  'ppi'  : True,
  'bti'  : True,
  'plot' : True
}

# -----------------------------------------------------------------------------
# Functions
def main(flags):
  # Run SVD test
  if flags['svd']:
    run_test_svd()

  # Run PPI test
  if flags['ppi']:
    run_test_ppi()

  # Run BTI test
  if flags['bti']:
    run_test_bti()

  # Plot the relevant results
  if flags['plot']:
    plot_inversions(flags)

  return


def run_test_bti():
  os.system('cp data/profile.txt profile_10.txt')
  os.system('cp data/profile.txt profile_11.txt')
  fmodules_test = '-fopenmp src/mod_kinds.f03 src/mod_logger.f03' \
      ' src/mod_regression.f03 src/mod_profile.f03' \
      ' src/mod_fm.f03 src/mod_svd.f03 src/mod_ensembles.f03' \
      ' src/mod_ghf.f03 src/mod_bti.f03'

  script_test = "test_bti.f03"
  run_fortran(fmodules_test, script_test)
  os.remove('profile_10.txt')
  os.remove('profile_11.txt')

  return

def run_test_ppi():
  os.system('cp data/profile.txt profile.txt')
  fmodules_test = 'src/mod_kinds.f03 src/mod_logger.f03' \
      ' src/mod_regression.f03 src/mod_profile.f03' \
      ' src/mod_fm.f03 src/mod_svd.f03 src/mod_ppi.f03'

  script_test = "test_ppi.f03"
  run_fortran(fmodules_test, script_test)
  os.remove('profile.txt')

  return

def run_test_svd():
  os.system('cp data/profile.txt profile.txt')
  fmodules_test = 'src/mod_kinds.f03 src/mod_logger.f03' \
      ' src/mod_regression.f03 src/mod_profile.f03' \
      ' src/mod_fm.f03 src/mod_svd.f03'

  script_test = "test_svd.f03"
  run_fortran(fmodules_test, script_test)
  os.remove('profile.txt')

  return

def plot_inversions(flags):
  tmin = 1600
  tmax = 2021

  if flags['svd']:
    file_svd = 'inv_svd.dat'
    color_svd = (0/255, 0/255, 204/255)
    alpha_svd = 0.3
    svd = read_inversion(file_svd,'temp','svd',tmin,tmax)
    psvd = create_polygon(svd,color_svd,alpha_svd)

  if flags['ppi']:
    file_ppi = 'inv_ppi.dat'
    color_ppi = (127/255, 127/255, 127/255)
    alpha_ppi = 0.3
    ppi = read_inversion(file_ppi,'temp','ppi',tmin,tmax)
    pppi = create_polygon(ppi,color_ppi,alpha_ppi)

  if flags['bti']:
    file_bti = 'inv_bti.dat'
    color_bti = (204/255, 0/255, 0/255)
    alpha_bti = 0.2
    bti = read_inversion(file_bti,'temp','bti',tmin,tmax)
    pbti = create_polygon(bti,color_bti,alpha_bti)


  fig = plt.figure()
  fig, ax = plt.subplots(figsize=(6,3))

  if flags['svd']:
    ax.add_patch(psvd)

  if flags['ppi']:
    ax.add_patch(pppi)

  if flags['bti']:
    ax.add_patch(pbti)

  if flags['svd']:
    ax.plot(svd[:,0],svd[:,2],color=color_svd,label='SVD')

  if flags['ppi']:
    ax.plot(ppi[:,0],ppi[:,2],color=color_ppi,label='PPI')

  if flags['bti']:
    ax.plot(bti[:,0],bti[:,2],color=color_bti,label='BTI')


  ax.set_xlabel('Years C.E.')
  ax.set_ylabel('Temperature Change ($^{\circ}$C)')

  ax.legend(loc='upper left')

  fig.tight_layout()
  fig.savefig('test_fig.pdf',format='pdf')


def run_fortran(fmodules,script):
  os.system("gfortran "+fmodules+" "+script+" -o script.out")
  os.system("./script.out > logfor.log")
  with open('logfor.log') as myfile:
    if 'Status = 1\n' not in myfile.read():
      print('STOP SCRIPT '+script)
      raise SystemExit

  os.remove('logfor.log')
  os.system('rm *.mod')
  os.remove('script.out')


def read_inversion(file,var,method="bti",ymin=0,ymax=0):
  if method == 'bti':
    key = '#g'+var
  else:
    key = '#g'

  data=open(file)
  data = data.readlines()
  temp1 = []
  temp2 = []
  temp3 = []
  time = []
  for i in range(len(data)):
    row=data[i].split()
    if len(row) > 0:
      if method == 'dat':
        time.append(float(row[0]))
        temp1.append(float(row[1]))
        temp2.append(float(row[2]))
        temp3.append(float(row[3]))
      else:
        if row[0] == key:
          time.append(float(row[1]))
          temp1.append(float(row[2]))
          temp2.append(float(row[3]))
          temp3.append(float(row[4]))

  time = np.array(time)
  temp1 = np.array(temp1)
  temp2 = np.array(temp2)
  temp3 = np.array(temp3)

  if ymin == 0 & ymax == 0:
    n = len(time)
    results = np.ndarray(shape=(n,4))
    results[:,0] = time[:]
    results[:,1] = temp1[:]
    results[:,2] = temp2[:]
    results[:,3] = temp3[:]
  else:
    ind = (time <= ymax) & (time >= ymin)
    n = sum(ind)
    results = np.ndarray(shape=(n,4))
    results[:,0] = time[ind]
    results[:,1] = temp1[ind]
    results[:,2] = temp2[ind]
    results[:,3] = temp3[ind]


  return results

def create_polygon(data,color,alpha):
  verts = np.ndarray(shape=(data.shape[0]*2,2))
  verts[0:data.shape[0],0] = data[:,0]
  verts[data.shape[0]:data.shape[0]*2,0] = data[::-1,0]
  verts[0:data.shape[0],1] = data[:,1]
  verts[data.shape[0]:data.shape[0]*2,1] = data[::-1,3]
  p = Polygon(verts,color=color,alpha=alpha,ec=None)

  return p

# -----------------------------------------------------------------------------

main(flags)



