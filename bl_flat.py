import numpy as np
import matplotlib.pyplot as plt
import struct
import Python.Gendata.gendata as gen

##====================================================##
##            >>--regridMaster.py--<<
## 
## Ryan Patmore 17/07/15
## Outputs all required for regridding in MITgcm
##====================================================##

case_name = 'BL_FLAT_001'

ps = None
pe = None
# Set grid
cartesian    = 1
res_multiplier = 1.0 
res          = 1/res_multiplier
ydim         = 180*res_multiplier
xdim         = 180*res_multiplier
#xdim         = int(144/res) 
#ydim         = int(72/res) 
zdim         = 180
latMax       = -1
latMin       = -75
z0           = 20

# Set bathymetry
depMax       = 180
wall_west    = 0
wall_east    = xdim
wall_south   = 0
wall_north   = 0 

# Plateau
plateau_height = 0


# Pressure
p_low  = 1013.5e2
p_high = 1013.4e2

# Temperature
t_low  = 1 
t_high = 2

cd = 0.0012
#cd = 0.00012
rowa = 1.2
       
bathy     = 1
wind      = 0
vels      = 0 
temp      = 1
pressure  = 1

ini_params = { 'xdim'  : xdim,
               'ydim'  : ydim,
               'zdim'  : zdim,
               'depMax': depMax,
               'p_low' : p_low,
               'p_high': p_high,
               't_low' : t_low,
               't_high': t_high }
# MITgcm binary file is saved in (x,y) format

state = gen.State(ini_params)
if temp:
    name = case_name + '_temp.bin'
    t = state.ini_field_linear_grad(state.ydim, state.gridz, t_low, t_high)
    print ('temp shape', t.shape)
    state.writeBin(t,name)
    y = state.readBin(name,int(xdim),int(ydim),int(zdim))
    fig = plt.figure(3)
    p = plt.pcolormesh(y[:,5,:])
    plt.axis('equal')
    plt.colorbar(p)
    plt.title('temp')

if pressure:
    name = case_name + '_pressure.bin'
    p = state.ini_field_linear_grad(state.xdim, state.gridx[:,:,0], p_low, p_high)
    state.writeBin(p,name)
    y = state.readBin(name,int(xdim),int(ydim))
    fig = plt.figure(4)
    p = plt.pcolormesh(y)
    plt.axis('equal')
    plt.colorbar(p)
    plt.title('pressure')


if vels:
    name = 'Binary/SPBC_303_uvel.bin'
    vels = Velocities('SPBC_303')
    vels.extract_velocities('SPBC_286')
    y = bathy.readBin(name,xdim,ydim)
    bathy.plot_single(y)
    
if bathy:
    name = case_name + '_bathy.bin'
    print (name)
    bathy = gen.Bathymetry(ini_params)
    print ('bathy', bathy)
    b = bathy.get_bathy()#[::-1]
    bathy.writeBin(b,name)
    y = bathy.readBin(name,int(xdim),int(ydim))
    fig = plt.figure(0)
    p = plt.pcolormesh(y)
    plt.axis('equal')
    plt.colorbar(p)
    plt.title('bathy')
    

if wind:
    wind_name = case_name + '_wind.bin'
    #print wind_name
    w = Wind()
    #w.ydim = 144
    #w.wind()
    w.shrunk_wind()
    wind = w.get_wind()
    #wind = np.hstack( (wind,np.zeros((288,144))) )
    #wind = np.hstack( (wind,-wind) )
    #wind = wind - wind.min()
    
    #wind = w.readBin(wind_name, xdim, ydim)
    #w.plot_single(wind)
    #wind[:,144:] = np.zeros((288,144))
    #w.plot_num = 1
    
    w.writeBin(wind, wind_name)
    wind = w.readBin(wind_name, int(xdim), int(ydim))
    w.case = '001'
    plt.figure(101)
    p = plt.pcolormesh(wind)
    plt.colorbar(p)
    #plt.figure(100)
    #plt.plot(wind[0,:])
    plt.show()

#plt.show()
