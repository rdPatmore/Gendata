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

case_name = 'SHIFT_208'

ps = None
pe = None
# Set grid
cartesian    = 1
res_multiplier = 1.0 
res          = 1/res_multiplier
zRes         = 1
ydim         = int(1*res_multiplier)
xdim         = int(10*res_multiplier)
zdim         = 5
latMax       = -1
latMin       = -75
z0           = 20

# Set bathymetry
depMax       = 4
wall_west    = 0
wall_east    = xdim
wall_south   = 0
wall_north   = 0 
 
iceDep = 10

# EOS
tAlpha = 3.9e-5
sBeta = 7.41e-4
rhoConst = 1030

g = 9.81

# Temperature
t_low  = 1.0000001 
t_high = 1

# Salinity
s_low  = 34.0 
s_high = 34.5

# Pressure
# suposed realistic values
#p_low  = 1000e2
#p_high = 4000e2
#etaGrad = 0.0001
p_low  = 1000e2
p_high = 4000e2
etaGrad = 0.000000001

# Velocity
v_low = 1
v_high = 2

cd = 0.0012
#cd = 0.00012
rowa = 1.2
       
bathy      = 0
wind       = 0
vels       = 0
temp       = 0
salt       = 0
pressure   = 1
shice      = 0
shice_0    = 0



ini_params = { 'xdim'  : xdim,
               'ydim'  : ydim,
               'zdim'  : zdim,
               'g'     : g,
               'sBeta' : sBeta,
               'tAlpha': tAlpha,
               'rhoConst': rhoConst,
               'depMax': depMax,
               'iceDep': iceDep,
               's_low' : s_low,
               's_high' : s_high,
               't_low' : t_low,
               't_high': t_high }
# MITgcm binary file is saved in (x,y) format

state = gen.State(ini_params)
if temp:
    name = case_name + '_ini_temp.bin'
    #t0 = state.add_heat_blob(state.gridx, state.gridz, 90, 170, 15)
    #t0 = state.ini_field_hill(state.xdim, state.gridx, 0.1)
    t = state.ini_field_linear_grad(state.zdim, state.gridz, t_low, t_high)
    print ('tshape', t[4,0,0])
    #t[:2,:,:] = t[0,0,0]
    print ('tshape', t[4,0,0])
    #t = t0 + t1
    print (state.shape)
    #t = np.full(state.shape, 1)
    
    #t[2,0,-3] = 2
    #t[0,0,-3] = 2
    #t[1,0,-3] = 2
    print ('state.shape', state.shape)
    print ('gridz.shape', state.gridx.shape)
    state.writeBin(t,name)
    y = state.readBin(name,int(xdim),int(ydim),int(zdim))
    fig = plt.figure(3)
    p = plt.pcolormesh(np.squeeze(y))
    plt.axis('equal')
    plt.colorbar(p)
    plt.title('temp')


if salt:
    name = case_name + '_ini_salt.bin'
    s0 = state.ini_field_hill(state.xdim, state.gridx, 0.05)
    s1 = state.ini_field_linear_grad(state.zdim, state.gridz, t_low, t_high)
    s = s0 + s1
    state.writeBin(s,name)
    y = state.readBin(name,int(xdim),int(ydim),int(zdim))
    fig = plt.figure(4)
    p = plt.pcolormesh(np.squeeze(y))
    plt.axis('equal')
    plt.colorbar(p)
    plt.title('temp')

if pressure:
    name = case_name + '_pForceX.bin'
    #pload = state.ini_p_force(state.xdim, state.gridx, p_low, p_high)
    #pload = pload * state.gridz[:,:,::-1] / state.gridz.max()
    pload = np.full(state.gridx.shape, rhoConst * g * etaGrad) # uniform pforce
    #pload = np.full(state.shape, 0)
    
    #pload[2,0,-2] = 0.001
    state.writeBin(pload,name)
    y = state.readBin(name,int(xdim),int(ydim),int(zdim))
    fig = plt.figure(4)
    print ('pressure shape', y.shape)
    p = plt.pcolormesh(y[:,0,:])
    plt.axis('equal')
    plt.colorbar(p)
    plt.title('pressure')

if shice_0:
    ShiceTopo  = case_name + '_ini_shice_topo.bin'
    #shice_topo = - np.full((state.xdim, state.ydim), iceDep)
    shice_topo = np.zeros((int(ydim),int(xdim)))
    shice_topo[0,:] = -1
    shice_topo[0,:2] = -2
    state.writeBin(shice_topo, ShiceTopo)
    y = state.readBin(ShiceTopo,int(xdim),int(ydim))
    fig = plt.figure(5)
    p = plt.pcolormesh(y)
    plt.axis('equal')
    plt.colorbar(p)
    plt.title('ini_shice')


if shice:
    # Check for initial state files
    if temp == 0:
        t = np.full(state.gridx.shape, 20)
    if salt == 0:
        s = np.full(state.gridx.shape, 30)

    # Set file names
    ShicePFile = case_name + '_ini_shice_p.bin'
    ShiceTopo  = case_name + '_ini_shice_topo.bin'
    
    # Make setup files
    iceProfile = np.zeros((ydim,xdim))
    iceProfile[0,:] = -1
    iceProfile[0,:2] = -2
    hFacC = state.readBin('../SHIFT_002_hFacC.data', x=int(xdim), y=int(ydim),
                                                     z=int(zdim))
    #hFacC = np.ones((zdim,ydim,xdim))

    # Create MITgcm input files
    shice_topo, shice_p = state.ini_shice(iceProfile, t, s, hFacC)
    state.writeBin(shice_p, ShicePFile)
    state.writeBin(shice_topo, ShiceTopo)

    # Plot binary
    y = state.readBin(ShicePFile,int(xdim),int(ydim))
    fig = plt.figure(5)
    p = plt.pcolormesh(y)
    plt.axis('equal')
    plt.colorbar(p)
    plt.title('shice')


if vels:
    name = case_name + '_ini_uvel.bin'
    #vels = state.add_heat_blob(state.gridx, state.gridz, 90, 100, 15)
    #vels = state.ini_field_linear_grad(state.zdim, state.gridz, v_low, v_high)
    vels = np.full(state.shape, 0.0000001)
    
    #vels[2,0,-3] = 1
    state.writeBin(vels,name)
    y = state.readBin(name,int(xdim),int(zdim))
    fig = plt.figure(7)
    p = plt.pcolormesh(y)
    plt.axis('equal')
    plt.colorbar(p)
    plt.title('vels')
    
if bathy:
    name = case_name + '_bathy.bin'
    print (name)
    bathy = gen.Bathymetry(ini_params)
    print ('bathy', bathy)
    b = bathy.get_bathy()#[::-1]
    b[:,8:] = -3
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
    wind = w.readBin(wind_name, int(ydim), int(xdim))
    w.case = '001'
    plt.figure(101)
    p = plt.pcolormesh(wind)
    plt.colorbar(p)
    #plt.figure(100)
    #plt.plot(wind[0,:])
    plt.show()

plt.show()
