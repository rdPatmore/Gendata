import numpy as np
import matplotlib.pyplot as plt
import struct
import config

##====================================================##
##            >>--regridMaster.py--<<
## 
## Ryan Patmore 17/07/15
## Outputs all required for regridding in MITgcm
##====================================================##

if __name__ == '__main__':
    case_name = 'RDP_032'
    
    ps = None
    pe = None

    # Set grid
    cartesian    = 1
    res_multiplier = 2.0 
    res          = 25000/res_multiplier
    ydim         = 288*res_multiplier
    xdim         = 288*res_multiplier
    latMax       = -1
    latMin       = -75
    z0           = 20
    
    # Set bathymetry
    depMax       = 5000
    wall_west    = 0
    wall_east    = xdim
    wall_south   = 1
    wall_north   = -1
    
    # Plateau
    plateau_height = 0
    
    # Ridge
    ridge        = 1 # 1 for ridge 0 for no ridge
    ridge_height  = 2000
    #ridge_width   = int(5/res)
    ridge_width   = 10*res_multiplier
    e_width = 10*res_multiplier # width of eastern ridge slope (sawtooth ridge)
    w_width = 10*res_multiplier # width of western ridge slope (sawtooth ridge)
    
    # Gap
    gap = 0
    #ridge_centre = 0	 
    v_centre = ydim/2.0 
    #v_centre = 0.0 #y centre of rotation for ridge
    ridge_centre = int(0.5*xdim)
    #ridge_centre = 80*res_multiplier
    print ('RIDGE CENTRE', ridge_centre)
    y_0          = (ydim*(2.0/5.0))- ridge_width
    grad = 0.0 
    
    # Wind
    windStrength = 2
    windStrength = 10
    curl         = 0 
    zero_curl    = 1 #shifts point of zero curl n & s (n x>1, s x<1)
    
    # Density
    sigma_0      = 999.8
    tAlpha       = 0.0002
    sigmaMin     = 26.1
    sigmaMax     = 28.35
    rhoDiff = 1
    
    cd = 0.0012
    #cd = 0.00012
    rowa = 1.2

class Binary(object):

    def writeBin(self, inputArray, outFileName):
        w = np.float64(inputArray).flatten()
        with open(self.writePath + outFileName,'wb') as f:
            for i in w:
                f.write(struct.pack('>d', i))   
    #def writeBin(self, inputArray, outFileName):
    #    print ('INPUT ARRAY SHAPE', inputArray.shape)
    #    inputArray.tofile(outFileName)

    def readBin(self, file, x=1, y=1, z=None, dtype='float64'):
        print ('PATH', self.writePath)
        print ('FILE', file)
        print ('dtype', dtype)
        file = np.fromfile(self.writePath + file,dtype=dtype).byteswap()
        print (file.shape)
        if z == None:
            matrix = file.reshape((y,x))
        else:
            matrix = file.reshape((z,y,x))
        print ('MATRICX SHAPE', matrix.shape)
        return matrix

class Grid(Binary):
    
    def __init__(self, ini_params):
        self.writePath = 'Binary/'
        self.ini_params = ini_params
        self.xdim = int(ini_params['xdim'])
        self.ydim = int(ini_params['ydim'])
        self.zdim = ini_params['zdim']
        #self.ps = ps
        #self.pe = pe
        self.grid()

    def grid(self):
        '''define coordinate axes'''
        print ('XDIM:', self.xdim)
        print ('YDIM:', self.ydim)
        self.gridz, self.gridy, self.gridx = np.meshgrid(np.arange(self.zdim), 
                                                         np.arange(self.ydim),
                                                         np.arange(self.xdim),
                                                         indexing='ij')
        print ('XDIM:', self.gridx.shape)
        print ('YDIM:', self.gridy.shape)


    def sheer_coordinates(self):
        '''coordinate transform that sheers the y-axis through angle theta'''
        
        theta = 15*np.pi/32

        xdash = self.gridx + (1 / np.tan(theta)) * self.gridy
        ydash = self.gridy

        self.gridx = xdash
        self.gridy = ydash


    def shift_coordinates(self):
        '''coordinate transform that shifts the x and y axes to be ridge 
           centred'''
         
        self.gridx -=  ridge_centre
        self.gridy -=  v_centre
        
  
    def latGrid(self):
        arr  = [-latMin for i in range(self.ydim+1)]
        DYC = [-latMin for i in range(self.ydim)]
        for i in range(self.ydim):
            DYC[i] = res*(lon2m(arr[i])/lat2m(arr[i]))
            arr[i+1] = arr[i]-DYC[i]
        self.DYC = DYC
        self.Y = arr

    def zGrid(self):
        gridSum = depMax+10.0
        depiter = 400.0 
        while gridSum > depMax:
            c         = (self.zdim-1)/np.log(depiter/z0)    
            zGrid     = [z0*np.exp(i/c) for i in range(self.zdim)]
            gridSum = np.sum(zGrid)
            depiter = depiter-1.0
        zGridNorm = (np.cumsum(zGrid)/np.sum(zGrid))*self.zdim
        zGridCum  = (np.cumsum(zGrid))
        self.drF  = np.array(zGrid),
        self.depth = zGridCum

class Bathymetry(Grid):

    def __init__(self, ini_params):
        Grid.__init__(self, ini_params)
        try:
            self.rw    = ridge_centre - ridge_width
            self.re    = ridge_centre + ridge_width
            self.rh    = ridge_height
        except:
            print ('no ridge')
        self.bathy = np.zeros((self.gridx[0,:,:].shape))
        self.bathy[:,:] = - ini_params['depMax'] 
        self.set_title = 'Bathymetry'
        self.set_save = 'Bathymetry'
           
    def get_bathy(self):
        return self.bathy


    def meridional_wall(self):
        self.bathy[self.gridx==0] = 0.0
    

    def ridge(self, north_atlantic=0):
        '''Generates a gaussian ridge'''

        msk = np.zeros((self.bathy.shape))

        msk[(self.gridx - ridge_width < 0) & (self.gridx + ridge_width > 0)] = 1
        self.bathy += msk * (
                      np.cos(self.gridx * np.pi/ridge_width) + 1) * self.rh / 2
        
        if north_atlantic: 
            for i in range(int(ydim)):
                self.bathy[np.argmax(self.bathy[:,i]),i] = 0

    def sawtooth_ridge(self, peak_wall=0):
        '''
        Generates a gaussian ridge with differeing gradients either side
           
        peak_wall add 0.0 metre bathymetry to ridge peak.
        '''

        msk1 = np.zeros((self.bathy.shape))
        msk2 = np.zeros((self.bathy.shape))
        print ('msk shape', msk1.shape)
        print ('', msk1.shape)

        msk1[(self.gridx - w_width < 0) & (self.gridx > 0)] = 1
        msk2[(self.gridx < 0) & (self.gridx + e_width > 0)] = 1

    
        self.bathy += msk1 * (
                      np.cos(self.gridx * np.pi/w_width) + 1) * self.rh / 2
        self.bathy += msk2 * (
                      np.cos(self.gridx * np.pi/e_width) + 1) * self.rh / 2

        peak = 0.0 if peak_wall == 1 else self.rh - depMax
        self.bathy[self.gridx==0] = peak

    def boundary_slope(self):
        '''
        Adds a ridge that spans a preiodic boundary.
        '''
        
        msk1 = np.zeros((self.bathy.shape))
        msk2 = np.zeros((self.bathy.shape))
        #gridy, gridx = np.meshgrid(np.arange(msk1.shape[1]),
        #                           np.arange(msk2.shape[0]))

        msk1[self.gridx - ridge_width < 0] = 1
        self.pre_bathy = self.bathy.copy()
        self.bathy= self.bathy+ msk1*((np.cos( self.gridx*np.pi/(ridge_width) )
                                    * (self.rh/2)) + (self.rh/2))
        msk2[self.gridx[:,::-1] - ridge_width < 0] = 1
        self.bathy= self.bathy+ msk2*((np.cos( self.gridx*np.pi/ridge_width )
                                    * (self.rh/2)) + (self.rh/2))[:, ::-1]
        print (self.bathy[:,10])
        

    def gap(self):
        y_0 = (ydim/2.0)-ridge_width
        y_1 = (ydim/2.0)+ridge_width
        self.bathy[self.rw:self.re,y_0:y_1] = -depMax

    def north_south_boundary(self):
        self.bathy[int(-2*res_multiplier):,:] = 0.0
        self.bathy[:int(2*res_multiplier),:] = 0.0

    def east_west_boundary(self):
        self.bathy[-1:,:] = 0.0
        self.bathy[:1,:] = 0.0
    
    def chop_ridge_end(self):


        shift = self.shift
        ridge_end = ydim/2.0
        self.ridge_end= ridge_end
        msk = np.zeros((self.bathy.shape))
        msk_end = np.zeros((self.bathy.shape))
        gridy, gridx = np.meshgrid(np.arange(msk.shape[1]),
                                   np.arange(msk.shape[0]))
        line      = gridx - (grad * (gridy + v_centre))
        #perp_line = (grad*gridx+v_centre) + gridy
        line = line - line[0,0]- line[ridge_centre, 0]
        perp_line = line[:,::-1]
        perp_line = perp_line -\
                  perp_line[np.argmin(np.abs(line)[:,ridge_end]), ridge_end]
        msk[(line+shift>0) & (line-shift<0) & (perp_line<0)] = 1
        msk_end[(line+shift>0) & (line-shift<0) &
                (perp_line>=0) & (perp_line-shift<0)] = 1
        bump = (1+np.cos( (line)*np.pi/shift )) *\
              ((1+np.cos( (perp_line)*np.pi/shift )))* (self.rh/4)

        self.bathy = self.pre_bathy + (self.bathy+5000.0)*msk
        self.bathy = self.bathy + (msk_end * bump)

    def plateau_slope(self):
        msk = np.zeros((self.bathy.shape))
        gridy, gridx = np.meshgrid(np.arange(msk.shape[1]),
                                   np.arange(msk.shape[0]))
        line = gridx - (grad * gridy) - ridge_centre
        xlimit = 240*res_multiplier
        ylimit = self.ridge_end/2.0
        msk[(line>=0) & (gridy>1) & (gridx<xlimit) & (gridy<ylimit)] = 1
        plt.figure(11)
        plt.pcolor(gridx, gridy, msk)
        
        pre_plateau = (res_multiplier*241.0/(gridx+1.0)) - 1.0
        plateau = (1500.0*pre_plateau/(msk*pre_plateau).max()) - depMax
        plateaued = plateau*msk
        plateaued[plateaued == 0] = -depMax
        self.bathy = np.where(self.bathy >= plateaued, self.bathy, plateaued)
   
    def smooth_gap(self, orientation = 1):     
        sigma =10
        y_0 = (ydim/2.0)-(orientation*ridge_width+sigma/2)
        w = ridge_width 
        indy = np.arange(sigma)[::orientation]
        xlen = ((self.re-self.rw)/2.0)
        indx = np.arange(-sigma,sigma)
        self.bathy[self.rw:self.re,y_0:y_0+sigma] =\
        ((self.rh-(5000-depMax))/4)*(1+np.cos(np.arange(-xlen,xlen)*np.pi/(xlen)))[:,None]*\
                    (1+np.cos(indy*np.pi/xlen))[None,:] - depMax
        #self.bathy[self.rw:self.re,y_0:y_0+sigma] =\
        #(self.rh/4)*(1+np.cos(np.arange(-xlen,xlen)*np.pi/ridge_width))[:,None]*\
        #            (1+np.cos(indy*np.pi/ridge_width))[None,:] - depMax
  
    def plateau(self):
        intersection = max(np.where(self.bathy[:,30] > -4500.0)[0])
        self.bathy[intersection:intersection + 50,2:(ydim/2.0) - (ridge_width+5)] = self.bathy[intersection,3]

class Wind(Grid):
    
    def __init__(self):
        Grid.__init__(self)
        self.curl = curl
        self.wind_strength = windStrength
        self.set_save = 'wind'
        self.set_title = 'Wind'

    def tau(self, U):
        tau = cd*rowa*U**2
        return tau

#        wind1d = 0.5*self.tau(self.wind_strength)*
#       (1-self.curl*np.cos(2*np.pi*np.power((np.arange(float(self.ydim))/
#             self.ydim),zero_curl)))[None].T
#        self.wind = np.tile(wind1d,self.xdim)
#        return self.wind.T
    def wind_reverse(self):
        wind1d    = self.tau(self.wind_strength) *\
                    -np.sin(np.pi * np.arange(self.ydim)/self.ydim)[None].T
        self.wind = np.tile(wind1d, self.xdim)
        return self.wind.T
    def wind(self):
        #wind1d    = self.tau(self.wind_strength) *\
        #       np.sin(((np.pi * np.arange(self.ydim))+ (self.ydim *np.pi/2.0))/self.ydim )[::-1][None].T
        print ('WINDDDDD', self.tau(self.wind_strength))
        wind1d    = self.tau(self.wind_strength) *\
               np.sin(2*np.pi * np.arange(self.ydim) /self.ydim )[None].T
        self.wind = np.tile(wind1d, self.xdim)
        return self.wind.T
  
    def constant_wind(self):
        wind1d = self.tau(self.wind_strength)*np.ones(self.ydim)
        self.wind = np.tile(wind1d, self.xdim)
        return self.wind.T
    
    def constant_wind_west_only(self):
        res = res_multiplier
        wind1d = self.tau(self.wind_strength) * np.ones((self.xdim/2)-(100*res))
        print (wind1d)
        transition2 = 0.5*self.tau(self.wind_strength) * \
                     (1+np.cos(np.pi * np.arange(50*res)/(50*res)))
        transition1 = transition2[::-1]
        zeros = np.zeros(self.xdim/2)
        xwind = np.concatenate((transition1, wind1d, transition2, zeros))
        print('TYPE',type(xwind))
        print('SGHAPE',xwind.shape)
        self.wind = np.tile(xwind, (self.ydim))
        return self.wind
 
    def constant_wind_plateau(self):
        res = int(res_multiplier)
        wind1d = self.tau(self.wind_strength) * np.ones(50*res)
        print (wind1d.shape)
        transition = (1 - (np.arange(10*res)/(10*res))
                     ) * self.tau(self.wind_strength)
        zeros = np.zeros(self.xdim - (120 * res))
        xwind = np.concatenate((wind1d, transition, zeros, transition[::-1],
                                                           wind1d))
        print('SHAPE', xwind.shape)
        north  = np.zeros((self.xdim, int(self.ydim * 1 / 2)))
        south  = np.zeros((self.xdim, int(self.ydim * 1 / 4)))
        plateau = np.tile(xwind, (int(self.ydim/4), 1)).T
        print ('south shape', south.shape)
        print ('north shape', north.shape)
        print ('plateau shape', plateau.shape)
        self.wind = np.concatenate((south, plateau, north), axis=1).T
        print('SHAPE', self.wind.shape)
        return self.wind
        

    def shrunk_wind(self, half_y=1):
        print ('ydim', self.ydim)
        if half_y:
            ydim_w = (self.ydim) - (4.0 * res_multiplier) 
        else:
            ydim_w = (self.ydim/2.0) - (4.0 * res_multiplier) 
        #ydim_w = self.ydim - (4.0 * res_multiplier) 
        #wind1d    = self.tau(self.wind_strength) *\
        #       np.sin(np.pi * np.arange(ydim_w) / ydim_w )
        wind1d    = 0.5 * self.tau(self.wind_strength) *\
               (1 + np.cos(np.pi + (2*np.pi * np.arange(ydim_w) / ydim_w )))
        pad1 = np.zeros((int(2 * res_multiplier)))
        if half_y:
            pad2 = np.zeros(( int(2 * res_multiplier) ))
        else:
            pad2 = np.zeros(( int((2 * res_multiplier) + (self.ydim/2.0)) ))
        wind1d = np.concatenate([pad1, wind1d, pad2], axis=0)
        print ('wind1d', wind1d.shape)
        self.wind = np.tile(wind1d, (self.xdim, 1))
        print ('wind', self.wind.shape)
        

    def shrunk_wind_test(self):
        ydim_w = self.ydim - (150.0 * res_multiplier) 
        wind1d    = self.tau(self.wind_strength) *\
               np.sin(np.pi * np.arange(ydim_w) / ydim_w )
        pad1 = np.zeros((75 * res_multiplier))
        pad2 = np.zeros(( (75 * res_multiplier) ))
        wind1d = np.concatenate([pad1, wind1d, pad2], axis=0)[None].T
        self.wind = np.tile(wind1d, self.xdim)

    def no_wind_north(self):
        self.wind[self.ydim/2.0:] = 0.0
    
    def get_wind(self):
        return self.wind.T

class Restoring(Grid):
    
    def __init__(self):
        Grid.__init__(self)
        self.zGrid()
        rhoDiff = 1 
        self.set_title= 'sal'
        self.set_save= 'sal'
        lat  = res*np.arange(ydim)
        self.DX, self.DY = np.meshgrid(lat,self.drF)
        
    def salt(self):
        self.zSal()
        northVariance = np.linspace(0,rhoDiff,self.ydim)
        self.salt = np.dstack([self.salt[:,None]-northVariance[None,:]]\
                         *self.xdim)
        return self.salt
    
    def zSal(self):
        alpha     = 0.0002 ;sBeta = 7.4e-4; sigma_0 = 999.8
        zMax  = self.depth[zdim-1]; zMin = self.depth[0]
        sigma = ((np.log(self.depth)-np.log(zMax))*(sigmaMax-sigmaMin)\
                / np.log(zMax/zMin)) + sigmaMax 
        self.salt  = sigma/(sBeta*sigma_0)
        return self.salt    
    
    def obc(self):
        obcts = self.salt[:,0,:]
        obctn = self.salt[:,-1,:]
        return obctn, obcts

class Velocities(Grid):
    
    def __init__(self, writeCase):
        Grid.__init__(self)
        self.writeCase = writeCase

    def extract_velocities(self, case):
        self.readPath = config.readPath(case)
        #vels_file = Dataset(self.readPath + 'vels.nc', mode='r')
        uvel = vels_file.variables['UVEL'][-1].T
        vvel = vels_file.variables['VVEL'][-1].T
        self.writeBin(uvel,self.writeCase + '_uvel.bin')
        self.writeBin(vvel,self.writeCase + '_vvel.bin')

class State(Grid):

    def __init__(self, ini_params):
        Grid.__init__(self, ini_params)
        self.ini_params = ini_params
        self.shape = (self.ini_params['zdim'],
                       self.ini_params['ydim'],
                       self.ini_params['xdim'])
    
    def ini_field_linear_grad(self, dimension, dimension_mesh, low, high):
        grad = (high - low) / (dimension -1)
        array = (dimension_mesh * grad) + low        
        print ('dimension', dimension)
        print ('line of array', array.shape)
        return array

    def ini_field_hill(self, dimension, dimension_mesh, amp):
        array = amp * np.sin(dimension_mesh * np.pi / dimension)
        return array
  
    def ini_shice(self, iceDep, t, s, hFacC): 
        ''' 
        initialise ice shelf's initial topography and pressure load along
        base 
        ''' 

        # Remove Bathy
        z = int(self.ini_params['zdim'])
        y = int(self.ini_params['ydim'])
        x = int(self.ini_params['xdim'])
        hFacC_new = np.ones((z,y,x))
        for i in range(x):
            for j in range(y):
                for k in range(z):
                    if hFacC[k,j,i] < 1:
                        hFacC_new[k,j,i] = hFacC[k,j,i]
                    else:
                        break
        print ('hgac', hFacC_new[1,0,0])                
        print ('hgac', hFacC_new[0,0,3])                
        #shice_topo = - np.ones((self.xdim, self.ydim)) * iceDep 
        shice_topo = iceDep 
        tRef = -1.9
        sRef = 34.4
        
        # Calcualte density
        rho = self.ini_params['rhoConst'] * \
                            (- self.ini_params['tAlpha'] * (t - tRef) 
                             + self.ini_params['sBeta']  * (s - sRef))

        # Depth integrate to ice base
        shice_ini_p = np.sum((1-hFacC_new) * self.ini_params['g'] * rho, axis=0)

        return shice_topo, shice_ini_p


    def add_heat_blob(self, x, y, x0, y0, r):
        amp = 0.001
        sigma = 1000
        array = amp * np.exp(-((x-x0)**2)/sigma - ((y-y0)**2)/sigma)
        return array

    def plot_shice():
        pass

    def ini_p_load(self, dimension, dimension_mesh, low, high):
        grad = (high - low) / (dimension -1)
        array = (dimension_mesh[0] * grad) + low
        return array

    def ini_p_force(self, dimension, dimension_mesh, low, high):
        resolution = 1
        rhoConst = 1030
        grad = (high - low) / ((dimension -1) * resolution * rhoConst)
        array = np.full(dimension_mesh.shape, grad)
        return array
     
       
bathy = 0
wind  = 0
vels  = 0 
# MITgcm binary file is saved in (x,y) format

if vels:
    name = 'Binary/SPBC_303_uvel.bin'
    vels = Velocities('SPBC_303')
    vels.extract_velocities('SPBC_286')
    y = bathy.readBin(name,xdim,ydim)
    bathy.plot_single(y)
    
if bathy:
    name = case_name + '_bathy.bin'
    print (name)
    bathy = Bathymetry()
    bathy.shift_coordinates()
    #bathy.sheer_coordinates()
    bathy.sawtooth_ridge(peak_wall=0)
    #bathy.meridional_wall()
    #bathy.boundary_slope()
    #bathy.bathy[210:220,192:288]  = 0.0
    #bathy.gap()
    #bathy.smooth_gap()
    #bathy.chop_ridge_end()
    ###bathy.plateau_slope()
    bathy.north_south_boundary()
    #bathy.east_west_boundary()
    #bathy.smooth_gap(orientation=-1)
    b = bathy.get_bathy()#[::-1]
    #b1 = b[xdim/2:]
    #b2 = b[:xdim/2]
    #b = np.concatenate((b1,b2))
    bathy.writeBin(b,name)
    y = bathy.readBin(name,int(xdim),int(ydim))
    bathy.case = '001'
    fig = plt.figure(0)
    p = plt.pcolormesh(y)
    plt.axis('equal')
    plt.colorbar(p)
    #plt.figure(1)
    #plt.plot(y[:,int(ydim/2.0)])
#    plt.show()
    

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

plt.show()
