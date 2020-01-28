import numpy as np
from numpy import cos, sin, arcsin, arctan2, sqrt

R_earth = 6378.1e3 
deg2rad = np.pi / 180. 

def gc_distance(x1,x2,y1,y2, r=6378.0):
    xr1 = x1*np.pi/180
    xr2 = x2*np.pi/180
    yr1 = y1*np.pi/180
    yr2 = y2*np.pi/180
    
    dx = np.abs(xr1-xr2)
    dsg = (sin((yr2-yr1)/2))**2 + cos(yr1)*cos(yr2)*(sin(dx/2))**2
    dsg = 2*arcsin(sqrt(dsg))
    return r*dsg

def dd(a, dim=0, dx = 1.0):
    # derivative along dimension dim
    d = np.zeros(a.shape)    
    if dim > 2:
        raise ValueError('Not so many dims!')
    if dim == 0:
        d[1:-1,...] = (a[2:,...]-a[0:-2,...]) / (2*dx)
        d[0,...] = (a[1,...]-a[0,...]) / dx
        d[-1,...] = (a[-1,...]-a[-2,...]) / dx
    if dim == 1:
        d[:,1:-1,...] = (a[:,2:,...]-a[:,:-2,...]) / (2*dx)
        d[:,0,...] = (a[:,1,...]-a[:,0,...]) / dx
        d[:,-1,...] = (a[:,-1,...]-a[:,-2,...]) / dx
    if dim == 2:
        d[:,:,1:-1,...] = (a[:,:,2:,...]-a[:,:,:-2,...]) / (2*dx)
        d[:,:,0,...] = (a[:,:,1,...]-a[:,:,0,...]) / dx
        d[:,:,-1,...] = (a[:,:,-1,...]-a[:,:,-2,...]) / dx
    return d

def divergence(u,v,w=None, d=None):
    if d == None:
        d = np.ones(u.ndim)
    else:
        d = np.array(d)
    # sanity check
    if w != None:
        assert(u.ndim == 3)
        fields = (u,v,w)
    else:
        assert(u.ndim == 2 and v.ndim == 2)
        fields = (u,v)

    div = np.zeros(u.shape)
    for i in range(0,len(fields)):
        div += dd(fields[i], i, d[i])
    return div

def grad(a, d = None):
    if d == None:
        d = np.ones(a.ndim)
    else:
        d = np.array(d)
    g = [dd(a,i,d[i]) for i in range(0,a.ndim)]
    return g

def curl(u, v, w=None, d=None):
    if u.ndim > 2:
        raise ValueError('Not implemented')

    if w is not None:
        raise ValueError('Not implemented')
    else:
        c = dd(v, 0, d[0]) - dd(u, 1, d[1])

    return c

def curl2d(psi, dx, dy):
    return dd(psi, 1, dy), -dd(psi, 0, dx)


def aave(F, lats, lons):
    # Area weighted average of F
    # lats, lons are assumed to represent a regular grid.
    
    dx = lons[1] - lons[0]
    dy = lats[1] - lats[0]

    if not all(lats[1:] - lats[0:-1] == dy) \
       or not all(lons[1:] - lons[0:-1] == dx):
        raise ValueError('Sorry, the grid has to be regular.')
    
    area = R**2 * dy*deg2rad * dx*deg2rad * np.cos(lats*deg2rad)
    area = np.outer(np.ones(lons.size), area)

    area_total = R**2 * deg2rad * (lons[-1] - lons[0])\
        * (sin(deg2rad*lats[1]) - sin(deg2rad*lats[0]))
    
    return np.sum(F * area) / area_total
    

def area(lats, lons):
    # Area of gridpoints defined by lats & lons.
    lats = np.asarray(lats)
    lons = np.asarray(lons)
    dx = abs(lons[1] - lons[0])
    dy = abs(lats[1] - lats[0])
    area = R_earth**2 * dy*deg2rad * dx*deg2rad * np.cos(lats*deg2rad)
    area = np.outer(area, np.ones(lons.size))
    return area.T

#
# Spherical geometry & rotations
#

def get_basis(latdeg, londeg):
    lat = latdeg * deg2rad
    lon = londeg * deg2rad
    
    B = np.empty((3,2))
    
    B[0,0] = -sin(lon)
    B[1,0] = cos(lon)
    B[2,0] = 0.
    B[0,1] = -cos(lon)*sin(lat)
    B[1,1] = -sin(lon)*sin(lat)
    B[2,1] = cos(lat)

    return B

def get_rotation(lat_pole, lon_pole):
    # Right multiplication -> world to rotated
    # Left multiplication -> rotated to world 
    pi = np.pi
    
    latr = deg2rad * lat_pole
    lonr = deg2rad * lon_pole

    alpha = (pi/2. + lonr)
    beta = pi/2. - latr
    gamma = 0.0#-alpha

    ca = cos(alpha)
    cb = cos(beta)
    cg = cos(gamma)
    sa = sin(alpha)
    sb = sin(beta)
    sg = sin(gamma)

    R = np.empty((3,3))

    R[0,0] = ca*cg - sa*cb*sg
    R[1,0] = sa*cg + ca*cb*sg
    R[2,0] = sb*sg
    R[0,1] = -ca*sg - sa*cb*cg
    R[1,1] = -sa*sg + ca*cb*cg
    R[2,1] = sb*cg
    R[0,2] = sb*sa
    R[1,2] = -sb*ca
    R[2,2] = cb

    return R

def get_rotation_v2(lat_pole, lon_pole, gamma=0.0):
    # Right multiplication -> world to rotated
    # Left multiplication -> rotated to world 
    pi = np.pi
    
    latr = deg2rad * lat_pole
    lonr = deg2rad * lon_pole

    alpha = lonr
    beta = pi/2 + latr
    gamma = deg2rad*gamma
    
    ca = cos(alpha)
    cb = cos(beta)
    cg = cos(gamma)
    sa = sin(alpha)
    sb = sin(beta)
    sg = sin(gamma)

    R = np.empty((3,3))

    R[0,0] = ca*cb*cg - sa*sg
    R[1,0] = sa*cb*cg + ca*sg
    R[2,0] = sb*cg
    R[0,1] = -ca*cb*sg - sa*cg
    R[1,1] = -sa*cb*sg + ca*cg
    R[2,1] = -sb*sg
    R[0,2] = -ca*sb
    R[1,2] = -sa*sb
    R[2,2] = cb

    return R
    

class Grid:
    def __init__(self, x0, dx, nx, y0, dy, ny, projection):
        self.x0 = x0
        self.y0 = y0
        self.dx = dx
        self.dy = dy
        self.nx = nx
        self.ny = ny
        self.proj = projection
        self.shape = (nx, ny)

    def grid_to_geo(self, X, Y):
        x = self.x0 + X*self.dx
        y = self.y0 + Y*self.dy
        lon, lat = self.proj.to_geo(x,y)
        return lon, lat
    
    def geo_to_grid(self, lon, lat):
        x, y = self.proj.to_proj(lon, lat)
        X = (x-self.x0) / self.dx
        Y = (y-self.y0) / self.dy
        return X, Y

    def proj_to_grid(self, x, y):
        X = (x-self.x0) / self.dx
        Y = (y-self.y0) / self.dy
        return X, Y

    def proj_to_geo(self, x, y):
        return self.proj.to_geo(x, y)

    def x(self):
        return self.x0 + np.arange(self.nx)*self.dx
    
    def y(self):
        return self.y0 + np.arange(self.ny)*self.dy
    
    def cell_dimensions(self):
        Y, X = np.meshgrid(np.arange(self.ny), np.arange(self.nx))
        lon, lat = self.grid_to_geo(X.ravel(), Y.ravel())
        dsx = self.proj.ds_dxi(lat, lon)*self.dx
        dsy = self.proj.ds_dnu(lat, lon)*self.dy
        return dsx, dsy

def gems_grid():
    import projections
    return Grid(-17.0, 0.2, 265, 33.0, 0.2, 195, projections.LatLon())

def gems_grid_ext():
    import projections
    return Grid(-25.0, 0.2, 351, 30.0, 0.2, 220, projections.LatLon())

def latLonFromPoints(lats, lons):
    import projections
    dx = lons[1]-lons[0]
    dy = lats[1]-lats[0]
    (nx, ny) = (len(lons), len(lats))
    return Grid(lons[0], dx, nx, lats[0], dy, ny, projections.LatLon())

class Rotation:
    def __init__(self, lat_pole, lon_pole, gamma=0.0):
        # SOUTH POLE
        self.R = get_rotation_v2(lat_pole, lon_pole, gamma)
        self.lat_pole = lat_pole * deg2rad
        self.lon_pole = lon_pole * deg2rad

    def _transform_latlon(self, latdeg, londeg, R):
        lat = latdeg * deg2rad
        lon = londeg * deg2rad
        if np.isscalar(latdeg):
            shape = ()
        else:
            shape = latdeg.shape

        
        r1 = np.empty((3,) + shape)
        r1[0,...] = cos(lon)*cos(lat)
        r1[1,...] = sin(lon)*cos(lat)
        r1[2,...] = sin(lat)

#        print "-------------"
#        print "R=", R 
#        print R.shape
#        print "r1=", r1
#        print r1.shape
        r2 = np.tensordot(R, r1, axes=(1,0))

        ###Numerics causes problems in poles (very rarely)...
        #Dirty hack to clean it
        nlatd = arcsin(r2[2,...]*0.9999999) / deg2rad

        nlond = arctan2(r2[1,...], r2[0,...]) / deg2rad
        
        return nlatd, nlond
 
    def vector3ToLatLon(self, r):
        nlatd = arcsin(r[2]) / deg2rad
        nlond = arctan2(r[1], r[0]) / deg2rad
        return nlatd, nlond

    def latLonToVector3(self, latdeg, londeg):
        lat = latdeg * deg2rad
        lon = londeg * deg2rad
        
        r1 = np.empty(3)
        r1[0] = cos(lon)*cos(lat)
        r1[1] = sin(lon)*cos(lat)
        r1[2] = sin(lat)

        return r1

    def latLonToWorld(self, latdeg, londeg):
        return self._transform_latlon(latdeg, londeg, self.R)

    def latLonToRotated(self, latdeg, londeg):
        return self._transform_latlon(latdeg, londeg, self.R.T)

    def vector3ToWorld(self, vector3):
        return np.dot(self.R, vector3)
    
    def vector3ToRotated(self, vector3):
        return np.dot(self.R.T, vector3)

    def transform_wind(self, u, v, latdeg, londeg, R):
        # This is from world to rotated!!
        # latdeg, londeg *already* in rotated crd
        lat = latdeg * deg2rad
        lon = londeg * deg2rad
        
        lat2, lon2 = self.transform_latlon(latdeg, londeg, R)
        # Basis in the reference frame
        B1 = get_basis(lat2, lon2)
        # Rotate B1 into modified frame
        B1 = dot(R.T, B1)
        # Basis in the rotated frame
        B2 = get_basis(latdeg, londeg)
        # Project the vector in modified B1-basis into B2 basis
        r = dot(B2.T, B1)
        V = np.array([u,v])
        return dot(r[0,:], V), dot(r[1,:], V)

    
    def windToRotated(self, u, v, latdeg, londeg):
        return transform_wind(u, v, latdeg, londeg, self.R)
    
    def windToWorld(self, u, v, latdeg, londeg):
        return transform_wind(u, v, latdeg, londeg, self.R.T)
