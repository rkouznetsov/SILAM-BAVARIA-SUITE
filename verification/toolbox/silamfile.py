"""
A module for handling SILAM output in grads and netcdf formats. The functions are as follows:
- the grid and vertical definition is parsed from the super-ctl/netcdf file
- a projection is defined using the grid information
- the variable descriptions are parsed.
"""

import numpy as np
import projections
from toolbox import namelist, gradsfile, gridtools, verticals, ncreader
from support import pupynere as netcdf
try:
    import netCDF4 as nc4
    from toolbox import nc4reader
except:
    nc4 = None
    
from os import path

class SilamfileError(Exception):
    pass

class SilamLatLonGrid(gridtools.Grid):
    """A gridtools.Grid subclass defining the SILAM latlon grid with
    interfaces to the namelist objects"""
    def __init__(self, *args):
        gridtools.Grid.__init__(self, *args)
        # print self.proj.south_pole
        try:
            self.proj.south_pole
            self.grid_type = 'LON_LAT'
        except AttributeError:
            raise ValueError('The projection does not appear to be lat/lon')        
                        
    def as_namelist(self):
        hash = {'grid_type' : self.grid_type,
                'grid_title' : '',
                'lon_s_pole' : self.proj.south_pole[0],
                'lat_s_pole' : self.proj.south_pole[1],
                'nx' : self.nx,
                'ny' : self.ny,
                'dx' : self.dx,
                'dy' : self.dy,
                'resol_flag' : 128,
                'ifReduced' : 0,
                'wind_component' : 0,
                'reduced_nbr_str' : 0,
                'lat_pole_stretch' : 0.0,
                'lon_pole_stretch' : 0.0,
                'lon_start' : self.x0,
                'lat_start' : self.y0
                }
        nl = namelist.Namelist('grid')
        for key, val in hash.iteritems():
            nl.put(key, val)
        return nl

    @staticmethod
    def from_namelist(nl):
        dx, dy, lon_start, lat_start = (float(nl.get_uniq(key))
                                        for key in ('dx', 'dy', 'lon_start', 'lat_start'))
        try:
            nx, ny = int(nl.get_uniq('nx')), int(nl.get_uniq('ny'))
        except KeyError:
            lon_end, lat_end = float(nl.get_uniq('lon_end')), float(nl.get_uniq('lat_end'))
            nx, ny = (lon_end-lon_start) / dx + 1, (lat_end-lat_start) / dy + 1
            nx, ny = int(nx), int(ny)
            
        sp_lon, sp_lat = float(nl.get_uniq('lon_s_pole')), float(nl.get_uniq('lat_s_pole'))
        proj = projections.LatLon(sp_lat, sp_lon)
        grid = SilamLatLonGrid(lon_start, dx, nx, lat_start, dy, ny, proj)
        return grid

LatLonGrid = SilamLatLonGrid
    
class SilamVertical:
    """Definitions for the verticals appearing in SILAM output. So far
    only the simple cases are considered - the vertical is assumed to be
    described by a 1d list of numbers."""
    def number_of_levels(self):
        return len(self.values)

    def as_namelist(self):
        nl = namelist.Namelist('vertical')
        if self.number_of_levels() > 0:
            nl.put('number_of_levels', self.number_of_levels())
            nl.put(self.__class__.values_label, ' '.join('%g' % x for x in self.values))
        nl.put('vertical_method', self.__class__.vertical_method.upper())
        if self.__class__.level_type:
            nl.put('level_type', self.__class__.level_type.upper())
        return nl

    @classmethod
    def from_namelist(cls, nl):
        # This function defines the vertical from a namelist. The
        # attributes taken depend on the subclass.
        levtype, method = nl.get_uniq('level_type'), nl.get_uniq('vertical_method')
        values = [float(x) for x in nl.get_uniq(cls.values_label).split()]
        return cls(values)

    def num_levels(self):
        return len(self.values)

    def unit(self):
        return self.__class__.unit

    def values(self):
        return self.values
    
class SilamHybrid: # Nothing to inherit from SilamVertical
    # To create a vertical class, one defines these three class
    # attributes. They are used for the common operations defined in
    # the SilamVertical class.
    vertical_method = 'custom_layers'
    level_type = 'hybrid'
    values_label = 'hybrid_coefficients_bottom'
    values_top_label = 'hybrid_coefficients_domain_top'
    unit = 'Pa'
    
    def __init__(self, a_half, b_half):
        self.a_half = tuple(a_half)
        self.b_half = tuple(b_half)

    def number_of_levels(self):
        return len(a_half) - 1  #Level boundaries stored

    def as_namelist(self):
        nlevs = self.number_of_levels()
        nl = namelist.Namelist('vertical')

        if nlevs > 0:
            nl.put('number_of_levels', nlevs)
            for i in range(nlevs):
                  nl.put(self.__class__.values_label, '%d %g %g' %(i,a_half[i], b_half[i]))
            nl.put(self.__class__.values_top_label, '%d %g %g' %(i,a_half[nlevs+1], b_half[nlevs+1]))

        nl.put('vertical_method', self.__class__.vertical_method.upper())
        if self.__class__.level_type:
            nl.put('level_type', self.__class__.level_type.upper())
        return nl

    @classmethod
    def from_namelist(cls, nl):
        # This function defines the vertical from a namelist. The
        # attributes taken depend on the subclass.
        levtype, method = nl.get_uniq('level_type'), nl.get_uniq('vertical_method')
        if levtype.lower() == "hybrid":
            print "Hybrid vertical"
        else:
            print "SilamHybrid initializing with not 'hybrid' level type"
            raise SilamfileError
            return None

        a_half = []
        b_half = []
        for l in nl.get(cls.values_label):
                n, a, b = l.split()
                a_half.append(a)
                b_half.append(b)
        l=nl.get_uniq(cls.values_top_label)
        a, b = l.split()
        a_half.append(a)
        b_half.append(b)
        return SilamHybrid(a_half,b_half)


    def num_levels(self):
        return len(a_half) - 1 

    def unit(self):
        return self.__class__.unit

    def values(self):
        return tuple([ (a_half[i] + a_half[i+1])*0.5 + (b_half[i] + b_half[i+1])*0.5*101325
                       for i in range(self.number_of_levels())])

    def thickness(self):
        return tuple([ (a_half[i] - a_half[i+1]) + (b_half[i] - b_half[i+1])*101325
                       for i in range(self.number_of_levels())])

    def boundaries(self):
        return tuple([ _half[i] + b_half[i]*101325 for i in range(self.number_of_levels()+1)])

    def midpoints(self):
        return self.values



class SilamHeightLevels(SilamVertical):
    # To create a vertical class, one defines these three class
    # attributes. They are used for the common operations defined in
    # the SilamVertical class.
    vertical_method = 'custom_levels'
    level_type = 'height_from_surface'
    values_label = 'levels'
    unit = 'm'
    
    def __init__(self, midpoints):
        self._midpoints = tuple(midpoints)
        self.values = self._midpoints

    def thickness(self):
        return verticals.thickness(self._midpoints)

    def boundaries(self):
        return verticals.boundaries(self.thickness())

    def midpoints(self):
        return self.values
    
class SilamHeightLayers(SilamVertical):
    vertical_method = 'custom_layers'
    level_type = 'height_from_surface'
    values_label = 'layer_thickness'
    unit = 'm'
    
    def __init__(self, thickness):
        self._thickness = tuple(thickness)
        self.values = self._thickness
        
    def boundaries(self):
        return verticals.boundaries(self._thickness)
        #raise SilamfileError('not implemented')

    def midpoints(self):
        return verticals.midpoint(self._thickness)

    def thickness(self):
        return self._thickness
      
class SilamSurfaceVertical(SilamVertical):
    vertical_method = 'surface_level'
    level_type = None
    values_label = None
    unit = ''
    def __init__(self, *args):
        self.values = 0.0,
    def thickness(self):
        raise ValueError()
    def boundaries(self):
        raise ValueError()
    @staticmethod
    def from_namelist(nl):
        return SilamSurfaceVertical()
    
def get_vertical_class(vertical_method, level_type):
    # Return the correct class for vertical definition based on the
    # two values (typically from a namelist).
    classes = SilamHeightLevels, SilamHeightLayers, SilamSurfaceVertical, SilamHybrid
    for cls in classes:
        if cls.vertical_method == vertical_method.lower() and (cls.level_type == level_type.lower()
                                                               or not cls.level_type):
            return cls
    raise SilamfileError('Unsupported vertical: %s %s' % (vertical_method, level_type))

class ProjectedReader(gradsfile.GriddedDataReader, gradsfile.GenericWrapper):
    """A GriddedDataReader with a nontrivial projection, defined by the
    gridtools.Grid instance. Can wrap any GriddedDataReader object,
    and overloads only the meshgrid and indices methods. These methods
    will now be based on the geographic coordinates, while x() and y()
    will return the grid-native coordinates. This serves the two
    common uses: plotting maps and extracting timeseries."""
    def __init__(self, reader, grid):
        methods = 'coordinates', 'indices', '_get_geo_coords', 'meshgrid'
        gradsfile.GenericWrapper.__init__(self, reader, methods)
        self._orig_descr = reader.describe()
        self.grid = grid
        self._get_geo_coords()
        
    def _get_geo_coords(self):
        # midpoints
        X, Y = gradsfile.GriddedDataReader.coordinates(self)
        shape = X.shape
        lon, lat = self.grid.proj_to_geo(X.ravel(), Y.ravel())
        self.lon = lon.reshape(shape)
        self.lat = lat.reshape(shape)
        # vertices
        X, Y = gradsfile.GriddedDataReader.coordinates(self, True)
        shape = X.shape
        lon, lat = self.grid.proj_to_geo(X.ravel(), Y.ravel())
        self.lon_vx = lon.reshape(shape)
        self.lat_vx = lat.reshape(shape)

    def coordinates(self, vertices=False):
        if vertices:
            return self.lon_vx, self.lat_vx
        else:
            return self.lon, self.lat
        
    def meshgrid(self, vertices=False):
        warnings.warn('meshgrid() should be replaced by coordinates()')
        if vertices:
            return self.lon_vx.T, self.lat_vx.T
        else:
            return self.lon.T, self.lat.T
        
    def indices(self, lon, lat):
        # Project lon, lat into the native coordinates and return the indices
        x, y = self.grid.geo_to_grid(lon, lat)
        #return int(x+0.5), int(y+0.5)
        # better when x, y can be negative:
        return int(round(x)), int(round(y))

    def describe(self):
        description = 'PROJ(%s)' % self._orig_descr
        return description


class BaseSilamfile:
    def _check_match(self, metadata, key, request, get_value):
        # the attribute matching rules are here:
        match = True
        try:
            request.__call__
            request_callable = True
        except AttributeError:
            request_callable = False
        testval = get_value(metadata, key).split()[0]
        if request_callable:
            try:
                match = match and request(testval)
            except ValueError:
                return False
        else:
            reqval = request
            valtype = reqval.__class__
            try:
                testval = valtype(testval)
            except ValueError:
                return False
            if isinstance(testval, float):
                match = match and np.abs(testval-reqval) < 1e-6
            else:
                match = match and testval == reqval

        return match

    def get_variable(self, **conditions):
        """ A convenience version of get_variables: return a single variable, and if
        more than one matches, raise an error. """
        
        selected = self.get_variables(**conditions)
        if len(selected) > 1:
            raise SilamfileError('More than one variable matches')
        return selected[0]

    
class Silamfile(BaseSilamfile):
    """A class for parsing and describing the information in a silam super
    ctl file: the grid, vertical and variable descriptions. Does not read
    the data, that is handled by requesting a reader object with a given
    grads variable expression.

    Arguments:
    superctl : the superctl file path
    relocate : modify the path of the ctl file. If True, assume that
               the ctl is in the same directory as the superctl. If a string,
               use the given directory. Also the binary files will be relocated
               accordingly."""
    
    def __init__(self, superctl, relocate=None):
        nlgrp = namelist.NamelistGroup.fromfile(superctl)
        nl_general = nlgrp.get('general')
        grid = SilamLatLonGrid.from_namelist(nl_general)
        vertical_class = {('custom_layers', 'height_from_surface') : SilamHeightLayers,
                          ('custom_levels', 'height_from_surface') : SilamHeightLevels}
        if not nl_general.has_key('level_type'):
            vertical = None
        else:
            vert_method, vert_type = (nl_general.get_uniq('vertical_method'),
                                      nl_general.get_uniq('level_type'))
            vertical_class = get_vertical_class(vert_method, vert_type)
            vertical = vertical_class.from_namelist(nl_general)

        ctlpath_original = nl_general.get_uniq('ctl_file_name')
        if ctlpath_original.startswith('^'):
            ctlpath_original = path.join(path.dirname(superctl), ctlpath_original[1:])
        if relocate is True:
            ctlpath = path.join(path.dirname(superctl), path.basename(ctlpath_original))
            self.gradsdescriptor = gradsfile.GradsDescriptor(ctlpath)
            self.gradsdescriptor.relocate()
        elif relocate:
            ctlpath = path.join(relocate, path.basename(ctlpath_original))
            self.gradsdescriptor = gradsfile.GradsDescriptor(ctlpath)
            self.gradsdescriptor.relocate(relocate)
        else:
            ctlpath = ctlpath_original
            self.gradsdescriptor = gradsfile.GradsDescriptor(ctlpath)
            
        self.grid = grid
        self.vertical = vertical
        self.nlgrp = nlgrp
        self.ctlpath = ctlpath
        
    def get_ctl(self):
        return self.ctlpath
        
    def get_variables(self, **conditions):
        """Find variables from the superctl by matching a set of values given
        as keyword arguments. If a given variable does not have the specified
        attribute, it is never selected.
        
        The value read from the namelist is coerced to the type of the
        respective argument, if possible. This means that if the
        argument is a string, the values are compared as strings, same
        for floats and integers. If the comparison is impossible, the
        variable is ignored. Only the first word of the value is taken
        - units are ignored.
        
        Examples:
        
        # Find all PM variables (depositions, concentration, all modes)
        variables = silamf.get_variables(substance_name='PM')
        # Find all PM concentration with a mean diameter of 0.5 micron
        variables = silamf.get_variables(substance_name='PM', fix_diam_mode_mean_diameter=0.5)
        # This will probably not work
        variables = silamf.get_variables(substance_name='PM', fix_diam_mode_mean_diameter='0.5')
        
        """
        selected = []
        for nl in self.nlgrp.lists.values():
            if nl.name == 'general':
                continue
            match = True
            for key, request in conditions.items():
                if not nl.has_key(key):
                    match = False
                    break
                match = match and self._check_match(nl, key, request, get_value=lambda m, k: m.get_uniq(key))
            if match:
                selected.append(nl.name)
        return selected
    

    def get_reader(self, grads_expression, ind_level=None, projected=True, mask_mode='first',
                   grads_class=gradsfile.MMGradsfile):
        """Return a ProjectedReader reading the given grads variable expression on a given level."""
        grfile = grads_class(self.gradsdescriptor, grads_expression, ind_level, mask_mode)
        if projected:
            return ProjectedReader(grfile, self.grid)
        else:
            return grfile

    def get_var_namelist(self, var):
        for nl in self.nlgrp.lists.values():
            if nl.name == var:
                return nl
        raise SilamfileError('Variable not found: %s' % var)

    def get_general_param(self, key):
        nl_general = self.nlgrp.lists['general']
        return nl_general.get(key)

    def get_attribute(self, variable, key):
        nl = self.get_var_namelist(variable)
        value = nl.get_uniq(key)
        return value
        
    
class SilamNCFile(BaseSilamfile):
    """
    Class for parsing variables from silam netcdf output. Netcdf-4 requires the NetCDF4
    python library and the nc4reader toolbox module. Variables can be selected by matching
    the attributes. For more info see Silamfile.
    """
    
    def __init__(self, netcdf_file=None, descriptor_file=None, ncversion=4):
        """
        Create a parser object - either from a netcdf file or a .ctl file. In the latter
        case, the first file pointed by the ctl will be examined, and if a reader is
        requested, it will include the whole dataset.
        """
        if not netcdf_file and not descriptor_file:
            raise ValueError('Not enough arguments given')
        if netcdf_file and descriptor_file:
            raise ValueError('Only one of netcdf_file and descriptor_file may be defined')
        try:
            if ncversion == 4:
                nc_class = nc4.Dataset
            elif ncversion == 3:
                nc_class = netcdf.netcdf_file
            else:
                raise ValueError('Invalid ncversion')
        except AttributeError:
            raise ValueError('Netcdf version %i not available' % ncversion)

        if descriptor_file:
            descr = ncreader.NCDescriptor.fromfile(descriptor_file)
            netcdf_file = descr.get_file(descr.times()[0])
            self.descr_path = descriptor_file
            self.reader_class = {3:ncreader.NCDataset, 4:nc4 and nc4reader.NC4Dataset}[ncversion]
        else:
            self.netcdf_path = netcdf_file
            self.descr_path = None
            self.reader_class = {3:ncreader.NCReader, 4:nc4 and nc4reader.NC4Reader}[ncversion]
            self.reader_class_expr = {3:ncreader.NCExpression, 4:nc4 and nc4reader.NC4Expression}[ncversion]
        try:
            nc_file = nc_class(netcdf_file, 'r')
        except:
            print 'Failed to open:', netcdf_file
            raise

        self._parse_grid_vert(nc_file)
        self._parse_attributes(nc_file)
        nc_file.close()
        
    def _parse_grid_vert(self, nc_file):
        rotated = False
        if 'lon' in nc_file.dimensions and 'lat' in nc_file.dimensions:
            lonname='lon'
            latname='lat'
        elif 'longitude' in nc_file.dimensions and 'latitude' in nc_file.dimensions:
            lonname='longitude'
            latname='latitude'
        elif 'rlon' in nc_file.dimensions and 'rlat' in nc_file.dimensions:
            lonname='rlon'
            latname='rlat'
            rotated = True
        else:
            raise ValueError('Cannot determine grid')

        try:
          if nc_file.grid_projection != 'lonlat':
            raise ValueError('Projection not supported: %s' % nc_file.grid_projection)
        except AttributeError:
            print "File has no grid_projection attribute, assuming lonlat and hope for the best"
            pass
            
        if rotated:
            try:
                sp_lon, sp_lat = nc_file.pole_lon, nc_file.pole_lat
            except AttributeError:
                rpvar=nc_file.variables['rp']
                assert rpvar.grid_mapping_name == "rotated_latitude_longitude"
                sp_lon = (rpvar.grid_north_pole_longitude + 360.) % 360.  - 180.
                sp_lat = - rpvar.grid_north_pole_latitude
            proj = projections.LatLon(sp_lat, sp_lon)
        else:
            proj = projections.LatLon()
        lon, lat = nc_file.variables[lonname][:], nc_file.variables[latname][:]
        nx, ny = len(lon), len(lat)
        # Try to recover the accuracy
        dx, dy = (lon[-1]-lon[0])/(nx-1), (lat[-1]-lat[0])/(ny-1)
        dx = 1e-5*round(dx*1e5)
        dy = 1e-5*round(dy*1e5)
        self.grid = SilamLatLonGrid(lon[0], dx, nx, lat[0], dy, ny, proj)

        if 'height' in nc_file.dimensions:
           hgt = nc_file.variables['height'][:]
           self.vertical = SilamHeightLevels(hgt)
        elif 'hybrid' in nc_file.dimensions:
           try:
               self.vertical = SilamHybrid(nc_file.variables['a_half'][:],nc_file.variables['b_half'][:])
           except Exception as e:
               print str(e)
               print "Failed to parse hybrid vertical pretending it is height"
               # Pretend that it is still height -- does not matter for the stations
               hgt = nc_file.variables['hybrid'][:]
               self.vertical = SilamHeightLevels(hgt)
               pass
        else:
           self.vertical=None

    _GLOB = 'glob'
    def _parse_attributes(self, nc_file):
        attr = {}
        for var_name, ncvar in nc_file.variables.items() + [(self._GLOB, nc_file)]:
            if not var_name in attr:
                attr[var_name] = {}
            for fldname in dir(ncvar):
                if fldname.startswith('__'):
                    continue
                fldval = getattr(ncvar, fldname)
                # take all float, int or character attributes.
                if not (isinstance(fldval, basestring) or isinstance(fldval, float) or isinstance(fldval, int)):
                    continue
                attr[var_name][fldname] = fldval
        self._attrs = attr
        
    def get_variables(self, **conditions):
        selected = []
        for var_name, var_attrs in self._attrs.items():
            match = True
            for key, request in conditions.items():
                if not key in var_attrs:
                    match = False
                    break
                match = match and self._check_match(var_attrs, key, request, get_value=lambda a, k: a[k])
            if match:
                selected.append(var_name)
        return selected

    def get_reader(self, expression, ind_level=None, mask_mode='first'):
        projected = True # Projected reader used always. 
        if self.descr_path is not None:
            reader = self.reader_class(self.descr_path, expression, ind_level, mask_mode=mask_mode)
        elif ncreader.NCExpression.is_expression(expression):
            reader = self.reader_class_expr(self.netcdf_path, expression, ind_level, mask_mode=mask_mode)
        else:
            reader = self.reader_class(self.netcdf_path, expression, ind_level, mask_mode=mask_mode)
        if projected:
            return ProjectedReader(reader, self.grid)
        else:
            return reader

    def get_general_param(self, key):
        return self._attrs[self._GLOB][key]

    def get_attribute(self, variable, key):
        return self._attrs[variable][key]
        
if __name__ == '__main__':
    superctl = '/home/vira/silam/output/blh/blh_ALL_SRCS_20101001.grads.super_ctl'
    superctl = '/home/vira/work/macc_ii/r-eda-env/wrk/test_cb4.02/iteration/gradient_96.grads.super_ctl'
    silamf = Silamfile(superctl)
    variables = silamf.get_variables(substance_name='O3')
    print variables
    #expr = '+'.join(variables[:-1])
    expr = '2*%s' % variables[0]
    reader = silamf.get_reader(expr, 0)
