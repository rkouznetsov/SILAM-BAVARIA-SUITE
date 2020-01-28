import numpy as np
import datetime as dt

from toolbox import verticals, gradsfile, util, gridtools
try:
    import figs
except ImportError:
    print 'Importing figs failed'
try:
    from mpl_toolkits import basemap
except ImportError:
    print 'Importing basemap failed'
#import figs

try:
    from toolbox import interp2d
    have_interp2d = True
except ImportError:
    print 'Importing interp2d failed, using basemap'
    have_interp2d = False
    
import datetime as dt

class GenericSectionReader(gradsfile.GenericWrapper):
    """
    Calculate cross sections from gridded 3D data, with interpolation to arbitrary points.
    """
    def __init__(self, reader, points_x, points_y):
        methods = 'x', 'y', 'meshgrid', 'ndims', 'nvars', 'indices', 'read'
        gradsfile.GenericWrapper.__init__(self, reader, methods)
        self._xout = points_x
        self._yout = points_y
        self.reader = reader
        
    def read(self, n=1):
        nz = len(self.reader.z())
        nx = len(self.reader.x())
        ny = len(self.reader.y())
        data3d = self.reader.read(n, squeeze=False).reshape((nx,ny,nz))
        return figs.make_cross_section(data3d,
                                       self.reader.x(), self.reader.y(), self.reader.z(),
                                       xout=self._xout, yout=self._yout)
    def x(self):
        return self._xout

    def y(self):
        return self._yout


class SectionReader(GenericSectionReader):
    """
    Calculate cross sections from gridded 3D data, with interpolation between start and end.
    """
    def __init__(self, reader, start_point, end_point, n_points):
        GenericSectionReader.__init__(self, reader,
                                      np.linspace(start_point[0], end_point[0], n_points),
                                      np.linspace(start_point[1], end_point[1], n_points))
        

    
class ColumnReader(gradsfile.GriddedDataReader):
    """
    Computes column integrals on the fly *assuming* that the input
    vertical dimension is the midpoint of each layer in meters, with
    the lowest layer next to the ground.
    
    """
    def __init__(self, reader, vertical_dimension=2, zmax=None, zmin=None):
        self.reader = reader
        if reader.ndims() != 3:
            raise ValueError('Reader must be 3d')
        z = self.reader.z()
        thickness = np.array(verticals.thickness(z))
        self._nz = len(z)
        self.vertdim = vertical_dimension
        self.thickness = thickness
        boundaries = np.array(verticals.boundaries(thickness))
        self.zmax = zmax
        self.zmin = zmin
        if zmin is None:
            zmin = 0.0
        if zmax is None:
            zmax = boundaries[-1]
        self.thickness = []
        for lower, upper in zip(boundaries[:-1], boundaries[1:]):
            bottom = max(lower, zmin)
            top = min(upper, zmax)
            self.thickness.append(max(top-bottom, 0.0))
        self.thickness = np.array(self.thickness)
        self._shape = None
        self._set_methods()
        
    def _get_shape(self, input_shape):
        ndims = len(input_shape)
        shape = self.vertdim*(1,) + (self._nz,)
        if ndims > self.vertdim + 1:
            shape += (ndims-self.vertdim-1)*(1,)
        self._weight_shape = shape

        shape = list(input_shape)
        shape[self.vertdim] = 1
        self._output_shape = shape
        
    def _set_methods(self):
        self.x, self.y = self.reader.x, self.reader.y
        self.t, self.dt = self.reader.t, self.reader.dt
        self.nvars = self.reader.nvars
        self.undef = self.reader.undef
        self.tell = self.reader.tell
        self.seek = self.reader.seek
        self.rewind = self.reader.rewind
        self.goto = self.reader.goto
        self.close = self.reader.close
        
    def ndims(self):
        return self.reader.ndims() - 1

    def z(self):
        return np.array([-99])
    
    def read(self, n=1, squeeze=True):
        data_3d = self.reader.read(n, False)
        if not self._shape:
            self._get_shape(data_3d.shape)
        weights = self.thickness.reshape(self._weight_shape)
        columns = np.sum(weights*data_3d, self.vertdim).reshape(self._output_shape)
        if squeeze:
            columns = columns.squeeze()
        return columns

    def describe(self):
        if self.zmax is None:
            descr = 'COLUMN(%s)' % self.reader.describe()
        else:
            descr = 'COLUMN(%s, %f)' % (self.reader.describe(), self.zmax)
        return descr

class VertInterpolator(gradsfile.GriddedDataReader, gradsfile.GenericWrapper):
    """
    Interpolate linearly between the vertical coordinates given by the reader.z() method.
    """
    def __init__(self, reader, interp_points):
        methods = 'read', 'z', 'describe'
        gradsfile.GenericWrapper.__init__(self, reader, methods)
        self.reader = reader
        self._intpt = np.array(interp_points)
        if reader.ndims() != 3:
            raise ValueError('Reader must be 3d')
        reader_points = reader.z()
        zmin, zmax = reader_points.min(), reader_points.max()
        if np.any(self._intpt < zmin):
            raise ValueError('Cannot interpolate request below reader min')
        if np.any(self._intpt > zmax):
            raise ValueError('Cannot interpolate request above reader max')
        self._ind_in = []
        self._weight_in = []
        for zz in interp_points:
            for ind_lev, lev_val in enumerate(reader_points[:-1]):
                next_lev_val = reader_points[ind_lev+1]
                if lev_val <= zz <= next_lev_val:
                    self._ind_in.append(ind_lev)
                    lev_dist = next_lev_val - lev_val
                    weight_up = (zz - lev_val) / lev_dist
                    weight_down = 1 - weight_up
                    #print lev_val, next_lev_val, weight_down, weight_up

                    self._weight_in.append((weight_down, weight_up))

        self._shape = nx, ny, nz = len(reader.x()), len(reader.y()), len(self._intpt)
        
    def z(self):
        return self._intpt

    def read(self, num_steps=1, squeeze=True):
        fld_nd = self.reader.read(num_steps, squeeze=False)
        shape_in = fld_nd.shape
        shape_out = (shape_in[0], shape_in[1], self._shape[2]) + shape_in[3:]
        interpolated = np.empty(shape_out, dtype=fld_nd.dtype, order='f')
        #mask = np.zeros(dtype=bool, order='f')
        for ind_out, (ind_down, (weight_down, weight_up)) in enumerate(zip(self._ind_in, self._weight_in)):
            interpolated[:,:,ind_out,...] = (weight_down * fld_nd[:,:,ind_down,...]
                                             + weight_up * fld_nd[:,:,ind_down+1,...])
            #mask[:,:,ind_out,...] = np.logical_or(fld_nd.mask[:,:,ind_down,...]
        if squeeze:
            interpolated = interpolated.squeeze()
        return interpolated

    def describe(self):
        return 'VERTICAL_INTERPOLATION(%s)' % self.reader.describe()

    
class HorizInterpolatorBasemap(gradsfile.GriddedDataReader, gradsfile.GenericWrapper):
    """
    A horizontal interpolator. The desired output points are given in
    either 1D or 2D arrays, but must be defined in the native grid of
    the given GriddedDataReader. To allow curvilinear output grids,
    the x() and y() methods return only grid indices regardless of
    x_out, y_out arguments. The meshgrid() method however uses the
    actual coordinates.
    
    """
    
    def __init__(self, reader, x_out, y_out):
        methods = 'x', 'y', 'coordinates', 'indices', 'read'
        gradsfile.GenericWrapper.__init__(self, reader, methods)
        self.xout, self.yout = x_out, y_out
        self.reader = reader
        #assert x_out.shape == y_out.shape
        if len(x_out.shape) > 1:
            self.x_interp, self.y_interp = x_out, y_out
        else:
            X, Y = np.meshgrid(x_out, y_out)
            self.x_interp, self.y_interp = X.T, Y.T
        self._xout, self._yout = x_out, y_out # need also these to get indices()
        if hasattr(reader, 'proj'):
            raise UserWarning('Non-trivial projections are not supported')
        # here we could project x_interp & y_interp...
        self.xin, self.yin = self.reader.x(), self.reader.y()
        self._shape2d = self.x_interp.shape

    # These return GRID COORDINATES because geographical coords require 2d arrays -> the
    # coordinates method!
    def x(self):
        return np.arange(self.x_interp.shape[0])
    def y(self):
        return np.arange(self.y_interp.shape[1])
    
    def meshgrid(self):
        raise ValueError('Not implemented use coordinates()')
        #return self.x_interp.T, self.y_interp.T

    def coordinates(self, vertices=False):
        if len(self.xout.shape) > 1 and vertices:
            raise ValueError("Can't have vertices and 2d output points")
        elif vertices:
            xvert, yvert = util.points_to_edges(x_out), util.points_to_edges(y_out)
            X, Y = np.meshgrid(xvert, yvert)
            return X.T, Y.T
        # sorry vertices not working: hard to do if x_out, y_out are 2d.
        return self.x_interp, self.y_interp
    
    def indices(self, lon, lat):
        # override the default method so that we get the indices wrt to the geographical coordinates
        # also appreciate that xout, yout might not be equidistant and not 1d.

        X, Y = self.coordinates()
        # now we assume some lon-lat semantics for the X and Y. Linear approximation
        # shouldn't be too bad though.
        metric_x = np.cos(lat/180. * np.pi)
        metric_y = 1.0
        dist = np.sqrt(((lon-X)*metric_x)**2 + ((lat-Y)*metric_y)**2)
        ind_min_1d = np.argmin(dist)
        ind_min_x, ind_min_y = np.unravel_index(ind_min_1d, dist.shape)
        return ind_min_x, ind_min_y
    
    def read(self, n=1, squeeze=True):
        dataread = self.reader.read(n, squeeze=False)
        # 5d array
        num_times = dataread.shape[4]
        num_vars = dataread.shape[3]
        num_levs = dataread.shape[2]
        #print num_times, num_vars, num_levs
        shape_out = (self._shape2d[0], self._shape2d[1], num_levs, num_vars, num_times)
        dataout = np.ma.MaskedArray(np.empty(shape_out, dtype=np.float32), False)
        for ind_time in range(num_times):
            for ind_var in range(num_vars):
                for ind_lev in range(num_levs):
                    # transpose -- interp assumes y in 1st index, x in 2nd
                    fld = dataout[:,:,ind_lev,ind_var,ind_time]
                    fldmask = dataout.mask[:,:,ind_lev,ind_var,ind_time]

                    #interp = basemap.interp(dataread[:,:,ind_lev,ind_var,ind_time].T,
                    #                     self.xin, self.yin, self.x_interp.T, self.y_interp.T,
                    #                     masked=True, order=1).T

                    interp = basemap.interp(dataread[:,:,ind_lev,ind_var,ind_time],
                                         self.yin, self.xin, self.y_interp, self.x_interp,
                                         masked=True, order=1)

                    
                    
                    #print np.any(interp.mask)
                    fld[...] = interp
                    fldmask[...] = interp.mask
        if squeeze:
            dataout = dataout.squeeze()
        #print np.any(dataout.mask)
        return dataout

    def describe(self):
        return 'HORIZONTAL_INTERPOLATION(%s)' % self.reader.describe()

class HorizInterpolatorFast(HorizInterpolatorBasemap):
    def __init__(self, *args, **kwargs):
        HorizInterpolatorBasemap.__init__(self, *args)

        if 'allow_outside' in kwargs and kwargs['allow_outside']:
            self.interpolator = interp2d.Interp2D(self.xin, self.yin, self.x_interp, self.y_interp,
                                                  kwargs['allow_outside'], kwargs['fill_value'])
        else:
            self.interpolator = interp2d.Interp2D(self.xin, self.yin, self.x_interp, self.y_interp)
        
    def read(self, n=1, squeeze=True):
        dataread = self.reader.read(n, squeeze=False)
        #if hasattr(dataread, 'mask'):
        #    raise NotImplementedError("doesn't work with masked arrays")
        # 5d array
        num_times = dataread.shape[4]
        num_vars = dataread.shape[3]
        num_levs = dataread.shape[2]
        #print num_times, num_vars, num_levs
        shape_out = (self._shape2d[0], self._shape2d[1], num_levs, num_vars, num_times)
        dataout = np.ma.MaskedArray(np.empty(shape_out, dtype=np.float32), False)
        for ind_time in range(num_times):
            for ind_var in range(num_vars):
                interpolated = self.interpolator.interpolate(dataread[:,:,:,ind_var, ind_time])
                dataout[:,:,:,ind_var,ind_time] = interpolated
        if squeeze:
            dataout = dataout.squeeze()
        #print np.any(dataout.mask)
        return dataout
    
if have_interp2d:
    HorizInterpolator = HorizInterpolatorFast
else:
    HorizInterpolator = HorizInterpolatorBasemap
    
class AveragedReaderV2(gradsfile.GriddedDataReader, gradsfile.GenericWrapper):
    def __init__(self, reader, averaging_window, averaging_start=None, end_time_is_valid=True):
        """
        Create an object which reads fields and averages them in time windows
        t...t+averaging_window, optionally starting from averaging_start.

        The valid time for the averaged fields is end of the interval, which is consistent
        with what is assumed in timeseries BUT means the average for a given day will have
        timestamp for day+1 midnight.
        """

        self.window = averaging_window
        my_methods = 'read', 't', 'goto', 'dt', 'seek', 'rewind', 'tell'
        gradsfile.GenericWrapper.__init__(self, reader, my_methods)
        self.times = []
        if not averaging_start:
            averaging_start = reader.t()[0]
        elif averaging_start < reader.t()[0]:
            raise ValueError('Averaging start not be before first time of reader')
        aver_end = averaging_start + self.window
        #print aver_end
        for time in reader.t():
            if time > aver_end - self.window:
                self.times.append(aver_end)
                aver_end = aver_end + self.window
        self._time = self.times[0]
        self._shape = len(reader.x()), len(reader.y()), len(reader.z()), reader.nvars()
        self._end_time_valid = end_time_is_valid
        
    def describe(self):
        return 'TIME_AVERAGE(%s, %s)' % (self.wrapped.describe(), str(self.window))

    def tell(self):
        if self._time > self.times[-1]:
            raise gradsfile.EndOfFileException()
        if not self._end_time_valid:
            return self._time - self.window
        else:
            return self._time

    def t(self):
        if self._end_time_valid:
            return self.times
        else:
            return [time - self.window for time in self.times]

    def dt(self):
        return self.window

    def goto(self, time):
        if not self._end_time_valid:
            time = time + self.window
        if not time in self.times:
            raise gradsfile.GradsfileError('Only exact times are allowed')
        self._time = time

    def seek(self, num_steps):
        ind = self.times.index(self._time)
        if ind + num_steps >= len(self.times):
            raise EndOfFileException()
        self._time = self.times[ind + num_steps]
        
    def rewind(self):
        self._time = self.times[0]
        self.wrapped.rewind()
        
    def read(self, n=1, squeeze=True):
        fld_aver = np.zeros(self._shape + (n,))
        fld_mask = np.zeros(self._shape + (n,), dtype=bool)
        #print 'ttt', self._time, self.t()[-1]
        if self._time > self.times[-1]:
            raise gradsfile.EndOfFileException()
        for step in range(n):
            count = 0
            aver_start = self._time - self.window
            aver_end = self._time
            print 'Averaging:', aver_start, 'to', aver_end, 'now', self.wrapped.tell()
            
            #if self.wrapped.tell() > aver_end:
            #    self.wrapped.rewind()
            while self.wrapped.tell() <= aver_start:
                self.wrapped.seek(1)
            try:
                while aver_start < self.wrapped.tell() <= aver_end:
                    print '...read:', self.wrapped.tell()
                    fld_inst = self.wrapped.read(squeeze=False)
                    fld_inst.shape = fld_inst.shape[:-1]
                    fld_aver[...,step] += fld_inst
                    fld_mask[...,step] = np.logical_or(fld_mask[...,step], fld_inst.mask)
                    count += 1
            except gradsfile.EndOfFileException:
                print 'eof reached'
                pass # The window is possibly not fully covered, but still valid.
            fld_aver[...,step] /= count
            self._time += self.window
        fld_aver = np.ma.MaskedArray(fld_aver, fld_mask)
        if squeeze:
            fld_aver = fld_aver.squeeze()
        return fld_aver
    
AveragedReader = AveragedReaderV2
class MonthlyMeanReader(gradsfile.GriddedDataReader, gradsfile.GenericWrapper):
    def __init__(self, reader, first_time=None, last_time=None, end_time_is_valid=False):
        my_methods = 'read', 't', 'goto', 'dt', 'seek', 'rewind', 'tell'
        gradsfile.GenericWrapper.__init__(self, reader, my_methods)
        self.aver_starts = []
        self.aver_ends = []
        for time in reader.t():
            if first_time and time < first_time:
                continue
            if last_time and time > last_time:
                break
            if not self.aver_starts or self.aver_starts[-1].month != time.month:
                self.aver_starts.append(dt.datetime(time.year, time.month, 1))
                if time.month < 12:
                    aver_end = dt.datetime(time.year, time.month+1, 1)
                else:
                    aver_end = dt.datetime(time.year+1, 1, 1)
                self.aver_ends.append(aver_end)

        
        self._shape = len(reader.x()), len(reader.y()), len(reader.z()), reader.nvars()
        self.windows = zip(self.aver_starts, self.aver_ends)
        
        self._window_ind = 0
        self.end_valid = end_time_is_valid
        if self.end_valid:
            self.times = self.aver_ends
        else:
            self.times = self.aver_starts
            
    def dt(self):
        return gradsfile.GradsDescriptor.MONTHLY

    def t(self):
        return self.times
    
    def describe(self):
        return 'MONTHLY_AVERAGE(%s)' % self.wrapped.describe()
    
    def tell(self):
        if len(self.windows) <= self._window_ind:
                raise gradsfile.EndOfFileException()
        window = self.windows[self._window_ind]
        return window[1] if self.end_valid else window[0]

    def rewind(self):
        self._window_ind = 0
        self.wrapped.rewind()
        
    def goto(self, when):
        found = False
        # looping over the windows could be slow, but for monthly it is probably not
        # significant
        for ind, (time_start, time_end) in enumerate(self.windows):
            found = time_end == when if self.end_valid else time_start == when 
            if found:
                self._window_ind = ind
                return
        if not found:
            raise gradsfile.GradsfileError('Time not available: %s' % str(when))

    def seek(self, num_steps):
        if self._window_ind + num_steps >= len(self.windows):
            raise gradsfile.EndOfFileException()
        self._window_ind += num_steps
        
    def _posit_wrapped(self, when):
        wrapped_times = self.wrapped.t()
        for time_first in wrapped_times:
            if time_first > when:
                break
        self.wrapped.goto(time_first)
        
    def read(self, n=1, squeeze=True):
        fld_aver = np.zeros(self._shape + (n,))
        fld_mask = np.zeros(self._shape + (n,), dtype=bool)
        
        for step in range(n):
            count = 0
            if len(self.windows) <= self._window_ind:
                raise gradsfile.EndOfFileException()
            aver_start, aver_end = self.windows[self._window_ind]
            print 'Averaging:', aver_start, aver_end, self.wrapped.tell()
            if self.wrapped.tell() > aver_end:
                self.wrapped.rewind()
            self._posit_wrapped(aver_start)
            try:
                if self.wrapped.tell() > aver_end:
                    self.wrapped.rewind()
                #while self.wrapped.tell() <= aver_start:
                #    self.wrapped.seek(1)
                assert self.wrapped.tell() > aver_start
                    
                while aver_start < self.wrapped.tell() <= aver_end:
                    print '...read:', self.wrapped.tell()
                    fld_inst = self.wrapped.read(squeeze=False)
                    fld_inst.shape = fld_inst.shape[:-1]
                    fld_aver[...,step] += fld_inst
                    fld_mask[...,step] = np.logical_or(fld_mask[...,step], fld_inst.mask)
                    count += 1
            except gradsfile.EndOfFileException:
                print 'eof reached'
                pass # The window is possibly not fully covered, but still valid.
            fld_aver[...,step] /= count
            self._window_ind += 1
        fld_aver = np.ma.MaskedArray(fld_aver, fld_mask)
        if squeeze:
            fld_aver = fld_aver.squeeze()
        return fld_aver

class CentreOfMassFilter(gradsfile.GenericWrapper, gradsfile.GriddedDataReader):
    """
    Compute the vertical centre of mass for 3D data. Returns masked arrays where mask is
    set if the column amount is less than a given threshold.
    """
    def __init__(self, reader, threshold=1e-15):
        my_methods = 'read', 'z', 'ndims'
        if reader.ndims() != 3:
            raise ValueError('Need 3D Reader')
        gradsfile.GenericWrapper.__init__(self, reader, my_methods)
        self._zin = np.asarray(reader.z())
        self._dzin = verticals.thickness(self._zin)
        self.reader = reader
        self.threshold = threshold
        
    def z(self):
        return np.array([np.nan])
    def ndims(self):
        return 2

    def read(self, *args, **kwargs):
        fld = self.reader.read(*args, **kwargs)
        # fld[x,y,z][nt,nvars]
        new_shape = (fld.shape[0], fld.shape[1]) + fld.shape[3:]
        moment = np.zeros(new_shape)
        mass = np.zeros(new_shape)
        for iz, z in enumerate(self._zin):
            layer_mass = fld[:,:,iz,...]*self._dzin[iz]
            moment += layer_mass*z
            mass += layer_mass

        not_mask = mass > self.threshold
        cm = np.zeros(moment.shape)
        cm[not_mask] = moment[not_mask] / mass[not_mask]
        return np.ma.MaskedArray(cm, np.logical_not(not_mask))

    def describe(self):
        return '%s(%s)'  % ('CENTRE_OF_MASS', self.reader.describe())
    
class TimeInterpolator(gradsfile.GenericWrapper, gradsfile.GriddedDataReader):
    def describe(self):
        return '%s(%s)'  % ('TIME_INTERPOLATION', self.reader.describe())
    
    def __init__(self, reader, times):
        my_methods = 'read', 't', 'tell', 'goto', 'seek', 'rewind', 'dt', 'describe'
        gradsfile.GenericWrapper.__init__(self, reader, my_methods)
        reader_times = reader.t()
        for time in times:
            if not reader_times[0] <= time <= reader_times[-1]:
                raise ValueError('Requested time not covered by reader')
        self._times = times
        self.reader = reader
        self._time = times[0]
        self._shape = reader.shape()
        self._eof = False
        
    def t(self):
        return self._times

    def tell(self):
        if self._eof:
            raise gradsfile.EndOfFileException()
        return self._time

    def goto(self, time):
        if not time in self._times:
            raise gradsfile.GradsfileError('Time not available')
        self._time = time
        self._eof = False
        
    def seek(self, num_steps):
        ind = self._times.index(self._time)
        if ind + num_steps >= len(self._times):
            raise EndOfFileException()
        self._time = self._times[ind + num_steps]

        
    def rewind(self):
        self._time = self._times[0]
        self._eof = False
        
    def read(self, n=1, squeeze=True):
        if self._eof:
            raise gradsfile.EndOfFileException()
        
        ind_time = self.t().index(self._time)
        if ind_time + n > len(self._times):
            raise Exception()
        
        ret_arr = np.empty(self._shape + (n,))
        mask = np.empty(self._shape + (n,))
        for step in range(n):
            time = self._times[ind_time+step]
            for ind, time_past in enumerate(self.reader.t()):
                time_future = self.reader.t()[ind+1]
                if time_past < time <= time_future:
                    break
            step_sec = float(util.dt2sec(time_future - time_past))
            weight_futr = util.dt2sec(time-time_past) / step_sec
            weight_past = 1 - weight_futr
            #print 'current, past, future', time, time_past, time_future, weight_past, weight_futr

            self.reader.goto(time_past)
            fld_past = self.reader.read(squeeze=False)
            self.reader.goto(time_future)
            fld_futr = self.reader.read(squeeze=False)
            fld_interp = weight_past*fld_past + weight_futr*fld_futr
            ret_arr[...,step] = fld_interp[...,0] # last dim is time since no squeeze
            mask[...,step] = fld_interp.mask[...,0]
        self._eof = (ind_time + n == len(self._times))

        ret_arr = np.ma.MaskedArray(ret_arr, mask)
        if n == 1:
            ret_arr = ret_arr.reshape(self._shape)
        if squeeze:
            ret_arr = ret_arr.squeeze()
        #print ret_arr
        return ret_arr
    
class DailyMaxFilter(AveragedReader, gradsfile.GenericWrapper):
    def __init__(self, reader, end_time_is_valid=False):#, averaging_window, averaging_start=None):
        """
        Create an object returning daily max fields for a given stream. Valid time is the
        midnight for the day maximum is taken.
        """
        window = dt.timedelta(hours=24)
        rd_start = reader.t()[0]
        start_time = dt.datetime(rd_start.year, rd_start.month, rd_start.day)
        AveragedReader.__init__(self, reader, window, start_time, end_time_is_valid)

    def describe(self):
        description = 'DAILY_MAX(%s)' % (self.wrapped.describe())
        return description
    
    def read(self, n=1, squeeze=True):
        outbuf = np.zeros(self._shape + (n,))
        mask = None
        #print 'ttt', self._time, self.t()[-1]
        if self._time > self.times[-1]:
            raise gradsfile.EndOfFileException()
        for step in range(n):
            count = 0
            aver_start = self._time - self.window
            aver_end = self._time
            print 'Maxxing:', aver_start, 'to', aver_end, 'now', self.wrapped.tell()
            
            #if self.wrapped.tell() > aver_end:
            #    self.wrapped.rewind()
            while self.wrapped.tell() <= aver_start:
                self.wrapped.seek(1)
            try:
                while aver_start < self.wrapped.tell() <= aver_end:
                    print '...read:', self.wrapped.tell()
                    fld_inst = self.wrapped.read(squeeze=False)
                    fld_inst.shape = fld_inst.shape[:-1]
                    if count == 0:
                        outbuf[...,step] = fld_inst
                    else:
                        outbuf[...,step] = np.maximum(outbuf[...,step], fld_inst)
                    if mask is not None:
                        mask[...,step] = np.logical_or(mask[...,step], fld_inst.mask)
                    elif hasattr(fld_inst, 'mask'):
                        mask = np.zeros(outbuf.shape)
                        mask[...,step] = fld_inst.mask
                    count += 1
            except gradsfile.EndOfFileException:
                print 'eof reached'
                pass # The window is possibly not fully covered, but still valid.
            self._time += self.window
        if mask is not None:
            outbuf = np.ma.MaskedArray(outbuf, mask)
        if squeeze:
            outbuf = outbuf.squeeze()
        return outbuf

class DensityReader(gradsfile.GenericWrapper):
    def __init__(self, reader, areas=None):
        methods = 'read',
        gradsfile.GenericWrapper.__init__(self, reader, methods)
        if areas is None:
            self.area = gridtools.area(reader.y(), reader.x())
        else:
            self.area = areas
        self.reader = reader
        if not self.reader.nvars() == 1:
            raise ValueError('Cannot handle reader.nvars != 1')
        
    def read(self, num_steps=1, squeeze=True):
        values = self.reader.read(num_steps, squeeze=False)
        output = np.empty(values.shape)
        for ind_time in range(num_steps):
            for ind_lev in range(values.shape[2]):
                output[:,:,ind_lev,ind_time,0] = values[:,:,ind_lev,ind_time,0] / self.area
        if squeeze:
            output = output.squeeze()
        return output

    def describe(self):
        descr = 'DENSITY(%s)' % self.reader.describe()
        return descr
    
class EnsembleOperation:
    def __init__(self, op, name):
        self.name = name
        self.op = op

    def __call__(self, field_iter):
        return self.op(field_iter)

def _fld_mean(fld_iter):
    acc = count = 0
    mask = None
    for fld in fld_iter:
        acc += fld
        count += 1
    return acc / count

def _fld_stdev(fld_iter):
    acc = count = 0
    mask = None
    flds = list(fld_iter)
    for fld in flds:
        acc += fld
        count += 1
    mean = acc / count
    acc = 0
    for fld in flds:
        acc += (fld - mean)**2
    stdev = np.sqrt(acc/count)
    return stdev
    
class EnsembleProcessor(gradsfile.GriddedDataReader):
    mean = EnsembleOperation(_fld_mean, 'ens_mean')
    stdev = EnsembleOperation(_fld_stdev, 'ens_stdev')

    def __init__(self, readers, op):
        self.op = op
        self.readers = readers
        self.leader = leader = readers[0]
        self.x = leader.x
        self.y = leader.y
        self.z = leader.z
        self.t = leader.t
        self.ndims = leader.ndims
        self.dt = leader.dt
        self.nvars = leader.nvars
        self.undef = op(reader.undef() for reader in readers)
        self._chk_conform()
        
    def _chk_conform(self):
        if not all(reader.shape() == self.leader.shape() for reader in self.readers):
            raise ValueError('Shapes differ')
        for reader in self.readers:
            if not reader.t() == self.leader.t():
                raise ValueError('Time axes differ')
        for method in 'x', 'y', 'z':
            output_leader = getattr(self.leader, method)()
            all_leader_isnan = all(np.isnan(output_leader))
            for reader in self.readers:
                output_reader = getattr(reader, method)()
                if not ((all_leader_isnan and all(np.isnan(output_reader)))
                        or all(util.almost_eq(output_leader, output_reader, 1e-2))):
                    print output_leader, output_reader
                    raise ValueError('Coordinates differ')
        
    def _operate(self, method, *args):
        for reader in self.readers:
            getattr(reader, method)(*args)
        
    def _read_iter(self, *args):
        for reader in self.readers:
            yield reader.read(*args)
        
    def tell(self):
        return self.leader.tell()

    def goto(self, when):
        for reader in self.readers:
            reader.goto(when)

    def read(self, count=1, squeeze=True):
        return self.op(self._read_iter(count, squeeze))

    def close(self):
        self._operate('close')

    def seek(self, n=1):
        self._operate('seek', n)

    def rewind(self):
        self._operate('rewind')

    def describe(self):
        if len(self.readers) < 3:
            reader_text = ','.join(reader.describe() for reader in self.readers)
        else:
            reader_text = ','.join(reader.describe() for reader in self.readers[:2]) + ',...'
        text = '%s(%s; n=%i)' % (self.op.name, reader_text, len(self.readers))
        return text
        
        
    
    
            
TIME_DIM = -2
class DummyReader(gradsfile.GriddedDataReader):
    """
    A mock reader object, which reads from an array it is initialized with.
    """
    def __init__(self, array, times=None, x=None, y=None, z=None):
        assert len(array.shape) > 2
        # add the variable dimension, although unity
        self.ar = array.reshape(array.shape + (1,))
        
        if x is not None:
            self._x = np.array(x)
        else:
            self._x = np.linspace(0, 1, self.ar.shape[0])
        if y is not None:
            self._y = np.array(y)
        else:
            self._y = np.linspace(0, 1, self.ar.shape[1])

        if z is not None:
            self._z = np.array(z)
        elif len(array.shape) > 3:
            self._z = np.linspace(0, 1, self.ar.shape[2])
        else:
            self._z = np.array([np.nan])

        self._num_dims = 3 if len(array.shape) > 3 else 2
            
        if times is not None:
            self.times = times
        else:
            first_time = dt.datetime(2000,1,1)
            step = util.one_hour
            self.times = [first_time + ii*step for ii in xrange(self.ar.shape[TIME_DIM])]
        self.step = 0
        
    def x(self):
        return self._x
    def y(self):
        return self._y
    def z(self):
        return self._z
    def t(self):
        return self.times
    def ndims(self):
        return self._num_dims
    def undef(self):
        return np.nan
    
    def seek(self, n):
        if self.step + n < self.ar.shape[TIME_DIM]:
            self.step += n
        else:
            raise gradsfile.EndOfFileException()
    def tell(self):
        try:
            return self.times[self.step]
        except IndexError:
            raise gradsfile.EndOfFileException()
        
    def close(self):
        pass

    def rewind(self):
        self.step = 0

    def goto(self, when):
        try:
            self.step = self.times.index(when)
        except ValueError:
            raise gradsfile.GradsfileError('Failed to goto')
        
    def describe(self):
        return 'DummyReader'
            
    def read(self, n=1, squeeze=True):
        end_step = self.step + n
        if self.ar.shape[TIME_DIM] < end_step:
            self.step = end_step
            raise gradsfile.EndOfFileException()
        if squeeze:
            ret_ar = self.ar[...,self.step:end_step,0:1].squeeze()
        else:
            ret_ar = self.ar[...,self.step:end_step,0:1]
        self.step = end_step
        return np.ma.MaskedArray(ret_ar, False)

    def nvars(self):
        return 1
    
