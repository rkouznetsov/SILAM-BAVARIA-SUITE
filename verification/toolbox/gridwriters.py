from support import pupynere as netcdf
from toolbox import util, ncreader, gradsfile
import datetime as dt
from os import path
import os

_max_sec_i32 = 2**31
sec_in_day = 3600*24

"""
Classes for writing netcdf files.

NCWriter -> single netcdf file
NCDatasetWriter -> time-split files + .ctl file

Contrary to the griddedDataReaders, the writer objects handle multiple variables at
time. The variables can be 2d or 3d, but must share the same time dimension.

Usage:

writer = NCWriter(file_to_write)
writer.add_var('stuff', [dt.datetime(...), dt.datetime(...)...], x, y, z)
writer.write('stuff', 0, stuff_data(t=0))
writer.write('stuff', 1, stuff_data(t=1))
...
writer.close()
For NCWriter the timesteps can be written in any order,
for NCDatasetWriter they must be written sequentially.

The Endpoint classes define and object which with the readwrite() method runs through set
of readers, reading the given timesteps and immediately writing into a file.

"""


class CFMD:
    """
    CF Metadata descriptor. Use by calling the factory method.
    """
    @staticmethod
    def geo_lonlat_height():
        md = CFMD()
        md.yattr = {'long_name' : 'latitude',
                    'standard_name' : 'latitude',
                    'units' : 'degrees_north'}
        md.xattr = {'long_name' : 'longitude',
                    'standard_name' : 'longitude',
                    'units' : 'degrees_east'}
        md.zattr = {'units' : 'meters', 'positive' : 'up'}
        return md
    
    def getattr(self, dim):
        try:
            dct = {'x':self.xattr, 'y':self.yattr, 'z':self.zattr}
            return dct[dim]
        except KeyError:
            raise ValueError('Strange dimension %s' % dim)
        
    def putattr(self, dim, obj):
        dct = self.getattr(dim)
        for key, val in dct.items():
            setattr(obj, key, val)
        
class NCWriter:
    def _nctime(self, time):
        delta = time - self.ref_time
        seconds = delta.days*sec_in_day + delta.seconds
        if not -_max_sec_i32 < seconds < _max_sec_i32:
            raise ValueError('Seconds overflow')
        return seconds
    
    def __init__(self, filename, cf_attr=None):
        self._ncfile = netcdf.netcdf_file(filename, 'w', version=2)
        self._xdim = 'x'
        self._ydim = 'y'
        self._zdim = 'z'
        self._times = None
        self._dims = {}
        self._varnames = []
        self._cf = cf_attr
        
    def _make_cf(self, dimname, dimvar):
        self._cf.putattr(dimname, dimvar)

    def set_cf_attr(self, cf_attr):
        if len(self._dims) > 0:
            raise ValueError('Cannot change CF attributes after defining variables')
        self._cf = cf_attr
        
    def _set_dim(self, dimname, dimdata):
        if dimname in self._dims:
            if not all(util.almost_eq(self._dims[dimname], dimdata)):
                raise ValueError('Inconsistent dimension: %s' % dimname)
        else:
            self._dims[dimname] = dimdata 
            self._ncfile.createDimension(dimname, len(dimdata))
            dimvar = self._ncfile.createVariable(dimname, 'f', (dimname,))
            dimvar[:] = dimdata
            dimvar.axis = dimname
            if self._cf:
                self._make_cf(dimname, dimvar)
                
    def _set_time(self, times):
        if self._times is not None:
            if self._times != times:
                raise ValueError('Inconsistent time dimension')
        else:
            self._times = times
            self._ncfile.createDimension('time', len(times))
            self.ref_time = times[0]
            timevar = self._ncfile.createVariable('time', 'i', ('time',))
            timevar[:] = [self._nctime(time) for time in times]
            timevar.units = ncreader.make_epoch(self.ref_time, 'seconds')
            self._step = 0
        
    def _makevar(self, varname, ndims):
        self._varnames.append(varname)
        if ndims == 2:
            new_var = self._ncfile.createVariable(varname, 'f', ('time', self._ydim, self._xdim))
        elif ndims == 3:
            new_var = self._ncfile.createVariable(varname, 'f',
                                                  ('time', self._zdim, self._ydim, self._xdim))
        return new_var
    
    def add_var(self, varname, times, x, y, z=None, attrs={}):
        """
        Add a variable. Must be called before writing. For 2D variables, set z=None.
        """
        if varname in self._varnames:
            raise ValueError('Variable %s already defined' % varname)
        self._set_dim(self._xdim, x)    
        self._set_dim(self._ydim, y)            
        self._set_time(times)
        if z is not None:
            self._set_dim(self._zdim, z)
        new_var = self._makevar(varname, 2 if z is None else 3)
        for attr, value in attrs.items():
            setattr(new_var, attr, value)
        return new_var
    
    def write(self, varname, step, data):
        self._ncfile.variables[varname][step,...] = data.T
        
    def close(self):
        self._ncfile.close()

    def add_attr(self, attr_name, attr_value):
        setattr(self._ncfile, attr_name, attr_value)
        
class NCEndpoint(NCWriter):
    def __init__(self, filename, times=None, attrs={}):
        NCWriter.__init__(self, filename)
        self._readers = {}
        self.times = times
        for key, val in attrs.items():
            self.add_attr(key, val)
        
    def add(self, varname, reader, attrs={}):
        z = None if reader.ndims() == 2 else reader.z()
        times = self.times or reader.t()
        new_var = NCWriter.add_var(self, varname, times, reader.x(), reader.y(), z, attrs)
        self._readers[varname] = reader
        new_var.describe = reader.describe()
        
    def readwrite(self, verbose=False):
        vars_readers = self._readers.items()
        # get the times from what the writer class has made of it, if otherwise undefined
        times = self.times or self._times
        for ind_time, time in enumerate(times):
            if verbose:
                print time
            for varname, reader in vars_readers:
                reader.goto(time)
                fld = reader.read()
                NCWriter.write(self, varname, ind_time, fld)
        NCWriter.close(self)
    

class NCDatasetWriter:
    def __init__(self, dirname, split, cf_attr=None):
        allowed_splits = 'daily monthly hourly'
        if not split in allowed_splits:
            raise ValueError('split must be one of: %s' % allowed_splits)
        self.dirname = dirname
        try:
            os.makedirs(dirname)
        except OSError:
            pass
        timetmpl = {'daily':'%y4%m2%d2', 'monthly':'%y4%m2', 'hourly':'%y4%m2%d2%h2'}[split]
        self.nctmpl = '%s_%s.nc' % (path.basename(dirname), timetmpl)
        self._ncwr = None
        self._ncnm = None
        self._var_arg_lst = []
        self._times = None
        self._files_done = set()
        self.cf_attr = cf_attr
        
    def _set_time(self, times):
        if self._times is not None:
            if any(time != time_new for (time, time_new) in zip(self._times, times)):
                raise ValueError('Inconsistent time dimension')
        else:
            self._times = times
            if len(self._times) > 1:
                timestep = self._times[1]-self._times[0]
                if not all(t2-t1 == timestep for (t2,t1) in zip(times[1:], times[:-1])):
                    print times[1:]
                    print times[:-1]
                    print timestep
                    raise ValueError('Non-constant timestep')
                self.timestep = timestep
            else:
                self.timestep = dt.timedelta(minutes=1)
            self._shift4file = {}
            # find the shift from global timestep indices to the per-file ones
            for ind_time, time in enumerate(self._times):
                filename = self._get_filename(time)
                if not filename in self._shift4file:
                    self._shift4file[filename] = ind_time
                    
    def add_var(self, varname, times, x, y, z=None, attrs={}):
        self._var_arg_lst.append((varname, x, y, z, attrs))
        self._set_time(times)

    def _get_filename(self, time):
        filename = path.join(self.dirname, gradsfile.expandGradsTemplate(self.nctmpl, time))
        return filename

    def _get_writer(self, step):
        time = self._times[step]
        filename = self._get_filename(time)
        if filename in self._files_done:
            raise ValueError('NCDatasetWriter requires sequential writes')
            
        if filename != self._ncnm:
            if self._ncwr is not None:
                #print 'done:', step
                self._files_done.add(self._ncnm)
                self._ncwr.close()
            self._ncwr = NCWriter(filename, self.cf_attr)
            self._ncnm = filename
            self._step_shift = self._shift4file[filename]
            print 'new writer:', filename
            times4file = [time for time in self._times if self._get_filename(time) == filename]
            #print filename
            #print times4file
            for varname,  x, y, z, attrs in self._var_arg_lst:
                self._ncwr.add_var(varname, times4file, x, y, z, attrs)
        return self._ncwr
    
    def write(self, varname, step, data):
        writer = self._get_writer(step)
        writer.write(varname, step-self._step_shift, data)

    def _makectl(self):
        ctl_path = path.join(self.dirname, '%s.ctl' % path.basename(self.dirname))
        with open(ctl_path, 'w') as output:
            output.write('dset ^%s\n' % self.nctmpl)
            output.write('options template\n')
            output.write('tdef time %i linear %s %s\n' % (len(self._times),
                                                          gradsfile.makeGradsTime(self._times[0]),
                                                          gradsfile.makeGradsInterval(self.timestep)))
            
        
    def close(self, make_ctl=True):
        self._ncwr.close()
        if make_ctl:
            self._makectl()
    
    
class NCDatasetEndpoint:
    def __init__(self, filename, times=None, split='daily'):
        self._readers = {}
        self.times = times
        self.writer = NCDatasetWriter(filename, split)
                
    def add(self, varname, reader, attrs={}):
        z = None if reader.ndims() == 2 else reader.z()
        times = self.times or reader.t()
        #print 'readt', reader.t()
        attrs['describe'] = reader.describe()
        self.writer.add_var(varname, times, reader.x(), reader.y(), z, attrs)
        self._readers[varname] = reader
        
    def readwrite(self, verbose=False):
        vars_readers = self._readers.items()
        # get the times from what the writer class has made of it, if otherwise undefined
        times = self.times or self.writer._times
        for ind_time, time in enumerate(times):
            if verbose:
                print time
            for varname, reader in vars_readers:
                reader.goto(time)
                fld = reader.read()
                self.writer.write(varname, ind_time, fld)
        self.writer.close()
