from os import path
import numpy as np
from toolbox import namelist, util, silamfile, gradsfile, gridtools, projections
import datetime as dt

class GradsWriter:
    """
    A class for writing grads files, including super ctls.

    Interface is similar to NCWriter in gridwriters module. However, the data must be
    written in following the order of variables and timesteps.

    The general part of super_ctl is generated automatically. The var-specific namelists
    are given at variable definition step.

    """
    
    def __init__(self, ctl_name, binary_name, super_ctl_name=None, undef=-1e15):
        self.super_ctl_name = super_ctl_name
        self.ctl_name = ctl_name
        self.binary_name = binary_name

        self._x = self._y = self._z = self._times = None
        self.var_nls = {}

        if super_ctl_name and self.ctl_name.startswith('^'):
            self._super_ctl_dir = path.dirname(super_ctl_name)
            self._ctl_expanded = path.join(self._super_ctl_dir, self.ctl_name.replace('^',''))
        else:
            self._ctl_expanded = self.ctl_name
        if self.binary_name.startswith('^'):
            ctl_dir = path.dirname(self._ctl_expanded)
            self._bin_expanded = path.join(ctl_dir, self.binary_name.replace('^',''))
        else:
            self._bin_expanded = self.binary_name
        self.undef = undef
        self._vars_3d = set()
        self._def_done = False
        self.var_list = []
        self._closed = False
        self._general_attrs = namelist.Namelist('gattr')
        
    def _set_xy(self, which, dimdata):
        dimdata = np.array(dimdata)
        step = dimdata[1]-dimdata[0]
        if not all(util.almost_eq(dimdata[1:]-dimdata[:-1], step, 1e-5)):
            print dimdata
            raise ValueError('Dimension not linear: %s' % which)
        if which == 'x':
            if self._x is None:
                self._x = (len(dimdata), dimdata[0], step)
            dim = self._x
        elif which == 'y':
            if self._y is None:
                self._y = (len(dimdata), dimdata[0], step)
            dim = self._y
        else:
            raise ValueError('Bad which')
        if dim[0] != len(dimdata):
            raise ValueError('Dimension length differs: %s' % which)
        if not util.almost_eq(dim[1], dimdata[0]):
            raise ValueError('Dimension start differs: %s' % which)
        if not util.almost_eq(dim[2], step):
            raise ValueError('Dimension step differs: %s' % which)
            
    def _set_z(self, dimdata):
        dimdata = np.array(dimdata)
        if self._z is None:
            self._z = dimdata
        else:
            if not all(util.almost_eq(self._z, dimdata)):
                print self._z
                print dimdata
                raise ValueError('Dimension differs: z')

    def _set_time(self, times):
        if self._times is None:
            if len(times) > 1:
                step = times[1] - times[0]
            else:
                step = dt.timedelta(0)
            prev_time = times[0]
            for time in times[1:]:
                if (time - prev_time) != step:
                    raise ValueError('Nonlinear dimension: time')
                prev_time = time
            self._timestep = step
            self._times = times
        else:
            if len(times) != len(self._times):
                raise ValueError('Dimension differs: time')
            if not all(t1==t2 for (t1, t2) in zip(self._times, times)):
                raise ValueError('Dimension differs: time')
        
    def add_var(self, varname, times, x, y, z=None, var_nl=None):
        """
        Add a variable. For 2D variables, set z = None. A namelist is given to be included
        in the super ctl.
        """
        if self._def_done:
            raise ValueError('Not in definition mode')
        if varname in self.var_nls:
            raise ValueError('Variable already defined: %s' % varname)
        self._set_xy('x', x)
        self._set_xy('y', y)
        self._set_time(times)
        if z is not None:
            self._set_z(z)
            self._vars_3d.add(varname)
        if var_nl is not None:
            self.var_nls[varname] = var_nl
        else:
            self.var_nls[varname] = namelist.Namelist(varname)
        self.var_list.append(varname)
        
    def _get_dim_str_xy(self, dim):
        return '%i linear %f %f' % dim

    def _get_dim_str_z(self):
        return '%i levels %s' % (len(self._z), ' '.join('%.3f' % lev for lev in self._z))

    def _get_dim_str_time(self):
        return '%i linear %s %s' % (len(self._times),
                                    gradsfile.makeGradsTime(self._times[0]),
                                    gradsfile.makeGradsInterval(self._timestep))

    def _enddef(self):
        self._def_done = True
        self._var_done = {}
        for var in self.var_list:
            self._var_done[var] = False
        self._step_curr = 0
        self._var_next = self.var_list[0]
        self._filename_curr = None
        self._bin_out = None
        if self._z is None:
            # only 2d
            self._z = np.array([0.0])
            
    def _chk_bin_file(self):
        time = self._times[self._step_curr]
        filename_new = gradsfile.expandGradsTemplate(self._bin_expanded, time)
        if self._filename_curr is None:
            self._bin_out = open(filename_new, 'w')
            self._filename_curr = filename_new
        elif filename_new != self._filename_curr:
            self._bin_out.close()
            self._bin_out = open(filename_new, 'w')
            self._filename_curr = filename_new
            
    def write(self, varname, step, data):
        """
        Write one step for one variable. 2D variables given as 2D array (x,y), 3D
        variables given as 3D array (x,y,z). The order must be the as add_var() has been
        called. Index of timestep must be given and either same or one greater than the
        previous.
        """
        
        if not self._def_done:
            self._enddef()
        if self._closed:
            raise ValueError('File closed')
        if varname != self._var_next:
            raise ValueError('Variable expected: %s; variable given: %s' % (self._var_next, varname))
        if len(self._times) <= step:
            raise ValueError('Writing too many steps')
        if step != self._step_curr:
            raise ValueError('Timestep expected: %i; timestep given: %i' % (self._step_curr, step))
        if varname in self._vars_3d:
            shape = self._x[0], self._y[0], len(self._z)
        else:
            shape = self._x[0], self._y[0]
        if shape != data.shape:
            raise ValueError('Shape expected: %s; shape given: %s' % (str(shape), str(data.shape)))
        self._chk_bin_file()
        if hasattr(data, 'mask'):
            data = data.filled(self.undef)
        data.astype(np.float32).ravel('f').tofile(self._bin_out)
        # all steps before curr_step done
        if varname == self.var_list[-1]:
            self._var_next = self.var_list[0]
            self._step_curr = step + 1
        else:
            ind_var = self.var_list.index(varname)
            self._var_next = self.var_list[ind_var+1]

    def close(self):
        self._bin_out.close()
        self.make_ctl()
        if self.super_ctl_name:
            self.make_super_ctl()
        self._closed = True
        
    def make_ctl(self):
        if not self._def_done:
            raise ValueError('Writer in define mode')
        
        out = open(self._ctl_expanded, 'w')
        out.write('dset %s\n' % self.binary_name)
        out.write('options template\n')
        out.write('undef %g\n' % self.undef)
        out.write('xdef %s\n' % self._get_dim_str_xy(self._x))
        out.write('ydef %s\n' % self._get_dim_str_xy(self._y))
        out.write('zdef %s\n' % self._get_dim_str_z())
        out.write('tdef %s\n' % self._get_dim_str_time())
        out.write('vars %i\n' % len(self.var_list))
        num_levs = len(self._z)
        for var in self.var_list:
            var_num_levs = num_levs if var in self._vars_3d else 0
            out.write('%s %i 99 99 0 %s\n' % (var, var_num_levs, var))
        out.write('endvars\n')
        out.close()
        
    def make_super_ctl(self):
        if not self._def_done:
            raise ValueError('Writer in define mode')
        out = open(self.super_ctl_name, 'w')
        grid = silamfile.SilamLatLonGrid(self._x[1], self._x[2], self._x[0],
                                         self._y[1], self._y[2], self._y[0],
                                         projections.LatLon())
        vertical = silamfile.SilamHeightLevels(self._z)

        out.write('LIST = general\n')
        out.write('ctl_file_name = %s\n' % self.ctl_name)
        grid.as_namelist().tofile(out)
        vertical.as_namelist().tofile(out)
        self._general_attrs.tofile(out)
        out.write('END_LIST = general\n\n')
        for var, nl in self.var_nls.items():
            nl.tofile(out, listname=var)
        out.close()

    def add_attr(self, key, val):
        # add to the general namelist
        self._general_attrs.put(key, val)

class GradsEndpoint:
    grads_class = GradsWriter
    def __init__(self, ctl_name, bin_name, super_ctl_name=None, times=None):
        self.writer = self.grads_class(ctl_name, bin_name, super_ctl_name)
        self._readers = {}
        self.times = times
        self.writer.undef = None
        
    def add(self, varname, reader, attrs={}):
        z = None if reader.ndims() == 2 else reader.z()
        times = self.times or reader.t()
        var_nl = namelist.Namelist(varname)
        var_nl.set('describe', reader.describe())
        for key, val in attrs.items():
            var_nl.set(key, val)
        new_var = self.writer.add_var(varname, times, reader.x(), reader.y(), z, var_nl)
        self._readers[varname] = reader
        if self.writer.undef is None:
            self.writer.undef = reader.undef()
        elif not util.almost_eq(self.writer.undef, reader.undef()):
            raise ValueError('Cannot have different undef for several varibles')
        
    def readwrite(self, verbose=False):
        times = self.times or self.writer._times
        for ind_time, time in enumerate(times):
            if verbose:
                print time
            for varname in self.writer.var_list:
                reader = self._readers[varname]
                reader.goto(time)
                fld = reader.read()
                self.writer.write(varname, ind_time, fld)
        self.writer.close()

                                    
        
class BufferedGradsWriter(GradsWriter):
    def _get_buffer(self):
        # figure out up to how many steps go to a file
        count_steps_max = 0
        filename_prev = None
        for time in self._times:
            filename = gradsfile.expandGradsTemplate(self._bin_expanded, time)
            if filename != filename_prev:
                count_steps = 0
            count_steps += 1
            filename_prev = filename
            if count_steps > count_steps_max:
                count_steps_max = count_steps
        num_vars_3d = len(self._vars_3d)
        num_vars_2d = len(self.var_list) - num_vars_3d
        shape_2d = self._x[0], self._y[0]
        num_fields = (num_vars_3d*len(self._z) + num_vars_2d)*count_steps_max
        buf_shape = self._x[0], self._y[0], num_fields
        self._buffer = np.empty(buf_shape, dtype=np.float32, order='f')
        self._num_fields_buffer = 0
        
    def _enddef(self):
        GradsWriter._enddef(self)
        self._get_buffer()

    def _flush(self):
        self._buffer[...,:self._num_fields_buffer].ravel('f').tofile(self._bin_out)
        self._num_fields_buffer = 0
        
    def _chk_bin_file(self):
        time = self._times[self._step_curr]
        filename_new = gradsfile.expandGradsTemplate(self._bin_expanded, time)
        #print filename_new, self._filename_curr
        if self._filename_curr is None:
            self._bin_out = open(filename_new, 'w')
            self._filename_curr = filename_new
        elif filename_new != self._filename_curr:
            self._flush()
            self._bin_out.close()
            self._bin_out = open(filename_new, 'w')
            self._filename_curr = filename_new

    def close(self):
        self._flush()
        GradsWriter.close(self)

    def write(self, varname, step, data):
        if not self._def_done:
            self._enddef()
        if self._closed:
            raise ValueError('File closed')
        if varname != self._var_next:
            raise ValueError('Variable expected: %s; variable given: %s' % (self._var_next, varname))
        if len(self._times) <= step:
            raise ValueError('Writing too many steps')
        if step != self._step_curr:
            raise ValueError('Timestep expected: %i; timestep given: %i' % (self._step_curr, step))
        if varname in self._vars_3d:
            shape = self._x[0], self._y[0], len(self._z)
            num_fields = shape[2]
            shape_3d = shape
        else:
            shape = self._x[0], self._y[0]
            is_3d = False
            num_fields = 1
            shape_3d = shape + (1,)
        if shape != data.shape:
            raise ValueError('Shape expected: %s; shape given: %s' % (str(shape), str(data.shape)))
        self._chk_bin_file()
        if hasattr(data, 'mask'):
            data = data.filled(self.undef)
        ind_buf_start = self._num_fields_buffer
        ind_buf_end = ind_buf_start + num_fields
        assert ind_buf_end <= self._buffer.shape[2]
        self._buffer[:,:,ind_buf_start:ind_buf_end] = data.astype(np.float32).reshape(shape_3d)
        self._num_fields_buffer += num_fields
        #data.astype(np.float32).ravel('f').tofile(self._bin_out)
        # all steps before curr_step done
        if varname == self.var_list[-1]:
            self._var_next = self.var_list[0]
            self._step_curr = step + 1
        else:
            ind_var = self.var_list.index(varname)
            self._var_next = self.var_list[ind_var+1]


class BufferedGradsEndpoint(GradsEndpoint):
    grads_class = BufferedGradsWriter
