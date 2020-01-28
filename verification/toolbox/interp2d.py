"""
Interpolation routines using the interpf and interpfsp (single precision) Fortran modules.
See class Interp2D for further info.
"""

import interpf, interpfsp
import numpy as np

from toolbox import util

def scale_to_unit(grid_crd, pt_crd, allow_outside=False):
    dxg = grid_crd[1:]-grid_crd[:-1]
    dx = dxg[0]
    if not np.all(util.almost_eq(dx, dxg, 1e-5)):
        print dxg.min(), dxg.max(), dx
        raise ValueError('Grid axis is uneven')
    if dx < 1e-15:
        raise ValueError('dx near zero')
    x0 = grid_crd[0]
    new_pt = (pt_crd-x0) / dx
    mask = _get_mask(len(grid_crd), new_pt, allow_outside)
    return new_pt, mask


def _get_mask(grid_size, crd_in, allow_outside):
    mask = np.logical_or(crd_in < 0, crd_in > grid_size-1)
    if not allow_outside and np.any(mask):
        raise ValueError('Outside values not allowed')
    return mask

class Interp2D:
    def __init__(self, xin, yin, xout, yout, allow_outside=False, fillvalue=None, precision=None):
        """
        Create a 2D linear interpolator.
        Input arguments:
        
        - xin, yin : the grid of input data. The data array will be given in (x,y) or (x,y,z) order.
        
        - xout, yout : the coordinates of output data. Either 2D or 1D arrays of same
          shape. The output will have this shape.

        - allow_outside [False], fillvalue [None] : if allow_outside=True, points outside
          the input grid will be set to fillvalue. Otherwise crash with out-of-grid
          requests.

        - precision : Select single or double precision. By default use double if any of
          the numeric arguments has double precision.

        Using the interpolator:
        interp = Interp2D(xin, yin, ...)
        outputdata = interp.interpolate(inputdata)

        If inputdata is 3D, the interpolation will happen along the first two dimensions.
          
        """
        
        self.xin, self.yin = xin, yin
        self.allow_outside = allow_outside
        if allow_outside and fillvalue is None:
            raise ValueError('Need fillvalue with allow_outside=True')
        
        xout = np.array(xout)
        yout = np.array(yout)
        assert yout.shape == xout.shape
        assert len(xout.shape) < 3 and len(yout.shape) < 3
        if len(xout.shape) == 1:
            assert len(yout.shape) == 1
            self.out_is_1d = True
        else:
            self.out_is_1d = False
        self.shape_out = xout.shape
        self.xout, self.yout = xout, yout

        xout1d = xout.ravel()
        self.ptx, mask_x = scale_to_unit(self.xin, xout1d, allow_outside)
        yout1d = yout.ravel()
        self.pty, mask_y = scale_to_unit(self.yin, yout1d, allow_outside)
        self.mask = np.logical_or(mask_x, mask_y)

        self.nx = len(xin)
        self.ny = len(yin)
        prec_allow = None, 'double', 'single'
        if not precision in prec_allow:
            raise ValueError('Only following precisions are allowed: %s' % prec_allow)
        if precision is None:
            if any(a.dtype == np.float64 for a in [xin, yin, xout, yout]):
                self.precision = 'double'
            else:
                self.precision = 'single'
        else:
            self.precision = precision
        
        self._interp_mdle = {'single':interpfsp.interp2d,
                             'double':interpf.interp2d, None:interpf.interp2d}[self.precision]
        ix_1base, iy_1base, self.weights = \
            self._interp_mdle.get_weights(self.ptx+1, self.pty+1, self.mask, self.nx, self.ny)
        #self._ix_1base = ix_0base #+ 1
        #self._iy_1base = iy_0base #+ 1
        self._ix_1base = ix_1base
        self._iy_1base = iy_1base
        
        assert np.all(np.logical_or(ix_1base >= 1, self.mask))
        assert np.all(np.logical_or(ix_1base <= self.nx, self.mask))
        assert np.all(np.logical_or(iy_1base >= 1, self.mask))
        assert np.all(np.logical_or(iy_1base <= self.ny, self.mask))

        if fillvalue is None:
            self.fillvalue = -1e30
        else:
            self.fillvalue = fillvalue
            
    def interpolate(self, data):
        """
        Perform interpolation. Data may be 2D or 3D.
        """
        if data.shape[:2] != (self.nx, self.ny):
            raise ValueError('Data has wrong shape')
        masked = hasattr(data, 'mask')
        #print 'Masked::', masked
        if len(data.shape) == 2:
            if masked:
                mask_curr = self._interp_mdle.chk_mask(self._ix_1base, self._iy_1base,
                                                       self.weights, data.mask)
                #print 'Mask curr', mask_curr
                mask = np.logical_or(self.mask, mask_curr)
            else:
                mask = self.mask
            val = self._interp_mdle.interp(self._ix_1base, self._iy_1base, self.weights, mask, data,
                                           self.fillvalue)
            if masked:
                val = np.ma.MaskedArray(val, mask)
            if not self.out_is_1d:
                val.shape = self.shape_out
            return val
        elif len(data.shape) == 3:
            return self._interpolate3d_dumb_noblk(data, masked)
            #return self._interpolate3d(data)
        else:
            raise ValueError('Data has wrong shape')

    # tried various arrangements for 3D data, the simplest is best.
    # - operating columnwise is not faster when data is in x-y-z order
    # - blocking the x-y loops is not faster/necessary (maybe with very large inputs).
        
    def _interpolate3d(self, data):
        assert data.shape[0:2] == (self.nx, self.ny)
        assert not self.out_is_1d
        blksz = 120/data.shape[2]
        #print 'dtype in', data.dtype
        val = self._interp_mdle.interp_xyz_blk(self._ix_1base, self._iy_1base, self.weights, self.mask,
                                              data, self.fillvalue, self.shape_out[0], self.shape_out[1], blksz)
        #print 'dtype out', data.dtype
        return val

    def _interpolate3d_dumb_blk(self, data):
        if self.out_is_1d:
            val_out = np.empty((self.shape_out[0], data.shape[2]), order='f', dtype=data.dtype)
        else:
            val_out = np.empty((self.shape_out[0], self.shape_out[1], data.shape[2]), order='f', dtype=data.dtype)
        blksz = 256
        for iz in xrange(data.shape[2]):
            val_out[...,iz] = self._interp_mdle.interp(self._ix_1base, self._iy_1base,
                                                       self.weights, self.mask,
                                                       data[:,:,iz], self.fillvalue).reshape(self.shape_out)
            #val_out[:,:,iz] = self._interp_mdle.interp_blk(self._ix_1base, self._iy_1base,
            #                                               self.weights, self.mask,
            #                                               data[:,:,iz], self.fillvalue,
            #                                               self.shape_out[0], self.shape_out[1], blksz)

        return val_out
    
    def _interpolate3d_dumb_noblk(self, data, masked):
        if self.out_is_1d:
            val_out = np.empty((self.shape_out[0], data.shape[2]), order='f', dtype=data.dtype)
        else:
            val_out = np.empty((self.shape_out[0]*self.shape_out[1], data.shape[2]), order='f', dtype=data.dtype)
        mask = self.mask
        if masked:
            val_out_mask = np.zeros(val_out.shape, dtype=np.bool, order='f')
        for iz in xrange(data.shape[2]):
            if masked:
                mask_curr = self._interp_mdle.chk_mask(self._ix_1base, self._iy_1base,
                                                       self.weights, data.mask[:,:,iz])
                mask = np.logical_or(self.mask, mask_curr)
            val_out[:,iz] = self._interp_mdle.interp(self._ix_1base, self._iy_1base,
                                                     self.weights, mask,
                                                     data[:,:,iz], self.fillvalue)
            if masked:
                val_out_mask[:,iz] = mask
        if masked:
            val_out = np.ma.MaskedArray(val_out, val_out_mask)
        if self.out_is_1d:
            val_out.shape = (self.shape_out[0], data.shape[2])
        else:
            val_out.shape = (self.shape_out[0], self.shape_out[1], data.shape[2])
        return val_out
        
        
            
    
