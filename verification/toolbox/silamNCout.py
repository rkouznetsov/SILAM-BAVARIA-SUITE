#
# Wrapper for pupyrene to handle temporary files
# "append" to existing netcdfs 
# easy create variables and dimensions filled with data and attributes
#
#   Roux



import os
import sys
#sys.path.insert(0, '/lustre/apps/silam/tools/python/')
import numpy as np




try:
        import pupynere as netcdf
except:
        from support import pupynere as netcdf



class silamNCout(netcdf.netcdf_file): # Class for our netcdf files...
#        targetfile=None
#        verbose=False

        # Create temporary file (not to destroy the prevous
        # version before writing succeeded) and copy there the
        # stuff from the previous file: a way to
        # append variables to existing file
        def __init__(self, filename, append=True, verbose=False, allowrecord=True):
                (head,tail) = os.path.split(filename)
                tmpfile=os.path.join(head,"0-%d-%s"%(os.getpid(),tail))
                netcdf.netcdf_file.__init__(self, tmpfile, "w")

                
                self.__dict__['_targetfile']=filename #Setattr writes also to self._attributes
                self.__dict__['_verbose']=verbose     # We do not want it...
                self.__dict__['_records']=allowrecord     # No record dimension
                oldF=None
                if append:
                      try: 
                        oldF = netcdf.netcdf_file(filename,"r")
                      except IOError as e:
                         if verbose:
                             print "Previous file can't be updated..."

                # Copy from the oldfile
                if oldF: # Copy all stuff from the old file

                        for k, v in oldF._attributes.items():
                                if hasattr(self, k): #Overwriting targetfile is harmful
                                                     # it should not appear in files anyhow
                                        continue
                                self.__setattr__(k,v)
                        for dim in oldF._dims: 
                                if (allowrecord or oldF.dimensions[dim]>0):
                                        value = oldF.dimensions[dim]
                                else:
                                        value = oldF.variables[dim].shape[0]
                                if value>-1:
                                        self.createDimension(dim, int(value))
                                else:
                                        self.createDimension(dim, None)
                        for key, var in oldF.variables.items():
                                 ovar = self.createVariable(key, var.typecode(), var.dimensions)
                                 for k, v in var._attributes.items():
                                         ovar.__setattr__(k,v)
                                 if len(var.data.shape):
                                        ovar[:]=var.data.copy()
                                 else:
                                        ovar=var.data.copy()
                        oldF.close()
                else:
                        self.history=""

        # 
        #
        def finalize(self):
                if self._verbose:
                        print "Closing tmpfile", self.filename
                self.close()
                if self._targetfile:
                    if self._verbose:
                          print "Renaming %s to %s"%(self.filename,self._targetfile)
                    os.rename(self.filename,self._targetfile)
        
        # Checks if var exists and creates or verifies it
        # intended to be used with dimensions 
        # They have to be same and have same values
        # Time is unlimited dimension
        def createVar(self, varname, dtype, dims, vals, attrs={}, epsilon = 0):
                if (varname,)==dims: # The variable is dimension itself
                        if (varname=='time') and self._records:
                                dimlen = None
                        else: 
                                dimlen = len(vals)
                        if  varname in self.dimensions:
                                if self.dimensions[varname] != dimlen:
                                        print "Dimension %s has size"%(varname,), self.dimensions[varname], "requested", dimelen
                                        self.close()
                                        exit(-1)
                        else:
                                self.createDimension(varname, dimlen)
                # values
                if  varname in self.variables:
                        var=self.variables[varname]
                        if var.dimensions != dims:
                                print "Wrong datatype for", varname, var.dimensions, ", requested", dims
                                self.close()
                                exit(-1)
                        if var.data.shape != vals.shape:
                                print "Wrong shape of", varname, var.data.shape, ", requested", vals.shape
                                self.close()
                                exit(-1)
                        
                        if dtype != 'c':
                            VallErr=np.any(abs(var.data - vals)>epsilon)
                        else:
                            VallErr= (var.data != np.array(vals))  #False  #FIXME Hack. Do not know how to handle char scalars
                            VallErr= False  #FIXME Hack. Do not know how to handle char scalars
                         
                        if VallErr:
                                print "Wrong data", varname,"in file",var.data, ", requested", vals
                                self.close()
                                exit(-1)


                        if ("units" in attrs) and (var._attributes['units'] != attrs['units']):
                                print "Wrong units for ", varname, var.units, ", requested", attrs['units']
                                self.close()
                                exit(-1)

                else:
                    var = self.createVariable(varname, dtype, dims)
                    if (vals != None): 
                       # Trick to handle scalar variables...
                        if vals.shape:
                                var[:] = vals
                        else:
                                var.assignValue(np.array(vals))

                # Set attributes
                for k,v in attrs.items():
                                var.__setattr__(k,v)

        
