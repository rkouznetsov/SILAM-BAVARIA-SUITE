#!/usr/bin/python

my_description="""
      Handles timeseries for bunch of stsations, calculates and plots timevars

      Enjoy! (Roux)
"""

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import sys # Not needed normally

import netCDF4 as nc4
from toolbox import timeseries as ts

import os

myhost = os.uname()[1]

if myhost.startswith("nid") or  myhost.startswith("voima") or myhost.startswith("teho") or  myhost.startswith("eslogin"):
    VoimaBug=True
else:
    VoimaBug=False


class TimeVarData:
        """
                Calculate daily, weekly and monthly quartiles
        """
        weekhours = np.arange(24*7)
        hours    = np.arange(24)
        seasonmonths = [[12, 1,2], [3,4,5], [6,7,8], [9,10,11]]
        weekdays = np.arange(7) 
        months    = np.arange(12)


        def  __init__(self,times, vals, timeFactors = None, axis=None):
                # Vals -- np.array [itime, istation]
                # or 
                #  np.array [itime]
                 self.weeklyvar =  np.nan * np.empty((24*7,4,))  # time, low,median,high,mean
                 self.hourvar   =  np.nan * np.empty((24,4))
                 self.seasonHourVar = np.nan * np.empty((4,24,4))
                 self.dowvar    =  np.nan * np.empty((7,4))
                 self.monvar    =  np.nan * np.empty((12,4))
                 self.defined = False
                 if  timeFactors == None :
                        timeHours, timeDOW, timeMon = zip(*[(t.hour, t.isoweekday(), t.month) for t in times])
                        timeHours = np.array(timeHours)
                        timeDOW   = np.array(timeDOW)
                        timeMon   = np.array(timeMon)
                 else:
                        timeHours, timeDOW, timeMon = timeFactors
#                 print np.sum(np.isfinite(vals))
                 
                 if len(vals.shape) == 2:
                    (self.ntimes,self.nstations)=vals.shape
                 elif len(vals.shape) == 1:
                    self.ntimes = len(vals)
                    self.nstations = 0
                 else:
                     print "TimeVarData got wrong-shaped value argument", vals.shape
                     raise ValueError 



                 if np.sum(np.isfinite(vals)) < 200:
                        print "Warning! Not enough vals for Timevar"
                        return None

                 for day in self.weekdays:
                         for hour in self.hours:
                                 select = np.logical_and(timeHours==hour, timeDOW==day+1)
                                 if (self.nstations !=0):
                                     v=vals[select,:].flat
                                 else:
                                     v=vals[select]
                                 v=v[np.isfinite(v)]
#                                 try:
                                 if len(v) > 0 :
                                    self.weeklyvar[day*24+hour,0:3] = np.percentile(v,[25,50,75],axis=axis)
                                    self.weeklyvar[day*24+hour,3] = np.mean(v)
#                                  except:
#                                    print  vals.shape,  sum(select)
#                                    print v
 #                                   exit()

                 for day in self.weekdays:
                         select = (timeDOW==day+1)
#                         print "dow select:", np.where(select)
                         if (self.nstations !=0):
                             v=vals[select,:].flat
                         else:
                             v=vals[select]
                         v=v[np.isfinite(v)]
                         if len(v) > 0 :
                                self.dowvar[day,0:3] = np.percentile(v,[25,50,75],axis=axis)
                                self.dowvar[day,3] = np.mean(v)
                 #sys.exit()

                 for hour in self.hours:
                         select = (timeHours==hour)
                         if (self.nstations !=0):
                             v=vals[select,:].flat
                         else:
                             v=vals[select]
                         v=v[np.isfinite(v)]

                         if len(v) > 0 :
                                 self.hourvar[hour,0:3] = np.percentile(v,[25,50,75],axis=axis)
                                 self.hourvar[hour,3] = np.mean(v)

                 for iseason, months in enumerate(self.seasonmonths):
                     for hour in self.hours:
                             select = np.logical_and(timeHours==hour, (timeMon==months[0])+(timeMon==months[1])+(timeMon==months[2])) 
                             if (self.nstations !=0):
                                 v=vals[select,:].flat
                             else:
                                 v=vals[select]
                             v=v[np.isfinite(v)]
                             if len(v) > 0 :
                                     self.seasonHourVar[iseason, hour,0:3] = np.percentile(v,[25,50,75],axis=axis)
                                     self.seasonHourVar[iseason, hour,3] = np.mean(v)

                 for month in self.months:
                         select = (timeMon==month+1)
                         if (self.nstations !=0):
                             v=vals[select,:].flat
                         else:
                             v=vals[select]
                         v=v[np.isfinite(v)]
                         if len(v) > 0 :
                                 self.monvar[month,0:3] = np.percentile(v,[25,50,75],axis=axis)
                                 self.monvar[month,3] = np.mean(v)
                 
                 self.defined = True
                 return None

class EmisCorrection:
        """
                Calculates emission correction factors from ratios of medians
        """
        
        def  __init__(self,obsVar, modVar):
                self.hourcorr = obsVar.hourvar[:,1]/modVar.hourvar[:,1]
                self.dowcorr  = obsVar.dowvar[:,1]/modVar.dowvar[:,1]
                self.moncorr  = obsVar.monvar[:,1]/modVar.monvar[:,1]
                return None
        def dict(self):
                return {"hourcorr":self.hourcorr, "dowcorr":self.dowcorr, "moncorr":self.moncorr }



def PlotTimevars(fig, timevars,  labels=None, title=None, units=None, plotNeg=False):
        """
                plots array of timevars (array of timevar instances)
        """
        axweek=fig.add_subplot("311")
        axweek.set_xlim([0,24*7])
        axweek.set_xticks(range(0,24*7,24))
        axweek.xaxis.set_major_formatter(ticker.NullFormatter())
        axweek.xaxis.set_minor_locator(ticker.FixedLocator(np.arange(7)*24+12))
        axweek.xaxis.set_minor_formatter(ticker.FixedFormatter('Mon Tue Wed Thu Fri Sat Sun'.split()))

        seasonnames = "djf mam jja son".split()
        if title:
                axweek.set_title(title)

        axhour = fig.add_subplot("345")
        axhour.set_xlim([0,24])
        axhour.set_xticks(range(0,24,6))
        axdow = fig.add_subplot("346")
        axdow.set_xlim([0,7])
        axdow.xaxis.set_major_formatter(ticker.NullFormatter())
        axdow.xaxis.set_minor_locator(ticker.FixedLocator(np.arange(7)+0.5))
        axdow.xaxis.set_minor_formatter(ticker.FixedFormatter('Mo Tu We Th Fr Sa Su'.split()))
        axmon = fig.add_subplot("347")
        axmon.set_xlim([0,12])
        axmon.xaxis.set_major_locator(ticker.FixedLocator(np.arange(11)+1))
        axmon.xaxis.set_major_formatter(ticker.NullFormatter())
        axmon.xaxis.set_minor_locator(ticker.FixedLocator(np.arange(12)+0.5))
        axmon.xaxis.set_minor_formatter(ticker.FixedFormatter('J F M A M J J A S O N D'.split()))
        axlegend = fig.add_subplot("348")
        axes = [axweek, axhour, axdow, axmon]
        for iseason,season in enumerate(seasonnames):
            axes.append(fig.add_subplot(3,4,9+iseason))
            axes[-1].set_title(season)
            axes[-1].set_xlim([0,24])
            axes[-1].set_xticks(range(0,24,6))

        colors="red blue green cyan  magenta black".split()
        maxy = max([np.nanmax(v.weeklyvar[:,:]) for v in filter(None, timevars) ])
	
	
        if plotNeg:
            miny = min([np.nanmin(v.weeklyvar[:,:]) for v in filter(None, timevars) ])
            for ax in axes:
                ax.plot([0,150],[0,0], color='k', ls='-') # Zero line
        else:
          if not maxy > 0 :
              print "Maxy is wrong ", maxy
              raise ValueError("Maxy is wrong for strictly positive quantity") 
          miny = 0


        nvars=len(timevars)
        if nvars > len(colors):
          raise ValueError("Too many variables for Timevar")

        if labels==None:
                labels = ["%d"%i for i in range(nvars)]
        for ivar in range(nvars):
                var=timevars[ivar]
                if not var:
                    continue
                myxvals = [var.weekhours,var.hours, var.weekdays, var.months, var.hours, var.hours, var.hours, var.hours,]  
                myvars  = [var.weeklyvar, var.hourvar, var.dowvar, var.monvar,  var.seasonHourVar[0, :,:],  var.seasonHourVar[1, :,:],  var.seasonHourVar[2, :,:],  var.seasonHourVar[3, :,:] ]
                for ax,x,v in zip(axes,myxvals,myvars): 
                    ax.plot(x+0.5,v[:,1], color=colors[ivar], label=labels[ivar]) #Median
                    ax.plot(x+0.5,v[:,3], color=colors[ivar],  ls = '--', dashes=(5, 2)) #mean
                    ax.fill_between(x+0.5,v[:,0], v[:,2], color=colors[ivar], alpha=0.2)
                    ax.set_ylim([miny,maxy])
                    ax.set_label(units)
                    ax.grid()

        ## Legend axis
        axlegend.set_xlim([-2,-1])
        axlegend.axis('off')
        for ivar in range(nvars):
            axlegend.plot(range(2), range(2), color=colors[ivar], label=labels[ivar])
        axlegend.plot(range(2), range(2), color='grey',  ls = '-', label="Median")
        axlegend.plot(range(2), range(2), color='grey',  ls = '--', dashes=(5, 2), label="Mean")
        axlegend.fill_between(range(2), range(2), color='grey',  alpha=0.2,  label="25th-75th prcntl")
        axlegend.legend( loc='center', frameon=True)                    
        fig.tight_layout()


class SpatialStats:
    """
       Contain Spatial statistics 
    """
    def __init__(self, times, valsobs, valsmod, flipBias = False):
#        print np.sum(np.isfinite(vals))
        (self.ntimes,self.nstations)=valsobs.shape
        if (valsmod.shape != valsobs.shape):
            print valsmod.shape, valsobs.shape
            print "SpatialStats got different-shape arrays"
            raise ValueError

        self.times = times
        self.bias = np.nanmean(valsmod - valsobs, axis=1)
        self.FracB = 2*np.nanmean((valsmod - valsobs)/(valsmod + valsobs + 0.1), axis=1)
        self.FGerr = 2*np.nanmean(abs(valsmod - valsobs)/(valsmod + valsobs + 0.1), axis=1)

        if flipBias:
            self.bias *= -1
            self.FracB *= -1
        self.RMS  = np.sqrt(np.nanmean((valsmod - valsobs)*(valsmod - valsobs), axis=1)) 

        self.N = np.nan * np.empty((self.ntimes))
        self.corr = np.nan * np.empty((self.ntimes))
        self.meanOBS = np.nan * np.empty((self.ntimes))
        self.meanMOD = np.nan * np.empty((self.ntimes))
        self.stdOBS = np.nan * np.empty((self.ntimes))
        self.stdMOD = np.nan * np.empty((self.ntimes))

        for it in range(self.ntimes):
            idx=np.isfinite(valsmod[it,:] - valsobs[it,:])
            N = np.sum(idx)
            self.N[it] = N
            if (N<5): continue


            modmean = np.mean(valsmod[it, idx]) 
            obsmean = np.mean(valsobs[it, idx]) 
            self.meanOBS[it] = obsmean
            self.meanMOD[it] = modmean

            self.stdOBS[it] = np.sqrt(np.mean( (valsobs[it, idx] - obsmean)**2 ))
            self.stdMOD[it] = np.sqrt(np.mean( (valsmod[it, idx] - modmean)**2 ))

            self.corr[it] = np.mean(( valsmod[it, idx] - modmean )*(valsobs[it, idx] - obsmean))
            self.corr[it] /= self.stdOBS[it]*self.stdMOD[it]


    def getTV(self,stat):
        if stat == "RMSE":
            return TimeVarData(self.times,self.RMS)
        elif stat == "FGerr":
            return TimeVarData(self.times,self.FGerr)
        elif stat == "Corr":
            return TimeVarData(self.times,self.corr)
        elif stat == "Bias":
            return TimeVarData(self.times,self.bias)
        elif stat == "FracB":
            return TimeVarData(self.times,self.FracB)
        elif stat == "N":
            return TimeVarData(self.times,self.N)
        elif stat == "meanOBS":
            return TimeVarData(self.times,self.meanOBS)
        elif stat == "meanMOD":
            return TimeVarData(self.times,self.meanMOD)
        elif stat == "stdOBS":
            return TimeVarData(self.times,self.stdOBS)
        elif stat == "stdMOD":
            return TimeVarData(self.times,self.stdMOD)
        else:
            print "timevar can be one of RMSE Corr Bias"
            raise ValueError


class TemporalStats:

    """
       Contain Spatial statistics 
    """
    def __init__(self, lons, lats, case, valsobs, valsmod, flipBias = False):

#        print np.sum(np.isfinite(vals))
        (self.ntimes,self.nstations)=valsobs.shape
        if (valsmod.shape != valsobs.shape):
            print valsmod.shape, valsobs.shape
            print "TemporalStats got different-shape arrays"
            raise ValueError

        assert (len(lons) == self.nstations)
        assert (len(lats) == self.nstations)
        self.lons  = lons
        self.lats  = lats
        self.case  = case ##To appear in the table header and in the title

        self.validstats="Bias RMS N Corr meanObs meanMod stdObs stdMod MORatio MOstdRatio FracB FGerr".split()

        self.Bias = np.nanmean(valsmod - valsobs, axis=0)
        self.RMS  = np.sqrt(np.nanmean((valsmod - valsobs)*(valsmod - valsobs), axis=0)) 
        self.FracB = 2*np.nanmean((valsmod - valsobs)/(valsmod + valsobs + 0.1), axis=0)
        self.FGerr = 2*np.nanmean(abs(valsmod - valsobs)/(valsmod + valsobs + 0.1), axis=0)

        self.N = np.nan * np.empty((self.nstations))
        self.Corr = np.nan * np.empty((self.nstations))
        self.meanObs = np.nan * np.empty((self.nstations))
        self.meanMod = np.nan * np.empty((self.nstations))
        self.stdObs = np.nan * np.empty((self.nstations))
        self.stdMod = np.nan * np.empty((self.nstations))

        for ist in range(self.nstations):
            idx=np.isfinite(valsmod[:,ist] - valsobs[:,ist])
            N = np.sum(idx)
            self.N[ist] = N
            if (N<5): continue
            modmean = np.mean(valsmod[idx, ist]) 
            obsmean = np.mean(valsobs[idx, ist]) 
            self.meanObs[ist] = obsmean
            self.meanMod[ist] = modmean

            self.stdObs[ist] = np.sqrt(np.mean( (valsobs[idx,ist] - obsmean)**2 ))
            self.stdMod[ist] = np.sqrt(np.mean( (valsmod[idx,ist] - modmean)**2 ))

            self.Corr[ist] = np.mean(( valsmod[idx,ist] - modmean )*(valsobs[idx,ist] - obsmean))
            self.Corr[ist] /= self.stdObs[ist]*self.stdMod[ist]

        self.MORatio    = self.meanMod/self.meanObs
        self.MOstdRatio = self.stdMod/self.stdObs

    def getstat(self, statname):
        if statname in self.validstats:
            return getattr(self, statname)
        else:
            print "getstat alled with statname = '%s'"%(statname)
            raise ValueError("statlist can be one of "+",".merge(self.validstats))



    def printStats(self, outf, statlist,  statproc, rawHead, ifHead, idxStations=None):

        #case -- string 
        # statlist list of strings of statistics
        if ifHead:
            outf.write("\n\n%s\n"%(self.case,))
            outf.write("%25s"%("",))
            for st in statlist:
                outf.write("%10s"%(st,))
            outf.write("\n")

        statfuns=dict(
                    mean = np.nanmean, 
                    median = np.nanmedian, 
                    prc95 = lambda x : np.nanpercentile(x,95),
                    first = lambda x : x[0],  ##First station from the list
                    )

        try:
            stfun = statfuns[statproc]
        except KeyError:
            raise KeyError("statproc argument can be one of "+",".merge(statfuns.keys()))

        Nstations=self.nstations
        if idxStations:
            Nstations = np.sum(np.isfinite(self.meanObs[idxStations]))

        rawtitle="%25s"%("%s %s N=%d"%(rawHead,statproc,Nstations))
        outf.write(rawtitle)
        for st in statlist:
            statvals = self.getstat(st)
            if idxStations:
                outf.write("%10.2f"%(stfun(statvals[idxStations])))
            else:
                outf.write("%10.2f"%(stfun(statvals)))
        outf.write("\n")

    def plot2Map(self, basemap,  stat, cmap, norm, idxStations):
        starr=self.getstat(stat)
        return basemap.scatter(self.lons[idxStations], self.lats[idxStations], 
                    c=starr[idxStations], s=30, cmap=cmap, norm=norm, linewidths=0)
        






class SeasonalStats:

    """
       Contain Temporal statistics over Seasons
    """
    def __init__(self, obsMatr, modMatr, case):
        stations = obsMatr.stations
        hrlist  =  obsMatr.times

        valsobs = obsMatr.vals
        valsmod = modMatr.vals

        lons = np.array([st.x for st in stations])
        lats = np.array([st.y for st in stations])



#        print np.sum(np.isfinite(vals))
        (self.ntimes,self.nstations)=valsobs.shape
        if (valsmod.shape != valsobs.shape):
            print valsmod.shape, valsobs.shape
            print "TemporalStats got different-shape arrays"
            raise ValueError


        self.seasons="all djf mam jja son".split()
        self.stats={}
        seasonmonths = dict(djf=[12, 1,2], mam=[3,4,5], jja=[6,7,8], son=[9,10,11], all=range(1,13))

        self.stats["all"] = TemporalStats(lons, lats, case, valsobs, valsmod)
        timeMon = np.array([t.month for t in hrlist])
        for season in "djf mam jja son".split():
            months = seasonmonths[season]
            idxtime = (timeMon==months[0])+(timeMon==months[1])+(timeMon==months[2])
            self.stats[season] = TemporalStats(lons, lats, case+", "+season, valsobs[idxtime,:], valsmod[idxtime,:])









#
# Read the data
#
class PollutantModobs:
        """
                Class to contain synchronous data, station info, time axis etc.
                with finctions to select stations by ID etc
        """
        def __init__(self, ncpol, ncnames, titles=None, airbasefile=""):
                """
                        ncnames -- list of filenames, 
                """

                
                if titles==None:
                        titles=ncnames

                # Read the data
                self.titles=[]
                self.data={}
                stnames = None
                self.ncpol=ncpol # Variable name in .nc files
                for infile, title in zip(ncnames,titles):
                       if stnames==None: # First file (presumably, observations)
                               self.titles.append(title)
                               nc=netcdf.netcdf_file(infile,"r")
                               stnames = np.array([''.join(nc.variables['stn_id'].data[ist,:]) for ist in  range(nc.dimensions["station"])])
                               mdlTime = netcdftime.utime(nc.variables['time'].units)
                               self.times = mdlTime.num2date(nc.variables['time'].data) 
                               #for i in range(1000):
                                #        print self.times[i].strftime("%a %F %T")
                                
                                # Valid satation indices
                               nodataidx = [ np.all(np.isnan(nc.variables[ncpol].data[:,i]))   for i in range(len(stnames)) ]
                               self.stselect = np.logical_not(np.array(nodataidx))
                               self.stid=np.array(stnames[self.stselect]) # station ID
                               self.units=nc.variables[ncpol].units
                               self.data[title]=nc.variables[ncpol].data[:,self.stselect]
                               nc.close()
                               
                       else:
                               self.AddData(infile,title)
                
                
                airbasetypes = getAirbaseTypes(airbasefile)

                stidx = np.array( [airbasetypes["idx"][name] 
                                        for name in self.stid]
                                        )
                self.sttype = np.array(["%s_%s"%(airbasetypes["station_type_of_area"][i],
                                             airbasetypes["type_of_station"][i].lower(),
                                            ) for i in stidx ]
                                  )
                self.countrycode = np.array([airbasetypes["country_iso_code"][i] for i in stidx])
                self.country = np.array([airbasetypes["country_name"][i] for i in stidx])
                #for i in range(10):
                #       print "%10s %10s %5.2f"%(self.stid[i], self.country[i],nc.variables[ncpol].data[5,i])
                
                timeHours, timeDOW, timeMon = zip(*[(t.hour, t.isoweekday(), t.month) for t in self.times])
                self.timeHours=np.array(timeHours)
                self.timeDOW=np.array(timeDOW)
                self.timeMon=np.array(timeMon)

        def AddData(self,ncname,title): # adds one more ncfile to the existing dataset
               self.titles.append(title)
               nc=netcdf.netcdf_file(ncname,"r")
               mdlTime = netcdftime.utime(nc.variables['time'].units)
               mytimes = mdlTime.num2date(nc.variables['time'].data)
               if (mytimes != self.times).any():
                        print "Error: incompatible times in %s"%(ncname)
                        print "OBS", self.times[0].strftime("%F %T"), self.times[-1].strftime("%F %T")
                        print "mod", mytimes[0].strftime("%F %T"), mytimes[-1].strftime("%F %T")
                        exit (-1)

               mystnames = np.array([''.join(nc.variables['stn_id'].data[ist,:]) for ist in  range(nc.dimensions["station"])])
               mystid = mystnames[self.stselect]
               if (np.any(mystid != self.stid)):
                        print "Error: incompatible stations"
                        print mystid
                        print stid
                        exit(-1)
               self.data[title]=nc.variables[self.ncpol].data[:,self.stselect]
               nc.close()

        def DelData(self,title): # adds one more ncfile to the existing dataset
               self.titles.remove(title)
               del self.data[title]

        
        def SelectStations(self,names=None, country=None, ccode=None, sttype=None, minstations = 1):
                """
                        return a bool array  selecting  stations
                """
                sel=np.ones((len(self.stid)),dtype=bool)
                if names != None:
                        sel *= np.in1d(self.stid, names)
                
                if ccode != None:
                        #sel *= np.array([c in ccode for c in   self.countrycode])
                        sel *= np.in1d(self.countrycode,ccode)

                if country != None:
#                        sel *= np.array([c in country for c in  self.country])
                        sel *= np.in1d(self.country, country)

                if sttype != None:
                        #sel *= np.array([t in sttype for t in  self.sttype])
                        sel *= np.in1d(self.sttype, sttype)

                if np.sum(sel) >= minstations:
                        return sel
                else:
                        return None



class TsMatrix:

        def __init__(self, times, stations, vals, units):        
    #        print np.sum(np.isfinite(vals))
            (self.ntimes,self.nstations)=vals.shape

            if ((len(times),len(stations)) != vals.shape):
                print (len(times),len(stations)), "!=", vals.shape
                print "TsMatrix got incompatibel-shape arrays"
                raise ValueError

            self.times = times
            self.stations = stations
            self.vals = vals
            self.units = units

        def to_nc(self, filename):

            nt,nst = len(self.times),len(self.stations)
            reftime=self.times[0]
            
            with nc4.Dataset(filename, "w", format="NETCDF4") as outf:

                outf.featureType = "timeSeries";
                strlen=64

                time = outf.createDimension("time", nt)
                slen = outf.createDimension("name_strlen", strlen)
                station    = outf.createDimension("station",nst)

                t = outf.createVariable("time","i4",("time",))
                t.standard_name="time"
                t.long_name="end of averaging interval"
                t.calendar="standard"
                t.units = reftime.strftime("minutes since %Y-%m-%d %H:%M:%S")
                t[:] = [ (h - reftime).total_seconds()/60 for h in self.times ]

                lon = outf.createVariable("lon","f4",("station",))
                lon.standard_name = "longitude"
                lon.long_name = "station longitude"
                lon.units = "degrees_east"

                lat = outf.createVariable("lat","f4",("station",))
                lat.standard_name = "latitude"
                lat.latg_name = "station latitude"
                lat.units = "degrees_north"

                alt= outf.createVariable("alt","f4",("station",))
                alt.standard_name = "surface_altitude"
                alt.latg_name = "station altitude asl"
                alt.units = "m"
                alt.positive = "up";
                alt.axis = "Z";

                stcode = outf.createVariable("station_code","c",("station","name_strlen"))
                stcode.long_name = "station code"
                stcode.cf_role = "timeseries_id";
                
                stname = outf.createVariable("station_name","c",("station","name_strlen"))
                stname.long_name = "station name"
                
                starea = outf.createVariable("area_type","c",("station","name_strlen"))
                starea.long_name = "area_type"
                
                stsource = outf.createVariable("source_type","c",("station","name_strlen"))
                stsource.long_name = "source_type"

                val= outf.createVariable("val","f4",("time","station"), zlib=True, complevel=4)
                val.coordinates = "lat lon alt station_code"
                val.units = self.units

                ## .01 absolute precision decimals=2
               # if not (decimals is None):
               #     outdat = np.around(self.vals,  decimals=decimals)
                outdat = self.vals
                
                #Output sorted stations
                for ist,st in enumerate(self.stations):

                    stcode[ist,:] = nc4.stringtoarr( st.code, strlen,dtype='S')
                    stname[ist,:] = nc4.stringtoarr( st.name, strlen,dtype='S')
                    lon[ist] = st.x
                    lat[ist] = st.y
                    alt[ist] = st.hgt
                    starea[ist,:] = nc4.stringtoarr( st.area, strlen,dtype='S')
                    stsource[ist,:] = nc4.stringtoarr( st.source, strlen,dtype='S')

                precision=500 ##Relative precision
                keepbits=np.int(np.ceil(np.log2(precision)))
                maskbits=20 -keepbits
                mask=(0xFFFFFFFF >> maskbits)<<maskbits
                b=outdat.view(dtype=np.int32)
                b &= mask

                val[:] = outdat

        @classmethod
        def fromNC(self,ncfile):

          with nc4.Dataset(ncfile) as nc:
              stlist = nc4.chartostring(nc.variables['station_code'][:])
              tvals = nc.variables['time'][:]
              tunits = nc.variables['time'].units
              #
              ## in Voima there is an issue:
              #https://jira.fmi.fi/browse/STU-10345

              if VoimaBug:
                  hrlist = np.array([ nc4.num2date(t,tunits) for t in tvals])
              else:
                  hrlist = nc4.num2date(tvals,tunits)


              valmatr = nc.variables['val'][:]
              try:
                  units = nc.variables['val'].units
              except AttributeError:
                  print "Unknown units in %s, ug/m3 assumed"%(ncfile)
                  units="ug/m3"



              #Get stations
              stnames = nc4.chartostring(nc.variables['station_name'][:])
              stcodes = nc4.chartostring(nc.variables['station_code'][:])
              stareas = nc4.chartostring(nc.variables['area_type'][:])
              stsources = nc4.chartostring(nc.variables['source_type'][:])
              stx = nc.variables['lon'][:]
              sty = nc.variables['lat'][:]
              stz = nc.variables['alt'][:]

              nst = len(stlist)
              stations = []
              for i in range(nst):
                  stations.append(ts.Station(stcodes[i], stnames[i], stx[i], sty[i], 
                          height=stz[i], area_type=stareas[i], dominant_source=stsources[i]))

          return self(hrlist,stations,valmatr, units)




