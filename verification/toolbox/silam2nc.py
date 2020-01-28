#!/usr/bin/env python

my_description="""
  Converts silam grads into netCDF trying to produce reasonable attributes
  understandable by thredds.
  Some semi-intelligent guess of units and other attributes is done.
  Also tries to bring some sense into variable naming convention
  (to avoid cnc_38 variables), understands rotated-pole grids.
 
  
  Can make temporal slicing, spatial subsetting, selection of variables, 
  split into separate files taking surface values only (with optional 
  keeping vertical dimension for fields that have it), 
  calculation of total PM, appending to earlier written files
 
  For now (09.2012) SILAM does not provide any reasonable way
  to tell from output files (1)if a quantity is cumulative or not,
  and (2)units. Thus, non-cumulative are assumed for everything.
 
  If the same variable is requested from several input files -- last wins.
 
  Nucleides are still TODO...
  Enjoy! (Roux)
"""


import datetime as dt
import time
import threading
import sys

import numpy as np
from toolbox import namelist
from toolbox import silamfile as sf
from toolbox import silamNCout
import projections as proj
import os
import re
import argparse

try:
    from support import ncepgrib2
    Have_ncepgrib2 = True
except:
    Have_ncepgrib2 = False
    
    


#default location for chemical files
chem_db='silam_chemicals.dat'
os.environ["TZ"]="UTC" # Should be UTC unless you want something weird...
totpmspecies="SO4,NH415SO4,NH4NO3,sslt,PM,PM_FRP" # Default species for total PM
#molar_mass = {'SO4':0.096, 'NH415SO4':0.123, 'NH4NO3':0.080, 'sslt':1.0, 'PM':1.0}

# Default attributes
epoch=dt.datetime(2000, 1, 1)                

def addContact(outF,attrlist):
        outF.source = "SILAM XXX"
        outF.institution= "FMI, Finnish Meteorological Institute"
        outF.title = "Air quality forecast" ;
        outF.Contact = "Mikhail Sofiev"
        outF.Email   = "firstname.lastname@fmi.fi"
        outF.Phone   = "+358 29 539 1"
        if attrlist:
          for a in attrlist.split(","): #FIXME Should allow for escaped commas
                o=re.match("([^=]+)=([^=]+)",a)
                if o:
                        val=re.sub("_"," ",o.group(2))
                        outF._attributes[o.group(1)]=val
        return


#class KeyboardInterruptError(Exception): pass

def ProcessCmdLine(): # returns hash with processed arguments

        # Parser for per-file arguments
        file_parser=argparse.ArgumentParser(add_help=False)
        file_parser.add_argument("-v", "--only-vars", help="coma-separated list of variables (regexps) to output. \
                        All vars are in output by default")
        file_parser.add_argument("-d", "--drop-vars", help="coma-separated list of variables (regexp) to skip in output list")
        
        pmargs=file_parser.add_mutually_exclusive_group()
        pmargs.add_argument("--pm-species",  help="coma-separated list of SPECIES (regexp) to include in PM\
                        leave empty to select automaticaly default: %s; non-pm species are ignored"%totpmspecies, default=totpmspecies) 
        pmargs.add_argument("--pm-vars", help="coma-separated list of VARIABLES (regexp), to include in PM", default=None)
        pmargs.add_argument("--tag", help="\"tag\" attribute to be added to each variable taken form this file", default=None)
        file_parser.add_argument("--skip-pm-vars", help="coma-separated list of VARIABLES (regexp), to Exclude from \
                                        selected by --pm-species or --pm-vars", default=None)


        file_parser.add_argument("super_ctl", help="SILAM super ctl file")
        
        # This should go after the last argument added to file_parser
        fpusage = " ".join((file_parser.format_usage().split())[2:]) # Fileparser usage
        fphelp = "\n".join((file_parser.format_help().splitlines())[2:]) # Remove Usage from help
        fphelp = re.sub("optional arguments:","per-file optional arguments:",fphelp)

        parser = argparse.ArgumentParser(
                formatter_class=argparse.RawDescriptionHelpFormatter,
                description=my_description,
                epilog=fphelp
                )
        
        srfargs=parser.add_mutually_exclusive_group()
#        srfargs.add_argument("-S", "--surface-keep-level",help="Surface data only, dummy height dimension",
#                                            action="store_true")
        srfargs.add_argument("-s", "--surface-only",help="Surface data only, no height dimension",
                                            action="store_true")
        
        parser.add_argument("-o", "--output-nc", 
                        help="name of the netcdf output, accepts strftime templates,  \
                           By default creates file in the current directory taking the name from .super_ctl")

        parser.add_argument("-g", "--output-grib", 
                        help="name of the grib output, accepts strftime templates,  \
                           creates file in the current directory taking the name from .super_ctl")

        parser.add_argument("--prod-status", help="Production status of data (GRIB2) \
                0: Operational Products,  1: Operational Test Products, 2: Research Products",  type=int, default=2)

        parser.add_argument("--no-record", 
                        help="Make time dimension non-record in NetCDF", 
                        action="store_true")

        parser.add_argument("--absolute-time", 
                        help="Store time as float YYYYmmdd.dayfraction  in NetCDF(not CF compatible, but readable by cdo)", 
                        action="store_true")

        parser.add_argument("--epoch", 
                        help="Give epoch for NetCDF file as \"YYYY-mm-dd\" or \"start\" for the beginning of first output", 
                        default="1970-01-01")
        
        parser.add_argument("--time-units", 
                        help="NetCDF time-units \"s[econds]\" \"h[ours]\" or \"d[ays]\"", 
                        default="seconds")
        
        bboxargs=parser.add_mutually_exclusive_group()
        bboxargs.add_argument("-b", "--bbox", help="BROKEN! Select bounding box \"rlomin,rlomax,rlamin,rlamax\"  \
                                in degrees (of rotated grid)", default=None)
        bboxargs.add_argument("-B", "--Bbox", 
                        help="Select bounding box \"lomin,lomax,lamin,lamax\" in degrees (of geographic grid)", 
                        default=None)

        parser.add_argument("-l", "--latlondims", help='Use "latitude" and "longitude" (or "rlat" and "rlon") as dimension names \
                                                        instead of "y" and "x"', action="store_true")
        parser.add_argument("-L", "--level", help='Use "level" as vertical dimension name \
                                                        instead of "height"', action="store_true")
        
        parser.add_argument("-P","--total-pm", help="Make PM2.5 and PM10", action="store_true")
        parser.add_argument("-n","--num-threads", help="Max number of threads to use \
                        (works with multiple output files)",  type=int, default=0)
        
        parser.add_argument("-t", "--times", help="Select time slice \"itmin[,itmax[,itstep]]\" in time indices", 
                        default=None)
        parser.add_argument("--latlon", help="Store geographical grid as 2D array",action="store_true")
        parser.add_argument("-c", "--chem-db", help="Use chemical database instead of default", default=False)
        parser.add_argument("-k", "--kilograms", help="Use kg instead of micrograms for cnc_, wd_, and dd_",action="store_true")
        parser.add_argument("-V", "--verbose",help="Be talkative", action="store_true")
        parser.add_argument("-N", "--steps-per-file",help="Split output by N time steps, out template  \
                                        (if any) is better to include \"%%d\"", type=int, default=0)
        parser.add_argument("-D", "--date",help="Use first date of current output file in its name \
                                                instead of the runs first date", action="store_true")
        parser.add_argument("-q", "--quiet",help="Be silent", action="store_true")
        parser.add_argument("-a", "--attr",help="Additional NetCDF global string attributes as par=val,par1=val1,...\
                                  no comas allowd yet in new attributes, underscores are replaced with spaces...")
        parser.add_argument("-A", "--append",help="Append variables to the existing netcdf. Checks for the dimensions (size and values only)",
                                        action="store_true")

        # To add standard_name to variables. standard_name is used by ncwms to chose default scaling...
        parser.add_argument("-p", "--pasodoble",help="Add standard name to pasodoble variables",action="store_true")

        # This should go after the last argument added to parser 
        usage = " ".join((parser.format_usage().split())[1:]) # Automatic usage
        parser.usage="%s\n        %s\n       [%s] \n        ..."%(usage, fpusage, fpusage) # Custom usage

        # Parse arguments
        args, unk = parser.parse_known_args()

        fargs=[] #per-file arguments
        n=len(unk)
        while True:
                i=0 
                for i in  range(n): # attempt parsing from start
                    Parsed, unkf = file_parser.parse_known_args(unk[i:n])
                    if not unkf: break #first fully parsed tail
                if i==n: # did not found parsable tail
                    print parser.format_help()
                    exit(-1) # Exit with error 
                else:
                    fargs.append(Parsed)
                    n=i
                if (n<=0): break
        fargs.reverse()     

        # process arguments
        if args.verbose:
                print "Arguments global:", args, "\n and per-file:"
                for argl in fargs:
                        print argl

        if args.output_grib and not Have_ncepgrib2:
            print "Grib2 output requested, but no ncepgrib2 module load...."
            exit(255)

    
        # By default output netcdf only -- generate template
        if not args.output_nc and not  args.output_grib:
                print fargs[0].super_ctl
                args.output_nc=os.path.basename(fargs[0].super_ctl)
                args.output_nc=re.sub("\.super_ctl$","",args.output_nc)
                if args.steps_per_file:
                        if args.date:
                           args.output_nc += "_%Y%m%d" 
                        else:     
                           args.output_nc += "_%%03d" 
                args.output_nc += ".nc"
        
        # Check if we have proper template
        if  args.steps_per_file:
            if args.output_nc and not  re.search("%%[0-9]*d",args.output_nc): # missing timestep
                if not args.quiet:
                        print "Warning! No room for timestep in the NetCDF output template:", args.output_nc
            if args.output_grib and not  re.search("%%[0-9]*d",args.output_grib): # missing timestep
                if not args.quiet:
                        print "Warning! No room for timestep in the GRIB2 output template:", args.output_grib

        setattr(args, "Write3D", not(args.surface_only))
        Write3D=not(args.surface_only)

        if not args.chem_db:
                if os.path.isfile(chem_db):
                        if args.verbose:
                                print "Using ", chem_db
                        setattr(args,"chem_db",chem_db)
                elif os.path.isfile("./ini/"+chem_db):
                        if args.verbose:
                                print "Using ", "./ini/"+chem_db
                        setattr(args,"chem_db","./ini/"+chem_db)
                else:
                        print "Can't find chemical database in %s or %s..."%(chem_db,"./ini/"+chem_db)
                        exit()

        setattr(args,"micrograms",not(args.kilograms))
       
        return args, fargs
        
        



# Thredds (ncwms) can choose palette basing on standard_name tag. Include it for pasodoble vars
standard_names={"cnc_co":"mass_concentration_co_in_air", 
                "cnc_no":"mass_concentration_no_in_air", 
                "cnc_no2":"mass_concentration_no2_in_air", 
                "cnc_o3":"mass_concentration_o3_in_air", 
                "cnc_pm2_5":"mass_concentration_pm25_in_air", 
                "cnc_pm10":"mass_concentration_pm10_in_air", 
                "cnc_so2":"mass_concentration_so2_in_air", 
                "cnc_co":"mass_concentration_co_in_air",
                "cnc_passive":"mass_concentration_passive_in_air",
               }

Grib2params={
#(discipline, category, number, vertical type, vertical value)
                "temperature":(0,0,0,None,None),
                "temp_2m":(0,0,0,103,2.0),
                "dew_p_temp":(0,0,6,None,None),
                "dew_p_temp2m":(0,0,6,103,2.0),
                "temp_2m_acc":(0,0,0,103,2.0),
                "U_wind_10m":(0,2,2,103,10.0),
                "V_wind_10m":(0,2,3,103,10.0),
                "windspeed_10m":(0,2,1,103,10.0),
                "pot_temp":(0,0,2,None,None),
                "cloud_cover":(0,6,22,None,None), # Fixme  -- to %
                "total_cloud":(0,6,1,1,None), # Fixme  -- to %
                "nwp_lat_hflux":(0,0,193,1,None),
                "nwp_sens_hflux":(0,0,193,1,None),
                "scav_coef":(0,7,197,None,None),
                "pot_temp_2m":(0,0,2,103,2.0), # Local
                "spec_humid":(0,1,1,None,None),
                "rel_humid":(0,1,0,None,None),
                "spec_humid_2m":(0,1,0,103,2.0),
                "rel_humid_2m":(0,1,0,103,2.0),
                'pressure':(0,3,0,None,None),
                "gr_surf_pr":(0,3,0,1,None),
                "albedo":(0,19,1,1,None),
                'acc_ls_rain':(0,1,9,1,None),
                'acc_conv_rain':(0,1,10,1,None),
                'ls_rain_inten':(0,1,54,1,None),
                'cnv_rain_inten':(0,1,37,1,None),
                'tot_prec':(0,1,8,1,None),
                'prec_rate':(0,1,52,1,None),
                'soil_moisture':(2,0,9,1,None),
                'ice_fract':(1,2,7,1,None),
                'land_fract':(1,2,8,1,None),
                #'land_fract':(2,0,0,1,None), # Alternative
                         # Obly one roughness in grib
                "srf_rough_disp":(2,0,1,1,None),  # #FIXME Roughness  will be written with no distinction
                "srf_rough_met":(2,0,1,1,None),  # #FIXME Roughness 
                "srf_rough":(2,0,1,1,None),
                'BLH':(0,3,18,1,None), #Defined at surface
                'nwp_BLH':(0,3,192,1,None), #Defined at surface
                "BL_top_pr":(0,3,0,192,None),
                "Kz_1m":(0,7,193,1,None),
                "MO_len_inv":(0,7,193,1,None),
                "SILAM_sens_hfl":(0,0,11,1,None),
                "SILAM_lat_hfl":(0,0,10,1,None),
                "SW_int_dn_srf":(0,4,192,1,None),
                "SW_int_up_srf":(0,4,193,1,None),
                "fric_vel":(0,2,30,1,None),
                "cnv_vel_scale":(0,7,194,1,None),
                "turb_temp":(0,7,195,1,None),
                "humid_scale":(0,7,196,1,None),
                "Ra_res_drydep":(0,7,197,1,None),


                "cnc":(0,20,None,None,None),
                "cnc2m":(0,20,None,103,2.0),
                "dd":(0,20,None,1,None),
                "wd":(0,20,None,1,None),
                "cnc":(0,20,None,None,None),

                "od":(0,20,105,None,None),#1/m Extinction
#                "od":(0,20,106,None,None),#1/m Absoprption
                "ocd":(0,20,102,None,None),

#                #0,20 stuff -- constituents
#                "cnc_mass":(0,20,0,None,None),#kg/m3
#                "cnc_molar":(0,20,51,None,None),#Mol/m3
#                "cnc_number":(0,20,59,None,None),#1/m3
#                #"cnc_bq":(0,20,210,None,None),#Local, Mol/m3
#                "cnc_2m_mass":(0,20,0,103,2.0),
#                "cnc_2m_molar":(0,20,51,103,2.0),
#                "cnc_2m_number":(0,20,59,103,2.0),
#                #"cnc_bq":(0,20,210,103,2,0),#Local, Mol/m3
#                "cnc_mass_column":(0,20,1,10,None),#Column integrated
#                "cnc_molar_column":(0,20,1,10,None),#Column integrated
#                "cnc_number_column":(0,20,1,10,None),#Column integrated
#                "vmr":(0,20,52,None,None),#Mol/Mol
#
#                "dd_mass":(0,20,6,1,None),
#                "dd_molar":(0,20,193,1,None),
#                "dd_number":(0,20,192,1,None),
#
#                "wd_mass":(0,20,7,1,None),
#                "wd_molar":(0,20,198,1,None),
#                "wd_number":(0,20,197,1,None),

                #Pollen related quantities
                "Poll_Rdy2fly":(0,20,220,1,None),
                "alrg_Rdy2fly":(0,20,221,1,None),
                "heatsum":(0,20,222,1,None),
                "poll_left":(0,20,223,1,None),
                "poll_tot_m2":(0,20,224,1,None),
                "pollen_corr":(0,20,225,1,None),
        }
Grib2number20={ # 4.2.0.20 Table -- correction to above 
        "cnc_kg":0,
        "cnc_ug":0,
        "cnc_#":59,
        "cnc_mole":51, #Concentration in air (mol m-3)
        "cnc_bq":210,
        "cnc2m_kg":0,
        "cnc2m_ug":0,
        "cnc2m_#":59,
        "cnc2m_mole":51, #Concentration in air (mol m-3)
        "cnc2m_bq":210,
        "wd_kg":7,
        "wd_ug":7,
        "wd_#":197,
        "wd_mole":198, #Concentration in air (mol m-3)
        "wd_bq":199,
        "dd_kg":6,
        "dd_ug":6,
        "dd_#":192,
        "dd_mole":193, #Concentration in air (mol m-3)
        "dd_bq":194,
            }

#
# Returns a hash of namelists for all species form silam_chemicals.dat
def getChemData(chem_db):
        nlfile = open(chem_db, 'r')
        chemicals={}
        for  chemical in namelist.Namelist.lists_from_lines(nlfile):
                  names = chemical.get("chemical_name")
                  chemicals[names[0].upper()]=chemical  #FIXME 
        nlfile.close()
        return chemicals

#
# Guesses unit and factor and set some other attributes 
# returns dictionary  of tuples (factor, gribfactor, params)
def getVarParams(SF, varlist, chemicals, nucleides=None):
        dimlessunit=""
        myunits={"Ra_res_drydep": "s/m",  
                         "grad_Ri_nbr":dimlessunit,
                         "bulk_Ri_nbr":dimlessunit,
                         "brunt_vais_f":"1/s",
                         "Kz_1m":"m2/s",
                         "land_fract":"fraction",
                         "total_cloud":"fraction",
                         "relpollen":dimlessunit,
                         }
        varparams={}
        for var in varlist: #SF.gradsdescriptor.vars:

                longvar=SF.gradsdescriptor.longVars[var]
                varnamelist=SF.get_var_namelist(var)
#                myspecs="PM2.5 PM10".split()
                attr={"long_name":longvar}

                # Set attributes from superctl
                for k,v in  varnamelist.hash.iteritems():
                        if k=="units": continue #FIXME Ignore units form superctl
                        attr[k]=v[-1]


                
                if args.pasodoble and var.lower() in standard_names:
                        attr["standard_name"]=standard_names[var.lower()]
                factor=1.
                gribfactor=1.

                if "substance_name" in attr:
                        subst=attr["substance_name"]
                else:
                        subst=None
#                print "var, subst", var, subst
#                print sorted(chemicals.keys())
                # Try grads longname []
                m = re.search(r'\[(?P<unit>.+)\]\s*$', longvar)
                try:
                        subst_upper = subst.upper()
                except:
                        subst_upper = None
                if m:
                        attr["units"]=m.group('unit')
                elif var in myunits: # Try defined units
                        attr["units"]=myunits[var]
                
                elif "mass_unit" in attr:
                        unit=attr["mass_unit"]
                        factor=1.
                        if unit.upper() == "MOLE" and subst.upper() in chemicals:
                             #Try to recalculate to mass...
                             c=chemicals[subst.upper()]
                             if "molar_mass" in c.hash:
                                mm=c.hash["molar_mass"][-1] # last is valid
                                [mmv, mmunit] = mm.split()
                                mmval=float(mmv)
                                if mmunit.lower() != "g/mole":
                                        print "chemical %s has wierd molar mass units %s!"%(var,mmunit)
                                        raise Exception("Can handle only g/mole units for molar mass")

                                if mm.startswith('99999'):
                                        if not args.quiet:
                                           print "Warning! Chemical %s has wierd molar mass %s! Leaving it in moles..."%(var,mm)
                                        unit='mole'
                                else:
                                   if (mmval > 0.): # small positive molar mass
                                        attr["molar_mass"]=mmval
                                        attr["molar_mass_unit"]= "g/mole"
                                   else: #mole unit with negative molar mass
                                        print "chemical %s has wierd molar mass %s g/mole!"%(var,mmval)
                                        raise Exception("Do not know what to do!")
                                   unit='kg'   
                                   factor=mmval*1e-3

                        gribfactor=factor 
                        if unit=='kg' and args.micrograms:
                             unit='ug'
                             factor *= 1e9



                elif subst_upper in chemicals:
#                elif subst.upper() in chemicals:
                        c=chemicals[subst.upper()]

                        if ("pollen_type" in c.hash) and (len(c.hash["molar_mass"][-1])>0):
                                        unit='#'  # v5_1 pollen
                        elif  subst.upper().startswith("POLLEN") :
                                #v5_2 pollen 
                                unit='#'
                        elif  subst.upper().startswith("PROB") :
                                unit=''
                        elif "molar_mass" in c.hash:
                #                print c.hash
                                mm=c.hash["molar_mass"][-1] # last is valid
                                [mmv, mmunit] = mm.split()
                                mmval=float(mmv)
                                if mmunit.lower() != "g/mole":
                                        print "chemical %s has wierd molar mass units %s!"%(var,mmunit)
                                        raise Exception("Can handle only g/mole units for molar mass")

                                if mm.startswith('99999'):
                                        if not args.quiet:
                                           print "Warning! Chemical %s has wierd molar mass %s! Leaving it in moles..."%(var,mm)
                                        unit='mole'
                                else:
                                     if (mmval > 0.) and (subst!="sslt"): # small positive molar mass
                                        attr["molar_mass"]=mmval
                                        attr["molar_mass_unit"]= "g/mole"
                                     else: #already in kg
                                     #   if  (subst=="sslt"):
                                     #           print "SSLT!!!"
                                        mmval=1e3 

                                     gribfactor = mmval*1e-3
                                     if args.micrograms:
                                        unit='ug'
                                        factor=mmval*1e6
                                     else:
                                        unit='kg'   
                                        factor=mmval*1e-3
                        else: 
                                unit= 'unknown'
                                raise Exception ("Unknown molar mass for %s!"%(var))
                elif nucleides:
                    #if  subst in nucleides:
                         raise Exception ("Not yet implemented for nucleides.... %s!"%(var))
                
#                elif re.search(r'(\S+)\s*$', longvar).group(1) in myspecs:
#                     gribfactor=1.
#                     if args.micrograms:
#                        unit='ug'
#                        factor=1e9
#                     else:
#                        unit='kg'   
#                        factor=1.
                else:

                     print "Don't know how to handle variable %s! Excluding from output..."%(var)
                     varparams[var]=(-1,-1,None)
                     continue


                adjParNumber=False
                #print v, subst
                if subst:
                        v=var.lower()
                        if v.startswith("cnc_") or  v.startswith("cnc2m_"):
                                attr["units"]="%s/m3"%(unit,)
                                adjParNumber = True
                        elif v.startswith("dd_"):
                                attr["units"]="%s/m2s"%(unit,)
                                adjParNumber = True
                        elif v.startswith("mass_in_air_"):
                                attr["units"]="%s"%(unit,)
                                adjParNumber = True
                        elif v.startswith("adv_moment_"):
                                attr["units"]="%s"%(unit,)
                                adjParNumber = True
                        elif v.startswith("wd_"):
                                attr["units"]="%s/m2s"%(unit,)
                                adjParNumber = True
                        elif v.startswith("ems_"):
                                attr["units"]="%s/m3s"%(unit,)
                                adjParNumber = True
                        elif v.startswith("od_"):
                                attr["units"]="1/m"
                                factor=1.
                        elif v.startswith("ocd_"):
                                attr["units"]=dimlessunit
                                factor=1.
                        elif (v.startswith("deg_days_") or v.startswith("heatsum_") or v.startswith("start_hsum_th_")):
                                attr["units"]="degree-day"
                                factor=1.
                        elif (v.startswith("poll_amt_m2_") or v.startswith("poll_tot_m2_") or  v.lower().startswith("poll_rdy2fl")):
                                attr["units"]="%s/m2"%(unit,)
                                factor=1.
                        elif v.startswith("clim_corr_to_") or v.startswith("pollen_corr_") or v.startswith("poll_left_"):
                                attr["units"]=dimlessunit
                                factor=1.
                        elif (v.startswith("day_")):
                                attr["units"]="day"
                                factor=1.
                        elif (v.startswith("landcover")):
                                attr["units"]="eff.fr."
                                factor=0.01
                        else:
                                raise Exception ("Don't know how to handle var %s wit substance %s!"%(var, subst))

                        if subst.upper() in chemicals:
                            if "GRIB2_constituent_type" in chemicals[subst.upper()].hash:
                                attr["GRIB2_constituent_type"] = chemicals[subst.upper()].hash["GRIB2_constituent_type"][-1]
                            else:

                                attr["GRIB2_constituent_type"] = -1


                shortname=attr["quantity_short_name"]
                try:
                    (discipline, category, par_number, vertical_type, vertical_value) = Grib2params[shortname]
                    if adjParNumber: #mass_unit related ting -- several parameters fit
                          par_number=Grib2number20["%s_%s"%(shortname,unit)]
                          #print var, "%s_%s"%(shortname,unit), par_number
                    attr["GRIB2_discipline"]=discipline
                    attr["GRIB2_category"]=category
                    if par_number!=None:
                       attr["GRIB2_par_number"]=par_number
                    if vertical_type != None:
                        attr["GRIB2_surface_type"] = vertical_type
                    if vertical_value != None:
                        attr["GRIB2_vertical_value"] = vertical_value
                except KeyError: #Can't make GRIB2
                    attr["GRIB2_discipline"]=-1
                
                varparams[var]=(factor,gribfactor,attr)
#                print "Make varpaprams:", var, attr
                #END loop over variables
#        exit();
        return varparams



def ParseCoordinates(SF, args): # gets indices from args
                                # and returns the dictionary of patameters of coordinates and their
                                # subset
                                #       Rotated_pole=None| rotated_pole_dictionary
                                #       truelat, truelon, x, y, z, t
                                #       timeidx=[tmin,tmax,tstep]
                                #       ibbox=[ixmin,ixmax,iymin,iymax]
                                #       
                                # 2D arrays  truelat and truelon
        params={}  # dictionary to return
        params["global"]={} # Global file attributes
        
        # Pole
        params['rotated_pole']={}
        myproj= SF.grid.proj

        if (myproj.south_pole[1] > -89.999): # Rotated pole
                if args.verbose:
                        print "Rotated pole."
                # params["global"]["grid_title"] =  myproj.grid_title
                rp=params['rotated_pole']
                rp["grid_mapping_name"] = "rotated_latitude_longitude" ;
                rp["grid_north_pole_latitude"] =  -myproj.south_pole[1] ;
                rp["grid_north_pole_longitude"] =  ( myproj.south_pole[0] + 180 +180)%360 -180 ;
                rp["_CoordinateTransformType"] = "Projection";
                rp["_CoordinateAxisTypes"] = "GeoX GeoY";
                params["lat_s_pole"]= SF.grid.proj.south_pole[0]
                params["lon_s_pole"] =  SF.grid.proj.south_pole[1]

        params["dx"] = SF.grid.dx 
        params["dy"] = SF.grid.dy


        # Grid indexes
        gfd=SF.gradsdescriptor

        x=gfd.x
        y=gfd.y

        plon,plat=np.meshgrid(x,y)
        truelon,truelat=myproj.to_geo(plon.flatten(),plat.flatten())
        truelon = truelon.reshape(plon.shape)
        truelat = truelat.reshape(plat.shape)

        ibbox = [0,len(x),0,len(y)]
        if args.Bbox: # Subset in rotated grid
                [lomin,lomax,lamin,lamax]= map( float, args.Bbox.split(',') )
                xind=np.nonzero( (x>=lomin)&(x<=lomax) )[0] # indices within the range
                yind=np.nonzero( (y>=lamin)&(y<=lamax) )[0]

                if len(xind)<1:
                        print "Can't select rotated lon range [%g,%g]; available range is [%g,%g] "%(lomin,lomax, x.min(),x.max()) 
                        exit(-1)
                if len(yind)<1:
                        print "Can't select rotated lat range [%g,%g]; available range is [%g,%g] "%(lamin,lamax, y.min(),y.max())
                        exit(-1)
                ibbox = [xind.min(),xind.max(),yind.min(),yind.max()]

        if args.bbox: # Subset in true lat-lon
                [lomin,lomax,lamin,lamax]= map( float, args.bbox.split(',') )
                select = (truelon>=lomin) & (truelon<=lomax) & (truelat>=lamin) & (truelat<=lamax)
                xind=np.nonzero( select.any(axis=1) )[0]
                yind=np.nonzero( select.any(axis=0) )[0]
                if len(xind)<1:
                        print "Can't select bbox"
                        print "Truelat:", truelat
                        print "Truelon:", truelon
                        exit(-1)
                ibbox = [xind.min(),xind.max(),yind.min(),yind.max()]
        params["ibbox"]=ibbox
        params["x"]=x[ibbox[0]:ibbox[1]]
        params["y"]=y[ibbox[2]:ibbox[3]]
        params["truelat"]=np.array(truelat[ ibbox[2]:ibbox[3], ibbox[0]:ibbox[1] ])
        params["truelon"]=np.array(truelon[ ibbox[2]:ibbox[3], ibbox[0]:ibbox[1] ])

        # Height
        params["z"]=gfd.z
        params["nz"]=len(gfd.z)
        params["vertical"]=SF.vertical

#        params["zlist"]=np.arange(gfd.nz)

        # Time
        nt=gfd.nt
        [itmin,itmax,itstep] = [0,nt,1]
        if args.times: # Time range
                trange = map( int, args.times.split(',') )
                if len(trange) == 1:
                        [itmin]=trange
                elif len(trange) == 2:
                        [itmin,itmax]=trange
                elif len(trange) == 3:
                        [itmin,itmax,itstep]=trange
                else:
                        print "Can't parse time specification \"%s\":"%( args.trange,)
                        exit(-1)
        if itmin>=nt:
                print "Only %d time steps available, %s requested"%(nt,trange,)
                exit(-1)
        itmin=max(itmin,0)
        itmax=min(itmax,nt)
        params["it"]=[itmin,itmax,itstep]
        params["alltimes"]=gfd.times()
        params["dt"]=gfd.dt

        return params
        


        
                
                
                
def initNC(outF, gridparams, args, itrange=None): #Initializes,verifies netcdf file dimensions and
                                               # global attributes

        if not (args.absolute_time): # JSBACH is not CF-compatible
            outF.Conventions = 'CF-1.6'
        addContact(outF,args.attr)

        outF._CoordinateModelRunDate =  gridparams["alltimes"][0].strftime('%Y-%m-%dT%H:%M:%SZ')

        # Time
        if itrange:
                [itmin,itmax,itstep] = itrange
        else:
                [itmin,itmax,itstep]=gridparams["it"]
        dtimes= gridparams["alltimes"][itmin:itmax:itstep]  #  [(gfd.t0+gfd.dt*it) for it in range(itmin,itmax,itstep)]
        if args.absolute_time:
                ttype = 'd'
                tarr = [(t.year*10000.+t.month*100. + t.day + (t.hour+(t.minute +t.second/60.)/60.)/24.) for t in dtimes]
                tattr= {"units":"UTC time as %Y%m%d.%f","calendar":"proleptic_gregorian"} ; #does not work with thredds
        else: 
                ttype = 'i'
				#tarr = [int(it.strftime('%s')) for it in dtimes]
                tarr = [ round((it-epoch).total_seconds()) for it in dtimes]
                epochstr=epoch.strftime("%Y-%m-%d %H:%M:%S UTC")
                tattr = {'units':'seconds since '+epochstr,"calendar":"gregorian"};
        tattr["axis"] = "T"
        outF.createVar('time', ttype, ('time',), np.array(tarr), tattr, 1)
        
        # Height
        if not(args.surface_only): # verical present
            if args.level: 
                setattr(args,"zname","level")
            else:
                setattr(args,"zname","height")

            attr={"units":'m',"long_name":"height above surface", "_CoordinateAxisType":"Height", "axis":"Z"}
            outF.createVar(args.zname, 'f', (args.zname,), gridparams["z"], attr, 0.1)

        #       X Y axes
        setattr(args,"xname","x")
        setattr(args,"yname","y")
        if gridparams['rotated_pole']:        
                attry={"axis":"Y","standard_name":"grid_latitude","long_name":"latitude in rotated pole grid",
                                "units":"degree_north", "_CoordinateAxisType": "GeoY"}
                attrx={"axis":"X", "standard_name":"grid_longitude", "long_name":"longitude in rotated pole grid",
                                "units":"degree_east", "_CoordinateAxisType":"GeoX"}
                if args.latlondims:
                        setattr(args,"xname","rlon")
                        setattr(args,"yname","rlat")
        else:
                attry={"axis":"Y","standard_name":"latitude","long_name":"latitude",
                                "units":"degree_north","_CoordinateAxisType":"Lat"}
                attrx={"axis":"X", "standard_name":"longitude", "long_name":"longitude",
                                "units":"degree_east","_CoordinateAxisType":"Lon"}
                if args.latlondims:
                        setattr(args,"xname","longitude")
                        setattr(args,"yname","latitude")
        outF.createVar(args.yname, 'f', (args.yname,), gridparams["y"], attry, 0.0001)
        outF.createVar(args.xname, 'f', (args.xname,), gridparams["x"], attrx, 0.0001)

        if gridparams["rotated_pole"]:
                outF.createVar('rotated_pole', 'c', (), np.array(chr(32)), gridparams["rotated_pole"], 1)
#                outF.variables["rotated_pole"].data = '1'


        if args.latlon:
                attr={"long_name":"geographic_latitude", "units":'degree_north'}
                outF.createVar('lat', 'f', (args.yname,args.xname,), gridparams["truelat"], attr, 0.0001)

                attr={"long_name":"geographic_longitude", "units":'degree_east'}
                outF.createVar('lon', 'f', (args.yname,args.xname,),  gridparams["truelon"], attr, 0.0001)
        
        return outF



#
# Generates idsect for FMI with given time
def GRIB2idSect(t0,opstatus):
        return [
            86, # 0,  86:Helsinki
            0, # 1, Subcentre FMI?
            10,  # 2,  GRIB Master Tables Version Number (Code Table 1.0) 
                 # http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table1-0.shtml
            1,  # 3, GRIB Local Tables Version Number (Code Table 1.1)  0 -- master, 255-missing
            1,     # idsect[4]=Significance of Reference Time (Code Table 1.2) 
                 #     0 Analysis, 1 Start of Forecast, 2 Verifying Time of Forecast
                 #     3 Observation Time, 4-191 Reserved, 192-254 Reserved for Local Use
                 #      255 Missing
            t0.year,    #idsect[5]=Reference Time - Year (4 digits)
            t0.month,    #idsect[6]=Reference Time - Month
            t0.day,    #idsect[7]=Reference Time - Day
            t0.hour,    #idsect[8]=Reference Time - Hour
            t0.minute,    #idsect[9]=Reference Time - Minute
            t0.second,    #idsect[10]=Reference Time - Second
                # idsect[11]=Production status of data (Code Table 1.3)
            opstatus,  # 0   Operational Products
                # 1   Operational Test Products
                # 2   Research Products
                # 3   Re-Analysis Products
                # 4   THORPEX Interactive Grand Global Ensemble (TIGGE)
                # 5   THORPEX Interactive Grand Global Ensemble (TIGGE) test
                # 6-191 Reserved
                # 192-254 Reserved for Local Use
                # 255 Missing

             2,   # idsect[12]=Type of processed data (Code Table 1.4) 
                  # 1 Forecast Products 
            ]


# 
# Parses gridparams
# Returns (gdsinfo, gdtmpl) tuple to be used by ncepgrib2.addgrid
def GRIB2idGds(gp): 
    nx=len(gp["x"]); ny=len(gp["y"])
    lomin=(gp["x"][0] + 360.)%360 #Positive 
    lamin=gp["y"][0] 
    dx=gp["dx"]; dy=gp["dy"]
    lomax=(lomin+(nx-1)*dx + 360.)%360 #Positive
    lamax=lamin+(ny-1)*dy

    if gp["rotated_pole"]:
        rp=gp['rotated_pole']
        splat = gp["lat_s_pole"] 
        splon = gp["lon_s_pole"] 
        gdtnum = 1 
    else:
         splat=0.
         splon=0.
         gdtnum = 0

       
        # Grid definition
    gdsinfo = [ 0, # gdsinfo[0] = Source of grid definition (see Code Table 3.0)
                     # 0 -- defined below
                nx*ny, #  gdsinfo[1] = Number of grid points in the defined grid.
                0, #  gdsinfo[2] = Number of octets needed for each additional grid points defn. 
                   #  Used to define number of points in each row for non-reg grids (=0 for regular grid).
                0, # gdsinfo[3] = Interp. of list for optional points defn (Code Table 3.11) 
                   # 0 No append list
                gdtnum, # gdsinfo[4] = Grid Definition Template Number (Code Table 3.1)
                   # 0 Latitude/Longitude (See Template 3.0)    Also called Equidistant Cylindrical or Plate Caree
                   # 1 Rotated Latitude/Longitude (See Template 3.1)
                   ]
    #GRID DEFINITION TEMPLATE 3.1 Rotate Latitude/Longitude (or equidistant cylindrical, or Plate Carree)
    gdtmpl = [
        6, # 15  Shape of the Earth (See Code Table 3.2)
            0xff, # 16  Scale Factor of radius of spherical Earth
            0xffffffff, # 17-20   Scale value of radius of spherical Earth
            0xff, # 21  Scale factor of major axis of oblate spheroid Earth
            0xffffffff, # 22-25   Scaled value of major axis of oblate spheroid Earth
            0xff, # 26  Scale factor of minor axis of oblate spheroid Earth
            0xffffffff, # 27-30   Scaled value of minor axis of oblate spheroid Earth
        nx, # 31-34   Ni-number of points along a parallel     
        ny, # 35-38   Nj-number of points along a meridian
            0, # 39-42   Basic angle of the initial production domain (see Note 1)
            0xffffffff, # 43-46   Subdivisions of basic angle used to define extreme longitudes and latitudes, and direction increments (see Note 1)
        int(lamin*1e6), # 47-50   La1-latitude of first grid point (see Note 1)
        int(lomin*1e6), # 51-54   Lo1-longitude of first grid point (see Note 1)
        int(0x30), # 55  resolution and component flags (see flag table 3.3) 
                    # 0x0c - uv  relative to easterly and northerly 
                    # 0x1c - uv  relative to grid
        int(lamax*1e6), # 56-59   La2-latitude of last grid point (see Note 1)
        int(lomax*1e6), # 60-63   Lo2-longitude of last grid point (see Note 1)
        int(abs(dx*1e6)) ,    # 64-67   Di-i direction increment (see Notes 1 and 5)
        int(abs(dy*1e6)),    # 68-71   Dj-j direction increment (see Note 1 and 5)
        0x40, # 72  Scanning mode (flags - see Flag Table 3.4)
            # These are only present if grid is rotated gdsinfo[4] == 1
        int(splat*1e6),   #73-76 Latitude of the southern pole of projection
        int(splon*1e6),    #77-80   Longitude of the southern pole of projection
#        0.    #81-84 Angle of rotation of projection
        ]
#      gdtmpl = [
#            6,#15  Shape of the Earth (See Code Table 3.2)
#              # 0 Earth assumed spherical with radius = 6,367,470.0 m 
#            0,  #16  Scale Factor of radius of spherical Earth
#            0,  #17-20   Scale value of radius of spherical Earth
#            0, #21  Scale factor of major axis of oblate spheroid Earth
#            0, #22-25   Scaled value of major axis of oblate spheroid Earth
#            0, #26  Scale factor of minor axis of oblate spheroid Earth
#            0, #27-30   Scaled value of minor axis of oblate spheroid arth
#            nx,#31-34   Ni-number of points along a parallel
#            ny, #35-38   Nj-number of points along a meridian
#            0, #39-42   Basic angle of the initial production domain (see Note 1)
#            0, #43-46   Subdivisions of basic angle used to define extreme longitudes and latitudes, and direction increments (see Note 1)
#            lamin*1e6,#47-50   La1-latitude of first grid point (see Note 1)
#            lomin*1e6,#51-54   Lo1-longitude of first grid point (see Note 1)
#            0x1c, # 55  resolution and component flags (see flag table 3.3) 
#                        # 0x0c - uv  relative to easterly and northerly 
#                        # 0x1c - uv  relative to grid
#            lamax*1e6,#56-59   La2-latitude of last grid point (see Note 1)
#            lomax*1e6,#60-63   Lo2-longitude of last grid point (see Note 1)
#            2e6,    # 64-67   Di-i direction increment (see Notes 1 and 5)
#            1e6,    # 68-71   Dj-j direction increment (see Note 1 and 5)
#            0x02, # 72  Scanning mode (flags - see Flag Table 3.4)
#  
#                # Below are only present if grid is rotated gdsinfo[4] == 1
#            splat*1e6,   #73-76 Latitude of the southern pole of projection
#            splon*1e6,    #77-80   Longitude of the southern pole of projection
#            0.    #81-84 Angle of rotation of projection
#            ]
    return (gdsinfo, gdtmpl)

#
# Returns PRODUCT DEFINITION TEMPLATE
def GRIB2pdtmpl(varparams, levelno, leveltype, levelvars, timeunit, timestep):
        # No layers yet...
        if "GRIB2_surface_type" in varparams: #Get level info from varparams
             verttype = varparams["GRIB2_surface_type"]
             if "GRIB2_vertical_value" in varparams:
                vertval  = int(varparams["GRIB2_vertical_value"])
             else:
                 vertval = 0

        else:
             verttype = leveltype
             vertval = levelvars[levelno]
         
        if "substance_name" in varparams:
             if ("mode_distribution_type" in varparams) and (varparams["mode_distribution_type"].upper() != "GAS_PHASE"):
                pdtnum = 48
                if varparams["GRIB2_constituent_type"] < 0:
                         #print "No GRIB2_constituent_type for substance %s. Will skip in GRIB output..."%(varparams[substance_name])
                         return (-1, None) #(pdtnum, pdtmpl)
                #Analysis or Forecast at a horizontal level or in a horizontal layer at a point in time for Optical Properties of Aerosol 
                if varparams["mode_distribution_type"].upper() != "FIXED_DIAMETER":
                    print varparams["mode_distribution_type"]    
                    print "Only FIXED_DIAMETER modes supported by now..."
                diam=float(( varparams["fix_diam_mode_mean_diameter"].split())[0])
                if "optical_wavelength" in varparams: # True optical property
                     wavelength = int(float(( varparams["optical_wavelength"].split())[0])*1000)
                     inttype = 11 # Equal to first limit
                     sf = -9       # Scalefactor -- nanometers
                else: # Just a replacement of buggy PDT 4.44
                      wavelength = 0
                      inttype = 255 # missing
                      sf = 0
                pdtmpl = [ # for pdtnum 48
                    varparams["GRIB2_category"],      # 10 Parameter category (see Code table 4.1)
                    varparams["GRIB2_par_number"],     # 11 Parameter number (see Code table 4.2)
                    varparams["GRIB2_constituent_type"], #12-13
                    11, #14 TYPE OF INTERVAL : 11 -- first size
                    9,     # 15 Scale factor of first size
                    diam * 1000,#16-19 Scale value of first size in meters
                    0,     #20 Scale factor of second size
                    0,     #21-24 Scale value of second size in meters
                    inttype, #25      Type of interval for first and second wavelength (see Code table 4.91) 
                    sf, #26      Scale factor of first wavelength
                    wavelength, #27-30   Scale value of first wavelength in meters
                    0, #31 Scale factor of second wavelength
                    0, #32-35   Scale value of second wavelength in meters
                    2, #36      Type of generating process (see Code table 4.3)
                    0, #37      Background generating process identifier (defined by originating centre)
                    0, #38      Analysis or forecast generating process identifier (see Code ON388 Table A)
                    0, #39-40 Hours of observational data cutoff after reference time (see Note)
                    0, #41      Minutes of observational data cutoff after reference time (see Note)
                    timeunit, #42      Indicator of unit of time range (see Code table 4.4)
                    int(timestep), #43-46   Forecast time in units defined by octet 42
                    verttype, #47 Type of first fixed surface (see Code table 4.5)
                    0, #48      Scale factor of first fixed surface
                    vertval, #49-52   Scaled value of first fixed surface
                    255,          #53 Type of second fixed surfaced (see Code table 4.5)
                    255,          #54 Scale factor of second fixed surface
                    255,          #55-58 Scaled value of second fixed surfaces
                ]

# 4.48 with missing wavelength to be used instead                
#             elif "mode_distribution_type" in varparams:
#                #Particles
#                pdtnum = 44 # Analysis or forecast at a horizontal level or in a horizontal layer at a point in time.
#                if varparams["mode_distribution_type"].upper() != "FIXED_DIAMETER":
#                    print "Only FIXED_DIAMETER modes supported by now..."
#                diam=float(( varparams["fix_diam_mode_mean_diameter"].split())[0])
#                pdtmpl = [ # for pdtnum 0
#                      varparams["GRIB2_category"],      # 10 Parameter category (see Code table 4.1)
#                      varparams["GRIB2_par_number"],     # 11 Parameter number (see Code table 4.2)
#                      varparams["GRIB2_constituent_type"],# 12-13 Aerosol Type (see Code table 4.233)
#                      11,                         #14 Type of interval for first and second size (see Code table 4.91)
#                      -9,      #15 Scale factor of first size  -9 -> nanometers
#                      diam * 1000,#16-19 Scale value of first size in meters (um to nm)
#                      0xff,     #20 Scale factor of second size
#                      0xffffffff,     #21-24 Scale value of second size in meters
#                      2,            #25 Type of generating process (see Code table 4.3)
#                      0,            #26 Background generating process identifier (defined by originating centre)
#                      0,            #27 Analysis or forecast generating process identifier (see Code ON388 Table A)
#                      0,            #28-29 Hours of observational data cutoff after reference time (see Note)
#                      0,            #30 Minutes of observational data cutoff after reference time (see Note)
#                      timeunit,      #31 Indicator of unit of time range (see Code table 4.4)
#                      int(timestep),      #32-33 Forecast time in units defined by octet 31
#                      verttype,      #34 Type of first fixed surface (see Code table 4.5)
#                      0,            #35 Scale factor of first fixed surface
#                      vertval,      #36-39 Scaled value of first fixed surface
#                      0xff,          #40 Type of second fixed surfaced (see Code table 4.5)
#                      0xff,          #41 Scale factor of second fixed surface
#                      0,          #42-45 Scaled value of second fixed surfaces
#                    ]

             else: # Gas, PM10 or PM2.5
                pdtnum = 40 # Analysis or forecast at a horizontal level or in a horizontal layer at a point in time.
                if varparams["GRIB2_constituent_type"] < 0:
                         #print "No GRIB2_constituent_type for substance %s. Will skip from output..."%(varparams["substance_name"])
                         return (-1, None) #(pdtnum, pdtmpl)
                pdtmpl = [ # for pdtnum 0
                      varparams["GRIB2_category"],      # 10 Parameter category (see Code table 4.1)
                      varparams["GRIB2_par_number"],     # 11 Parameter number (see Code table 4.2)
                      varparams["GRIB2_constituent_type"],# 12-13 Atmospheric Chemical Constituent Type (see Code table 4.230)
                      2,      # 14 Type of generating process (see Code table 4.3) 2-- forecasst
                      0,      # 15 Background generating process identifier (defined by originating centre)
                      0,      # 16 Analysis or forecast generating process identified (see Code ON388 Table A)
                      0,      # 17-18 Hours of observational data cutoff after reference time (see Note)
                      0,      # 19 Minutes of observational data cutoff after reference time (see Note)
                      timeunit, #hour      # 20 Indicator of unit of time range (see Code table 4.4)
                      int(timestep), #N of hours       # 21-24 Forecast time in units defined by octet 18
                      verttype,  # 25 Type of first fixed surface (see Code table 4.5)
                      0,      # 26 Scale factor of first fixed surface
                      vertval,# 27-30 Scaled value of first fixed surface
                      255,    # 31 Type of second fixed surfaced (see Code table 4.5)
                      255,    # 32 Scale factor of second fixed surface
                      255,    # 33-36 Scaled value of second fixed surface
                   ]
                 

             
        else:  #No substance
            pdtnum = 0 # Analysis or forecast at a horizontal level or in a horizontal layer at a point in time.
#            print "varparams:", varparams
            pdtmpl = [ # for pdtnum 0
                      varparams["GRIB2_category"],      # 10 Parameter category (see Code table 4.1)
                      varparams["GRIB2_par_number"],     # 11 Parameter number (see Code table 4.2)
                      2,      # 12 Type of generating process (see Code table 4.3) 2-- forecasst
                      0,      # 13 Background generating process identifier (defined by originating centre)
                      0,      # 14 Analysis or forecast generating process identified (see Code ON388 Table A)
                      0,      # 15-16 Hours of observational data cutoff after reference time (see Note)
                      0,      # 17 Minutes of observational data cutoff after reference time (see Note)
                      timeunit, #hour      # 18 Indicator of unit of time range (see Code table 4.4)
                      int(timestep), #N of hours       # 19-22 Forecast time in units defined by octet 18
                      verttype,  # 23 Type of first fixed surface (see Code table 4.5)
                      0,      # 24 Scale factor of first fixed surface
                      vertval,# 25-28 Scaled value of first fixed surface
                      255,    # 29 Type of second fixed surfaced (see Code table 4.5)
                      255,    # 30 Scale factor of second fixed surface
                      255,    # 31-34 Scaled value of second fixed surface
            ]

        return (pdtnum, pdtmpl) 
          







#
# Creates and writes data variables 
# Prepare VarParams should be called before to set fargs
def writeDataWithParams(outF, outB, gridparams, args, fargs, itrange=None):

        if itrange:
                [itmin,itmax,itstep] = itrange
        else:
                [itmin,itmax,itstep]=gridparams["it"]
        
        ntout=len(range(itmin,itmax,itstep))
        [ixmin,ixmax,iymin,iymax] = gridparams["ibbox"]
        FillValue=np.nan

        nz=gridparams["nz"]
        #Grib related stuff
        if outB!= None:
          t0=gridparams["alltimes"][0]
          idsect=GRIB2idSect(t0,args.prod_status)
          (gdsinfo, gdtmpl)=GRIB2idGds(gridparams)
          layertemplates=[]
#           if (params["vertical"].vertical_method == "custom_levels"):

            #
            # Layers should be treated properly
            # FIXME Now only levels processed
            #

          vert=gridparams["vertical"]
          if vert.level_type == "height_from_surface":
                   leveltype=103
                   levelvals=vert.midpoints()
                   coordlist=None
          elif vert.level_type == "hybrid":
                   leveltype=119
                   levelvals=arange(len(gridparams["z"]))
                   #FIXME  Midpoints... Should be layers
                   coordlist= [(vert.a_half[i] + vert.a_half[i+1])*0.5 for i in range(nz)]
                   coordlist= coordlist+ [(vert.b_half[i] + vert.b_half[i+1])*0.5 for i in range(nz)]
          else: 
             print "Unknown level_type", vert.level_type
             exit(-1)
#          print  "Vertical", leveltype, levelvals, coordlist

          tau = gridparams["dt"]
          print  "gridparams['dt']", gridparams["dt"]
          if tau == "__mo__":
                gribtimestep=3 #Month
                tfactor=1 #gribunits per grads step
          else:
              tau=int(tau.seconds)
              gribtimestep=13; tfactor=tau; #Seconds
              if (tau % 60 == 0): 
                  gribtimestep=0; tfactor=tau/60; #Minutes
              if (tau % 3600 == 0):
                  gribtimestep=1; tfactor=tau/3600; #Hour
              if (tau % (24*3600) == 0):
                  gribtimestep=2; tfactor=tau/24/3600; #Day
#          print  gribtimestep, tfactor, tau
#          exit()


        
        if "PM" in args.outparams: # Create variables for PM and zero them
#            if args.verbose:
#                    print "Creating PM vars"
            pmvars={}  # Dispatcher for PM accumulators
            attrlist='long_name units components substance_name quantity_short_name'.split()
            if outF!= None: 
                outF.history += "Creating PM variables:"+",".join(sorted(args.outparams["PM"].keys())) + "\n"
                for pmvar,pmpar in args.outparams["PM"].iteritems():
                    if (pmpar["nz"] == 1) or args.surface_only : #Surface variable
                            #args.yname and args.xname are set by initNC
                            pmvars[pmvar] = outF.createVariable(pmvar, 'f', ('time',args.yname,args.xname,))
                            pmvars[pmvar][:] = np.zeros((ntout,iymax-iymin,ixmax-ixmin,), dtype=float) 
                    else:
                            pmvars[pmvar] = outF.createVariable(pmvar, 'f', ('time',args.zname,args.yname,args.xname,))
                            pmvars[pmvar][:] = np.zeros((ntout,nz,iymax-iymin,ixmax-ixmin,), dtype=float)
                    
                    pmvars[pmvar]._attributes.update({k: pmpar[k] for k in attrlist})
                    if args.pasodoble and pmvar.lower() in standard_names: # Happen to be here...
                            setattr(pmvars[pmvar], "standard_name",standard_names[pmvar.lower()])
                    pmvars[pmvar]._FillValue = FillValue
            else:
                # no Netcdf variable -- use temporary arrays for  GRIB
                for pmvar,pmpar in args.outparams["PM"].iteritems():
                    if (pmpar["nz"] == 1) or args.surface_only : #Surface variable
                       pmvars[pmvar]=np.zeros((ntout,iymax-iymin,ixmax-ixmin,), dtype=float)
                    else:
                       pmvars[pmvar]=np.zeros((ntout,nz,iymax-iymin,ixmax-ixmin,), dtype=float)




        for farg in fargs:
               gfile={}
               outlist=[] # ncnames  will appear as such in output 
               readlist=[] # gradsnames need to read (Can appear as PM only)
               SF=farg.SF
               gfd=SF.gradsdescriptor
               if args.verbose:
                      print "processing", farg.super_ctl
               for var in gfd.vars: # Create variables
                    (factor, gribfactor, attr) = farg.varparams[var]

                    if farg.tag: # mark the origin of this variable
                            attr["tag"]=farg.tag

                    if gridparams['rotated_pole']:
                            attr["grid_mapping"] = "rotated_pole"

                    outp = farg.outparams[var]
#                    print "outp"
#                    print outp
#                    print outp['ncvar']
                    if outp["need"]: readlist.append(var)
                    if not outp["out"]:
                       if  args.verbose:
                           print "Skipping variable %s"%(ncvar,)
                       continue
                    ncvar = outp['ncvar'] # netcdf variable name

                    outlist.append(ncvar)
                    if outF:
                        if ncvar in outF.variables:
                                out=outF.variables[ncvar] 
                                if not args.quiet:
                                        print "Variable %s exists. Will be overwritten"%(ncvar,)
                                        # Can happen if appending to the file
                        else: # create it
                           if (gfd.nz[var] == 1) or args.surface_only : #Surface variable
                                #args.yname and args.xname are set by initNC
                                #print "creating 3D", var
                                out = outF.createVariable(ncvar, 'f', ('time',args.yname,args.xname,))
                           else:
                                #print "creating 4D", var
                                out = outF.createVariable(ncvar, 'f', ('time',args.zname,args.yname,args.xname,))
                        out._attributes.update(attr)
                        out._FillValue = FillValue

               if outF:
                  outF.history += "Adding variables from "+farg.super_ctl+":" + ",".join(sorted(outlist)) + "\n"

               for var in readlist: # Get readers
                      #  if args.verbose:
                      #      print "Getting reader", var 
                        if args.surface_only and gfd.nz[var] > 1:
                                gfile[var]=SF.get_reader(var, ind_level=0,  projected=True, mask_mode='none')
                        else:
                                gfile[var]=SF.get_reader(var, ind_level=None,  projected=True, mask_mode='none') 
                        gfile[var].seek(itmin) # goto fist needed time
                        
               itout=0
               for it in range(itmin,itmax,itstep): #Write!

                   if outF:
                      # Make sure we assign the fields to proper time 
                      nctime =  outF.variables["time"].data[itout]
                      gradstime = round( (gridparams["alltimes"][it]-epoch).total_seconds())
                      #if args.verbose:
                      if not args.quiet:
                           print "%5d %5d "%(it, itout), 
                           print gridparams["alltimes"][it].strftime('%Y-%m-%d %H'),
                           print gradstime, nctime,
                           print dt.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC")
                      if (gradstime != nctime ):
                           print "Gradstime is different from nctime :((("
                           print "Gradstime:", gradstime
                           print "NCtime:", nctime 
                           exit()
                   else:
                      # if args.verbose:
                      if not args.quiet:
                           print "%5d %5d "%(it, itout), 
                           print gridparams["alltimes"][it].strftime('%Y-%m-%d %H'),
                           print "now:",
                           print dt.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC")


                   for var in readlist:
                        outp = farg.outparams[var] # netcdf variable name
                        (factor, gribfactor, attr) = farg.varparams[var]
                        outp = farg.outparams[var]
                        data = gfile[var].read(1)

#                        print var, attr
                        ifPutToGrib =  (outB != None and outp["out"] and attr["GRIB2_discipline"]>=0)
                        ifPutToNC   =  (outF != None and outp["out"])
                        if (outB != None and outp["out"] and attr["GRIB2_discipline"]<0):
                                print "Variable %s will be skipped from GEIB2 output: No discipline.."%(var)
                        
                        if ifPutToGrib: #Prepare GRIB2 message
                                m=ncepgrib2.Grib2Encode(attr["GRIB2_discipline"],idsect)
                                m.addgrid(gdsinfo, gdtmpl, deflist=None)

                        if (gfd.nz[var] == 1) or args.surface_only: #data is 2d
                             a = data.transpose() 
                             a[np.abs(a - gfd.undef) < 1e-5] = FillValue
                             if ifPutToNC:
                                   out=outF.variables[outp['ncvar']] # netcdf variable name
                                   out[itout,:,:] = a[iymin:iymax,ixmin:ixmax]*factor
                             if ifPutToGrib:
                                 (pdnmpnum,pdtmpl) = GRIB2pdtmpl(attr, 0, leveltype, levelvals, gribtimestep, it*tfactor)
 #                                print "gribtimestep, it, tfactor:",  gribtimestep, it, tfactor
                                 if pdnmpnum >=0 :
                                         m.addfield(pdnmpnum,pdtmpl,4,[1],
                                                 a[iymin:iymax,ixmin:ixmax]*factor,
                                                 coordlist=coordlist)
                                         m.end()
                                         outB.write(m.msg)
                                 else:
                                         print "Skipping %s from GRIB2 output: No constituent.."%(var,)

                             for acc in outp["pm_acc"]:
                                   pmvars[acc][itout,:,:] += a[iymin:iymax,ixmin:ixmax]*factor
                        else:   #  3d  -> 4d  
                             a = data.transpose(2,1,0)
                             a[np.abs(a - gfd.undef) < 1e-5] = FillValue
                             if ifPutToNC:
                                out=outF.variables[outp['ncvar']] # netcdf variable name
                                out[itout,:,:,:] = a[:,iymin:iymax,ixmin:ixmax]*factor
                             if ifPutToGrib:
                                 for iLev in range(gfd.nz[var]):
                                   (pdnmpnum,pdtmpl) = GRIB2pdtmpl(attr, iLev, leveltype, levelvals, gribtimestep, it*tfactor)
#                                   print "iz, gribtimestep, it, tfactor:", iLev, gribtimestep, it, tfactor
                                   if pdnmpnum >= 0 :
                                     m.addfield(pdnmpnum,pdtmpl,4,[1],
                                        a[iLev,iymin:iymax,ixmin:ixmax]*factor, 
                                         coordlist=coordlist)
                                 if pdnmpnum >= 0 :
                                   m.end()
                                   outB.write(m.msg)
                                 else:
                                   print "Skipping %s from GRIB2 output: No constituent.."%(var,)

                             for acc in outp["pm_acc"]:
                                     pmvars[acc][itout,:,:,:] += a[:,iymin:iymax,ixmin:ixmax]*factor

                        gfile[var].seek(itstep-1) # skip if needed...
                        #END of variable loop
                   itout += 1
               
               for var in gfile.keys():
                        gfile[var].close()
        #Loop over input files over
        #Write  PM to GRIB 
        if outB != None and "PM" in args.outparams:
           itout=0
           for it in range(itmin,itmax,itstep): #Write!
              for pmvar,pmpar in args.outparams["PM"].iteritems():
                m=ncepgrib2.Grib2Encode(pmpar["GRIB2_discipline"],idsect)
                m.addgrid(gdsinfo, gdtmpl, deflist=None)
                if (pmpar["nz"] == 1) or args.surface_only : #Surface variable
                   (pdnmpnum,pdtmpl) = GRIB2pdtmpl(pmpar, 0, leveltype, levelvals, gribtimestep, it*tfactor)
                   if pmpar["units"].startswith("ug"):
                      m.addfield(pdnmpnum,pdtmpl,4,[1], pmvars[pmvar][itout,:,:]*1e-9, coordlist=coordlist)
                   else:
                      m.addfield(pdnmpnum,pdtmpl,4,[1], pmvars[pmvar][itout,:,:], coordlist=coordlist)
                   m.end()
                   outB.write(m.msg)
                else:
                   for iLev in range(gfd.nz[var]):
                      (pdnmpnum,pdtmpl) = GRIB2pdtmpl(pmpar, iLev, leveltype, levelvals, gribtimestep, it*tfactor)
                      if pmpar["units"].startswith("ug"):
                         m.addfield(pdnmpnum,pdtmpl,4,[1], pmvars[pmvar][itout,iLev,:,:]*1e-9, coordlist=coordlist)
                      else:
                         m.addfield(pdnmpnum,pdtmpl,4,[1], pmvars[pmvar][itout,iLev,:,:], coordlist=coordlist)
                   m.end()
                   outB.write(m.msg)
              itout += 1

# To be used in PrepareOutParams only
# Make sure that units are the same, generate longname etc..
def AssertPM(outpars, acc, accSP,  varname, attr, nz): #global outparams, accumulator, 
                                        #PMname, varname, variable attr, nz in original var
        myunits=attr["units"]
        if not "PM" in outpars:
                outpars["PM"]={}
        if not acc in outpars["PM"]:
                longname = re.sub(" \S+$", " "+accSP, attr["long_name"])# replace species with proper PM 
                outpars["PM"][acc]={"units":myunits, 
                                   "components":varname,
                                   "long_name":longname, 
                                   "nz":nz, 
                                   "substance_name":"PM",
                                    "quantity_short_name":attr["quantity_short_name"],
                                    "GRIB2_discipline":0,
                                    "GRIB2_category":20,
                                    "GRIB2_par_number":attr["GRIB2_par_number"],
                                    }
                try:
                  outparams["PM"][acc]["GRIB2_surface_type"] = attr["GRIB2_surface_type"]
                except:
                   pass
                try:
                  outparams["PM"][acc]["GRIB2_vertical_value"] = attr["GRIB2_vertical_value"]
                except:
                   pass
        else:
                outpars["PM"][acc]["components"] += " + " + varname

        pmunits=outpars["PM"][acc]["units"]
        if pmunits != myunits: # Should be the same units...
                print "Error! Units of %s [%s] do not match the units of %s [%s]!"%(varname,myunits,acc,pmunits)
                exit(-1)

        pmnz=outpars["PM"][acc]["nz"]
        if pmnz != nz: # Should be the same vertical dim...
                print "Error! nz of %s (%d) do not match nx of %s (%s)!"%(varname,nz,acc,pmnz)
                exit(-1)



def PrepareOutParams(args,fargs): # arguments namespace, array of fargs namespace
        alloutvars=[] # all variables for the output
        allpmvars=[] # all variables for total PM
        chemdata=getChemData(args.chem_db)
        nucleides=None
        
        diam=re.compile("^\s*([\.0-9]+)\s+(mkm|um)\s*$") # re to extract mode diameter

        setattr(args,"outparams",{}) # global out parameters


        for myargs in reversed(fargs): #last wins
                SF=sf.Silamfile(myargs.super_ctl)
                gfd=SF.gradsdescriptor
                outparams={}
                varlist = gfd.vars
                
                onlyre=[];  skipre = []; 
                if myargs.only_vars:
                    onlyre = [ re.compile("^"+v+"$") for v in myargs.only_vars.split(",") ]
                if myargs.drop_vars:
                    skipre = [ re.compile("^"+v+"$") for v in myargs.drop_vars.split(",") ]

                # To match substances
                pmspre=[];  pmvarre=[]; pmskipre=[];
                if  myargs.pm_species:
                    pmspre  = [ re.compile("^"+v+"$") for v in myargs.pm_species.split(",") ]
                if  myargs.pm_vars:
                    pmvarre = [ re.compile("^"+v+"$") for v in myargs.pm_vars.split(",") ]
                if  myargs.skip_pm_vars:
                    pmskipre = [ re.compile("^"+v+"$") for v in myargs.skip_pm_vars.split(",") ]
                
                varparams   = getVarParams(SF, varlist, chemdata, nucleides)
                
                
                for var in varlist: # Select vars
                     outparams[var]={"need":False, "out":False, "pm_acc":[]}
                     (factor,gribfactor,attr)=varparams[var]
#                     print "ATTR", attr 
                     if not attr:  # Ignore unknown variable
                         continue

                     # get full quantity name (to be used as netcdf variable)
    #                 print attr
                     shortname=attr["quantity_short_name"]
                     varname=shortname
                     if "substance_name" in attr:
                             subst=attr["substance_name"] 
                             varname += "_"+subst
                     else:
                         subst=None

                     dp=None
                     if "mode_distribution_type" in attr and attr["mode_distribution_type"].lower() != "gas_phase": # particle
                        outparams[var]['pm']=True
                        # Should compalin if something is wrong with the diameter string
                        if "moving_diam_mode_mean_diameter" in attr:
                            dp=float(diam.match(attr["moving_diam_mode_mean_diameter"]).group(1))
                        elif "fix_diam_mode_mean_diameter" in attr:
                            dp=float(diam.match(attr["fix_diam_mode_mean_diameter"]).group(1))
                        else:
                            print "Warning! Variable %s (%s) from %s: "%(var, varname, myargs.super_ctl)
                            print "Can't figure out mean diameter....."
                            print attr


                        if dp<1.:
                                mode="m_%02d"%(int(dp*100)) #m_05
                        elif dp<10.:
                                mode="m%d_%.0f"%(int(dp),(dp%1)*10) #m2_5
                        else: #should not happen
                                mode="m%.0f"%(dp)  # m125
                        varname += "_"+mode
                     # FIXME Wavelength-related naming should be treated here

                     # Select requested
                     if (( myargs.only_vars and not any((r.match(var) or r.match(varname)) for r in onlyre)   )  #Not requested
                           or  any(r.match(var) or r.match(varname) for r in skipre)   # Requested to remove
                           or  (varname in alloutvars)                             # Will be overwritten in output
                            ):
                            #print onlyre, myargs.only_vars.split(","), bool(onlyre)
                            #print  any((r.match(var) or r.match(varname)) for r in onlyre)   
                            
                            if not args.quiet and (varname in alloutvars):
                                 print "Warning! Variable %s (%s) from %s will be overwritten: "%(var, varname, myargs.super_ctl)
                     else:
                                outparams[var]['need']=True 
                                outparams[var]['out']=True  
                                alloutvars.append(varname)
                    
                     if args.total_pm and subst and dp: # More flags to make total PM
                            
                             if False:
                                     print "PM SPECIES:", var
                                     print [r.match(subst) for r in pmspre]
                                     print myargs.pm_species
                                     print [r.match(subst) for r in pmvarre]
                                     print (varname in allpmvars)
                                     print [r.match(var) or r.match(varname) for r in pmskipre]
                                     print myargs.skip_pm_vars,   [r.pattern for r in pmskipre]

                             if  (  not  ( any(r.match(subst) for r in pmspre) or any(r.match(varname) for r in pmvarre))
                                             # not requested
                                   or any(r.match(var) or r.match(varname) for r in pmskipre)         
                                   or  (varname in allpmvars)                               # Already accounted in output
                                 ):
                                 if not args.quiet and (varname in allpmvars):
                                      print "Warning! Variable %s from %s will not be accounted: "%(var, myargs.super_ctl)
                             else: # PM to account
                                if dp < 10:
                                    pmname="PM10"; accname="%s_%s"%(shortname,pmname) # Account as pm10
                                    if (( myargs.only_vars and not any(r.match(accname) for r in onlyre)   )  #Not requested
                                            or  any(r.match(accname) for r in skipre) ):  # Requested to remove
                                        if args.verbose:
                                            print "Skipping ",accname
                                    else:
                                        outparams[var]['pm_acc'].append(accname)
                                        AssertPM(args.outparams, accname, pmname,  varname, attr, gfd.nz[var]) #  Global accumulator
    #                                   args.outparams["PM"][accname]["mode_distribution_type"]="MAX_DIAMETER"
    #                                   args.outparams["PM"][accname]["fix_diam_mode_max_diameter"]="10 mkm"
                                        args.outparams["PM"][accname]["GRIB2_constituent_type"]=40008
                                        outparams[var]['need']=True #Have to read
                                        allpmvars.append(var)
                                if dp < 2.5: # Account also as pm2.5
                                    pmname="PM2_5"; accname="%s_%s"%(shortname,pmname)
                                    if (( myargs.only_vars and not any(r.match(accname) for r in onlyre)   )  #Not requested
                                            or  any(r.match(accname) for r in skipre) ):  # Requested to remove
                                        if args.verbose:
                                            print "Skipping ",accname
                                    else:
                                        outparams[var]['pm_acc'].append(accname) 
                                        AssertPM(args.outparams, accname, pmname,  varname, attr, gfd.nz[var])
    #                                   args.outparams["PM"][accname]["mode_distribution_type"]="MAX_DIAMETER"
    #                                   args.outparams["PM"][accname]["fix_diam_mode_max_diameter"]="2.5 mkm"
                                        args.outparams["PM"][accname]["GRIB2_constituent_type"]=40009
                                if (dp > 10) and not args.quiet: # bark
                                        print "Warning! %s(%s) from %s is too coarse: %f um. Not accounting as PM10."%(var, varname, myargs.super_ctl, dp)
                      
                     outparams[var]['ncvar']=varname

                #end variable loop
                setattr(myargs,"SF",SF)
                setattr(myargs,"varparams",varparams)
                setattr(myargs,"outparams",outparams)
        
        if args.total_pm and (not "PM" in args.outparams) and (not args.quiet):
                print "Warning! PM requested but no PM components found..."
                
        if args.verbose:
            print "Arranged output:"
            for myargs in fargs:
                    print "  From %s taking: "%(myargs.super_ctl,)
                    for var in sorted(myargs.outparams.keys()):
                       print "        ",var, myargs.outparams[var]
            if "PM" in args.outparams:
                    print "Generating total PM:"
                    for var in sorted(args.outparams["PM"].keys()):
                            print "        ",var, args.outparams["PM"][var]

# 
#
def MakeNC(outname_nc, outname_grib, gridparams, args, fargs, histappend, itrange=None):
    if outname_nc != None:
        if not args.quiet:
              print "Writing to:", outname_nc
        outF = silamNCout.silamNCout(outname_nc, append=args.append, 
                        verbose=args.verbose, allowrecord = (not args.no_record))  
                        # Create and (opt) copy the stuff from
                        # the previous file 
        outF.history +=  histappend 
        initNC(outF,  gridparams, args, itrange=itrange) # Crate/verify dimensions 
    else:
       outF = None

    if outname_grib != None:
       if not args.quiet:
              print "Writing to:", outname_grib
       outB = open(outname_grib,"wb")
    else:
       outB = None

    writeDataWithParams(outF, outB, gridparams,  args, fargs, itrange=itrange) # Create/write data vars
#                print "Commit"
        
    if outname_nc != None:
        outF.finalize() # finalize!
    if outname_grib != None:
        outB.close()


if __name__ == "__main__":

        args,fargs=ProcessCmdLine() # Parse arguments and pre-process them 
        
        PrepareOutParams(args,fargs)
        
        gridparams  = ParseCoordinates(fargs[0].SF, args)
        gfd=fargs[0].SF.gradsdescriptor

#        for var in varparams.keys():
#                print var, varparams[var]
#        exit(0)
        outname_nc=None
        outname_grib=None

        if args.output_nc:
           outtemplate_nc = gfd.t0.strftime(args.output_nc) # add time template
           outname_nc =  outtemplate_nc
        if args.output_grib: 
           outtemplate_grib = gfd.t0.strftime(args.output_grib) # add time template
           outname_grib = outtemplate_grib
        if args.verbose:
                print "Out template nc:", outname_nc
                print "Out template grib:", outname_grib

        currentTimestr=dt.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC")
        histappend = currentTimestr+': Commandline:'+' '.join(sys.argv)+"\n"
        if  args.steps_per_file:
                [itmin,itmax,itstep]=gridparams["it"]
                timeinds=range(itmin,itmax,itstep)
                if args.verbose:
                    print "requested time idx:",gridparams["it"]
                    print timeinds, range(0,len(timeinds),args.steps_per_file)

                for ii in range(0,len(timeinds),args.steps_per_file): 
                        myitmin=timeinds[ii]
                        myitmax=min(myitmin+args.steps_per_file, itmax) # Not beyond the file
                        itrange = list([myitmin,myitmax,itstep]) # This range overrides gridparams["it"]
                        if args.date:  # datetime template -- reset datetime
                            if outname_nc:
                                outname_nc=gridparams["alltimes"][myitmin].strftime(args.output_nc)
                            if outname_grib:    
                                outname_grib=gridparams["alltimes"][myitmin].strftime(args.output_grib)
                        else: # timestep template
                            if outname_nc:   
                               outname_nc=outtemplate_nc%(myitmin,) # timestep template
                            if outname_grib:   
                               outname_grib=outtemplate_grib%(myitmin,) # timestep template
                                                # Need to copy gridparams...
                        if   args.num_threads > 0:
                           print "ITRANGE", itrange

                           t=threading.Thread(group=None, target=MakeNC, name=None, 
                                       args=(outname_nc, outname_grib, gridparams, args, fargs, histappend), 
                                       kwargs={"itrange":itrange}) 
                           t.start()
                           while threading.active_count() > args.num_threads:
                                time.sleep(.01)
                        else:
                            MakeNC(outname_nc, outname_grib, gridparams, args, fargs, histappend, itrange=itrange) # With no threads
                if   args.num_threads > 0:
                        t.join()

        else:
                MakeNC(outname_nc, outname_grib, gridparams, args, fargs, histappend)
        print "Done!"
        exit(0)


