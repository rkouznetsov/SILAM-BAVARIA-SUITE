from __future__ import print_function
import numpy as np
import pytz
import  os, sys, argparse
from toolbox import MyTimeVars
import matplotlib.pyplot as plt
import memstations
from pylab import rcParams
import matplotlib.dates as mdt
from matplotlib import dates, use
widthP=int(24/2.54)
heightP=int(16/2.54)
rcParams['figure.figsize'] = (widthP, heightP)
rcParams.update({'figure.autolayout': True})
rcParams.update({'figure.autolayout': True})
species=['CO', 'NO', 'NO2', 'O3', 'PM10']
for j in range (0, len(species)):
    if (j<4):
     tsmatr1=MyTimeVars.TsMatrix.fromNC('/home/yury/Oper_out/ts/20200116/cnc_'+species[j]+'_gas.nc')
     tsmatr2=MyTimeVars.TsMatrix.fromNC('/home/yury/Oper_out/ts/20200117/cnc_'+species[j]+'_gas.nc')
     tsmatr3=MyTimeVars.TsMatrix.fromNC('/home/yury/Oper_out/ts/20200118/cnc_'+species[j]+'_gas.nc')
     tsmatr4=MyTimeVars.TsMatrix.fromNC('/home/yury/Oper_out/ts/20200119/cnc_'+species[j]+'_gas.nc')
     tsmatr5=MyTimeVars.TsMatrix.fromNC('/home/yury/Oper_out/ts/20200120/cnc_'+species[j]+'_gas.nc')
     tsmatr6=MyTimeVars.TsMatrix.fromNC('/home/yury/Oper_out/ts_MEM/'+species[j]+'MEM.nc')
    else:
     tsmatr1=MyTimeVars.TsMatrix.fromNC('/home/yury/Oper_out/ts/20200116/cnc_'+species[j]+'.nc')
     tsmatr2=MyTimeVars.TsMatrix.fromNC('/home/yury/Oper_out/ts/20200117/cnc_'+species[j]+'.nc')
     tsmatr3=MyTimeVars.TsMatrix.fromNC('/home/yury/Oper_out/ts/20200118/cnc_'+species[j]+'.nc')
     tsmatr4=MyTimeVars.TsMatrix.fromNC('/home/yury/Oper_out/ts/20200119/cnc_'+species[j]+'.nc')
     tsmatr5=MyTimeVars.TsMatrix.fromNC('/home/yury/Oper_out/ts/20200120/cnc_'+species[j]+'.nc')
     tsmatr6=MyTimeVars.TsMatrix.fromNC('/home/yury/Oper_out/ts_MEM/'+species[j]+'MEM.nc')

    fig=plt.figure()
    stcodes = [station.code for station in tsmatr6.stations]
    for i, code in enumerate(stcodes): 
     ax = fig.subplots()
     ax.plot(tsmatr6.times, tsmatr6.vals[:, i], color='black', linewidth=3.0, label='MEM')
     ax.plot(tsmatr2.times, tsmatr2.vals[:, i], color='red', linewidth=3.0)
     ax.plot(tsmatr3.times, tsmatr3.vals[:, i], color='gray', linewidth=3.0)
     ax.plot(tsmatr4.times, tsmatr4.vals[:, i], color='blue', linewidth=3.0)
     ax.plot(tsmatr5.times, tsmatr5.vals[:, i], color='royalblue', linewidth=3.0)
     ax.plot(tsmatr1.times, tsmatr1.vals[:, i], color='firebrick', linewidth=3.0, label='SILAM')
     ax.set_xlim(tsmatr1.times[0], tsmatr5.times[-1])
     xfmt = mdt.DateFormatter('%b %d %H Z')
     ax.xaxis.set_major_formatter(xfmt)
     ax.xaxis.set_major_locator(mdt.HourLocator(interval=24))
     ax.set_title('SILAM automatic forecasts and ' +code+' data for Moscow')
     ax.set_ylabel('Concentration'+species[j]+' mg/m3')
     ax.set_xlabel("Time UTC")
     name='Plot'+species[j]+'_ts_oper_SILAM_at_'+str(code)+'.png'
     fig.savefig(name, dpi=500)
     fig.clf()
     print(name)
