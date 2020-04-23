function silamaor(args)
# Drawing the SILAM output 
# USAGE:  silamaor ctlFName AnalString taxon OutputFile

InputCtl   = subwrd(args,1)
AnalString = subwrd(args,2)
taxon      = subwrd(args,3)
OutputFile =  subwrd(args,4)
mpdset = subwrd(args,5)

'reinit'



#some colors
'set rgb 22 254 254 254'
'set rgb 21 190 190 190'
'set rgb 20 220 220 220'
'set rgb 23 0 205 254'
'set rgb 24 0 220 0'
'set rgb 25 120 120 120'
'set rgb 26 60 60 60'
'set rgb 27 30 30 30'
#greens
'set rgb 31 50 254 50'
'set rgb 32 0 220 0'
'set rgb 33 0 190 0'
'set rgb 34 0 160 0'
'set rgb 35 0 130 0'
'set rgb 36 0 110 0'
'set rgb 37 0 90 50'
#oranges
'set rgb 41 254 200 0'
'set rgb 42 254 170 0'
'set rgb 43 254 140 0'
'set rgb 44 254 110 0'
'set rgb 45 254 80 0'
'set rgb 46 254 50 0'
'set rgb 47 254 0 0'
#light blue
'set rgb 50 150 240 200' 

say "mpdset:" mpdset
bigmapcmd=''
smallmapcmd=''
if (strlen(mpdset) >0)
   say "Setting:" mpdset
   'set mpdset 'mpdset
*        set mpt type color style thickness   
* coaslines+lakes+islands
   'set mpt 1 1 1 0.5'
* grid dashed
   'set mpt 2 1 5 0.5'
** Inland borders thick dark green
*   'set rgb  250  0  50  0'
*   'set mpt 3 15 1 2'
*   'set mpt 3 250 1 2'
* Inland borders thick dashed
   'set mpt 3 1 3 2'
* rivers -- blue 
   'set rgb  250 140 140 230'
*   'set mpt 4 4 1 0.5'
   bigmapcmd='set mpt 4 250 1 0.5'
   smallmapcmd='set mpt 4 off'
endif

#
# Variables, all possible
#
tax=substr(taxon,1,3)
vCnc = 'cnc_p'tax
vHeatsum='heatsum_p'tax
vPollAmt = 'poll_left_p'tax
vCorr = 'pollen__p'tax
vPRdy2Fly = 'poll_rdy2f_p'tax
v.2 = vPollAmt'*100'
vTit.2 = 'Pollen left in catkins, %'
ccols_v.2 = '0     35   10   7   12   8'
clevs_v.2 =   '0.01  25   50  75  99'
#
# Specifics of the taxa
#
if (taxon = "alder")
    ccols_cnc = "0 20 3 10 7 12 8 2 6 9"
    clevs_cnc = "0.1 1 5 10 25 50 100 500 1000"
    v.1 = vHeatsum
    vTit.1 = "Heat sum, degree-hour"
    ccols_v.1 = '0 21 10 31 32 33 34 35 41 42 43 44 45'
    clevs_v.1 = "0 200 500 750 1000 1300 1600 2000 2500 3000 4000 5000"
    
else
  if (taxon = "birch")
    ccols_cnc = "0 20 3 10 7 12 8 2 6 9"
    clevs_cnc = "1 5 10 25 50 100 500 1000 5000"
    v.1 = vHeatsum
    vTit.1 = "Heat sum, degree-day"
    ccols_v.1 = '0 21 10 31 32 33 34 35 36 37 41 42 43 44 45 47 6 '
    clevs_v.1 = '0 30 40 50 60 70 80 90 100 120 140 160 180 200 220 240'
    
  else
    if (taxon = "grass")
      ccols_cnc = "0 20 3 10 7 12 8 2 6 9"
      clevs_cnc = "1 5 10 25 50 100 500 1000 5000"
      v.1 = "(temp_2m-273.15)"
      vTit.1 = "Temperature 2m, C"
      ccols_v.1 = '15 4 11 5 13 3 10 7 12 8 2 6 9'
      clevs_v.1 = '8 10 12 14 16 18 20 22 24 26 28 30'
*      v.2 = "((ls_rain_inten+cnv_rain_inten)*3600)"
      v.2 = "(prec_rate*3600)"
      vTit.2 = "Precipitation rate, mm/h"
      ccols_v.2 = '0 20 11 5 3 10 7 12 8 2 6 '
      clevs_v.2 = "0.1 0.5 1 2 3 4 5 6 8 10"

    else
      if (taxon = "mugwort")
        ccols_cnc = '0 20 3 10 7 12 8 2 6 9'
        clevs_cnc = '0.1 1 5 10 25 50 100 500 1000'
        v.1 = vPRdy2Fly % "/100"
        vTit.1 = 'Pollen ready to fly, 100#/m2'
        ccols_v.1 = '0  31  32  33  34  35   41   42   43   44   45   46'
        clevs_v.1 = '0 1 2 5 10 20 50 100 200 500 1000'

      else
        if (taxon = "olive")
          ccols_cnc = "0 20 3 10 7 12 8 2 6 9"
          clevs_cnc = "0.1 1 5 10 25 50 100 500 1000"
          v.1 = vHeatsum
          vTit.1 = "Heat sum, degree day"
          ccols_v.1 = '0  31  32  33  34  35   41   42   43   44   45   46'
          clevs_v.1 = '0 200 400 600 800 1000 1200 1400 1600 1800 2000'

        else
          if (taxon = "ragweed")
            ccols_cnc = "0 20 3 10 7 12 8 2 6 9"
            clevs_cnc = "0.1 1 5 10 25 50 100 500 1000"
            v.1 = vPRdy2Fly % "/100"
            vTit.1 = "Pollen ready to fly, 100#/m2"
            ccols_v.1 = "0  31  32  33  34  35   41   42   43   44   45   46"
            clevs_v.1 = "0 1 2 5 10 20 50 100 200 500 1000"
          else
            say 'Wrong taxon ' taxon
            'quit'
          endif
        endif
      endif
    endif
  endif
endif

say InputCtl
'xdfopen  'InputCtl

'set gxout grfill'
*'set mpdset lowhir'
*'set mpdset world_map'
*'set mproj nps'
'set z 1'
'set xlab off'
'set ylab off'


'q file'
sizeLine = sublin(result,5)
timeSize = subwrd(sizeLine,12)

tStep=1
 
while (tStep <= timeSize)  
  'c'
  'set t 'tStep
  'query dims'
  rec5 = sublin(result,5)
  DateTime = subwrd(rec5,6)

#------------------------------- Concentration
#
  'set vpage 0 8.5 3.6 11'

  'set ccols 'ccols_cnc
  'set clevs 'clevs_cnc
  'set grads off'
  'set frame off'
  bigmapcmd

  'd 'vCnc
  'draw title SILAM model forecast: 'taxon' pollen\ (#/m3) 'DateTime
  'cbar'

#-------------------------- Low panels. Note switch of vpage
# 
  smallmapcmd
  'set vpage 0 4.25 0 3.6'
  iVar = 1
  while(iVar < 3)
    'set ccols 'ccols_v.iVar
    'set clevs 'clevs_v.iVar
    'set grads off'
    'set frame off'
    'd 'v.iVar
    'draw title 'vTit.iVar
    'cbar'
    'set vpage 4.25 8.5 0 3.6'
    iVar = iVar + 1
  endwhile

  outf=OutputFile'_'math_format('%03.0f',tStep-1)'.png'
  'printim 'outf' x655 y850 white'
  say outf ' Done.'

  tStep = tStep+1
*  break
endwhile

say '  'OutputFile
'quit'

