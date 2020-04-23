function prsrm(args)
*
* The function draws stuff from .nc files
*
* 'reinit'
ncFNm = subwrd(args,1)
substNm = subwrd(args,2)
substNm_out = subwrd(args,3)
factor = subwrd(args,4)
massUnit = subwrd(args,5)
AnalTimeString = subwrd(args,6)
colour = subwrd(args,7)
scale_legend = subwrd(args,8)
scale_type = subwrd(args,9)
scale_offset = subwrd(args,10)
outDir = subwrd(args,11)

* Forecast offset in hours of the first record, can be 'daily'
timeOff = subwrd(args,12)

mpdset = subwrd(args,13)

* say "mpdset:" mpdset
if (strlen(mpdset) > 0)
   'set mpdset 'mpdset
*   say result
*        set mpt type color style thickness   
* coaslines+lakes+islands
   'set mpt 1 1 1 1'
* grid dotted
   'set mpt 2 1 5 0.5'
** Inland borders thick dark green
*   'set rgb  250  0  50  0'
**   'set mpt 3 15 1 2'
*   'set mpt 3 250 1 2'
* Thin short dash
   'set mpt 3 1 3 2'
* rivers -- blue 
*   'set mpt 4 4 1 0.5'
   'set xlab off'
   'set ylab off'
else
   'set mpdset world_map'
*   say 'mpdset='mpdset',' result
endif




* Hack for retarded grads var name handling
gradsvar=substr(substNm,1,15)

*  say 'ncFNm ' ncFNm
*  say 'substNm 'substNm 
*  say 'substNm_out 'substNm_out 
*  say 'outDir ' outDir
*  say 'AnalTimeString 'AnalTimeString
*  say 'colour 'colour
*  say 'scale_legend 'scale_legend
*  say 'scale_type 'scale_type
*  say 'scale_offset 'scale_offset
*  
*  say 'timeOff ' timeOff

*
*----------------------------------- Open both ctl files and set the drawing size and common title
*
*'set cachesf 100'
* say ncFNm
'sdfopen 'ncFNm
'q file'
*say result
sizeLine = sublin(result,5)
xSize = subwrd(sizeLine,3)
ySize = subwrd(sizeLine,6)
timeSize = subwrd(sizeLine,12)
*'set x 5 'xSize-5
*'set y 5 'ySize-5
* 'set gxout shaded'
'set gxout grfill'
**'set dbuff on'
'run colors.gs 'colour


'set rgb 129 40 130 240'
'set rgb 130 102 229 102'
'set rgb 131 255 240 85'
'set rgb 132 255 187 87'
'set rgb 133 255 68 68'
'set rgb 134 182 70 139'
'set rgb 140 200 200 200'
'set rgb 141 100 100 100'



  title.1 = 'srf'



  'q dims'
  xline = sublin(result, 2)
  xmin=subwrd(xline, 6)
  xmax=subwrd(xline, 8)
  yline = sublin(result,3)
  ymin=subwrd(yline,6)
  ymax=subwrd(yline,8)

  if (xmax - xmin > 355)
   'set lon -180 180'
    IFLONGLOBAL=TRUE
  else
    IFLONGLOBAL=FALSE
  endif
*  say 'IFLONGLOBAL = 'IFLONGLOBAL
*
*--------------------------------- Hourly cycle
*
  time = 1
  if (timeOff = 'daily')
     filesuff='d'
     nbr=1
  else
     filesuff=''
     nbr = timeOff
  endif
  while(time < timeSize+1)
    'set t 'time
    'q dims'
    timeLine = sublin(result,5)
    dateTime = subwrd(timeLine,6)

* No hours for daily files -- less confusion
    if (timeOff = 'daily')
       hhmm= ' ' 
    else
       hhmm = substr(dateTime,1,2)':00'
    endif
    dd = substr(dateTime,strlen(dateTime)-8,strlen(dateTime))
*    say dateTime " " hhmm " " dd 

    iLev =1
      'set vpage 0 8.5 10.6 11'
      'set string 1 c '
      'set strsiz 0.15'
      'draw string 4.25 0.1 Forecast for POLLEN INDEX. Last analysis time: 'AnalTimeString
        

*
*---------------- concentration
*       #xmin xmax ymin ymax
*
      'set vpage 0.25 8.25 5.25 10.75'

      'set grads off'

      'set clevs  1.5 2.5 3.5 4.5'
      'set rbcols 130 131 132 133 134'

      'd POLI'
*      'set annot 1 1'
*      'set font 0'
indPol
      'draw title POLind 'dateTime' UTC'
*      'ccbar Good Fair Moderate Poor VeryPoor'
      'ccbar VeryLow Low Moderate High VeryHigh'

      'set vpage 0.25 8.25 0 5.5'
      'set grads off'

*http://cola.gmu.edu/grads/gadoc/colorcontrol.html
*'set clevs 1.5 2.5 3.5 4.5 5.5'
*  'set clevs 1.5 2.5 3.5 4.5 5.5'
*Reverse order to make the vertical colorbar alphabetical
  'set clevs  -5.5 -4.5 -3.5 -2.5 -1.5'
*  'set rbcols 40 41  3  4  12  34'
* stuff from draw_6_pics
*clr.4 = "60 61 62 63 64 65 66 67 68 69""RAGWEED"
*clr.6 = "70 71 72 73 74 75 76 77 78 79""MUGWORT"
*clr.3 = "40 41 42 43 44 45 46 47 48 49""OLIVE"
*clr.2 = "30 31 32 33 34 35 36 37 38 39""GRASS"
*clr.1 = "50 51 52 53 54 55 56 57 58 59""BIRCH"
*clr.5 = "20 21 22 23 24 25 26 27 28 29""ALDER"
  'set rbcols 65 75 45 35 55 25'



  'd -1*POLISRC'
  'draw title Main contributor to POLind'

  'ccbar Ragweed Mugwort Olive Grass Birch Alder'



*
*---- Printing itself
*
      outfname = outDir'/'substNm_out'_'filesuff''math_format('%03.0f',nbr)

     'printim 'outfname'.png x800 y1000 white'
      say 'Making ' outfname

      'clear'


    time=time+1
    nbr = nbr + 1
*say OOOpppppssss
*break 
  endwhile


*say END

