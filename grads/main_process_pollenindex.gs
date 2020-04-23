function drawemep(args)
LastAnalStr = subwrd(args,1)
ncFnm = subwrd(args,2)
outDir = subwrd(args,3)
hroff = subwrd(args,4)
mpdset = subwrd(args,5)


say ncFnm
say outDir
drawScr_POLI = "draw_POLI_nc4.gs"
*
*      ... and draw it
*



substNm = "POLI"

* These do not matter. Just to keep the same interface
factor = "1"
massUnit = "1"
scale_legend = "10.0 log 0"
colour = 'def_lowwhite' 
say 'run 'drawScr_POLI' 'ncFnm' 'substNm' 'substNm' 'factor' 'massUnit' 'LastAnalStr' 'colour' 'scale_legend' 'outDir' 'hroff' 'mpdset
'run 'drawScr_POLI' 'ncFnm' 'substNm' 'substNm' 'factor' 'massUnit' 'LastAnalStr' 'colour' 'scale_legend' 'outDir' 'hroff' 'mpdset
'quit'

'quit'

