import gradsfile
from gradsfile import Gradsfile
import numpy as np
from numpy import array as ar
import kml2 as kml
import matplotlib as mpl
from matplotlib.backends import backend_agg
import figs
from os import system, tmpnam, unlink

def make_loader(urls_to_link, kml_out, expires):
    doc = kml.Document()
    for url in urls_to_link:
        link = kml.KMLObject(name='Link')
        link.add('href', url)
        nwlink = kml.KMLObject(name='NetworkLink')
        nwlink.adopt(link)
        doc.adopt(nwlink)
    
    doc.toFile(kml_out)

class Linkset:
    def __init__(self, output_kml_file, name):
        self.kmlout = output_kml_file
        self.name = name
        self.links = []
        
    def addLink(self, url, name):
        link = kml.KMLObject(name='Link')
        link.add('href', url)
        nwlink = kml.KMLObject(name, 'NetworkLink')
        nwlink.adopt(link)
        self.links.append(nwlink)

    def create(self):
        if not self.links:
            return
        doc = kml.Document(self.name)
        for link in self.links:
            doc.adopt(link)
        doc.toFile(self.kmlout)
    
class GoogleEarthFig:
    def __init__(self, ctlFile, variable, figTemplate, 
                 KMLTemplate, figTemplateUrl=None):
        self.figTemplate = figTemplate

        if (figTemplateUrl is None):
            self.figTemplateUrl = figTemplate
        else:
            self.figTemplateUrl = figTemplateUrl

        self.KMLOutFile = KMLTemplate
        self.ctlFile = ctlFile
        self.variable = variable
        self.t0 = None
        self.t1 = None
        self.t = None
        self.iz = 0
        self.cbarFile = None
        self.tmpfile = None
        self.cm = mpl.cm.jet
        self.conversion = 1.0
        self.opacity = None
        self.title = None
        self.dpi = 72
        self.expires = None
        self.jpegQuality = 70
        self.processed = False
        self.files_created = []
        
    def setTimeRange(self, t0=None, t1=None):
        self.t0 = t0
        self.t1 = t1
        
    def setFigTemplateUrl(self, url):
        self.figTemplateUrl = url
        
    def setTime(self, t):
        self.t = t
    
    def setKMLFile(self, file):
        self.KMLOutFile = file
    
    def setColorbarOutFile(self, file):
        self.cbarFile = file
        self.cbarUrl = file

    def setColorbarUrl(self, url):
        self.cbarUrl = url

    def setTmpfile(self, file):
        self.tmpfile = file

    def setColormap(self, cm):
        self.cm = cm

    def setConversion(self, c):
        self.conversion = c
        
    def setOpacity(self, q):
        self.opacity = q
        
    def setTitle(self, title):
        self.title = title

    def setDpi(self, dpi):
        self.dpi = dpi

    def setExpiry(self, expires):
        self.expires = expires

    def setQuality(self, quality):
        self.jpegQuality = quality
        
    def process(self):
        if isinstance(self.ctlFile, gradsfile.GriddedDataReader):
            gf = self.ctlFile
        else:
            gf = Gradsfile(self.ctlFile, self.variable)

        kmldoc = kml.Document()

        times = ar(gf.t())
        if self.t0: times = times[times > self.t0]
        if self.t1: times = times[times < self.t1]
        if self.t: times = [self.t]

        # check the tmpfile
        if self.tmpfile is None:
            tmpfile = tmpnam()+'.png'
        else:
            tmpfile = self.tmpfile
            
        count = 0

        for t in times:
            count += 1
            gf.goto(t)
            dat = gf.read(1)

            if dat.ndim > 2:
                field = dat[:,:,self.iz]

            else:
                if self.iz > 0:
                    gf.close()
                    raise Exception("Figure is requested on a specific level, but the data are 2D")
                field = dat[:,:]

            field = self.conversion*field
            outFile = t.strftime(self.figTemplate)
            outFile = outFile.replace('#','%i'%count)
            fig = figs.contourf_raw(gf.x(), gf.y(), field, figs.clevs(self.variable), 
                                    figs.ccols(self.variable, baseCM=self.cm))

            cv = backend_agg.FigureCanvasAgg(fig)
            cv.print_figure(tmpfile, dpi=self.dpi)

            if outFile.endswith('jpg'):
                arguments = '-shave 1x0 -quality %i%%' % self.jpegQuality
            else:
                arguments = '-shave 1x0'
            command = 'convert %s %s %s' % (tmpfile, arguments, outFile)
            print command
            #system('convert '+tmpfile+' -shave 1x0 '+outFile)
            system(command)
            self.files_created.append(outFile)
            
            x = gf.x()
            y = gf.y()
            box = kml.GroundOverlay.LatLonBox(x[0], x[-1], y[0], y[-1])
            icon = kml.Icon(t.strftime(self.figTemplateUrl.replace('#', '%i'%count)))
            tspan = kml.TimeSpan(t, t+gf.dt())
            over = kml.GroundOverlay(self.variable, box, tspan, icon)

            if self.opacity is not None:
                over.add('color','%xffffff' % int(self.opacity*255))

            kmldoc.adopt(over)
            
        gf.close()

        if self.cbarFile:
            cbfig = figs.colorbarFig(self.variable, self.cm, text=self.title)
            cv = backend_agg.FigureCanvasAgg(cbfig)
            cv.print_figure(self.cbarFile, dpi=60)
            cbar = kml.ScreenOverlay('colorbar', kml.Icon(self.cbarUrl))
            cbar.add('visibility', '0')

            if self.opacity is not None:
                cbar.add('color','%xffffff' % int(self.opacity*255))

            cbar.makeBottomLeft()
            kmldoc.adopt(cbar)

        kmlfile = kml.KMLFile()
        kmlfile.adopt(kmldoc)

        if self.expires:
            control = kml.KMLObject(name='NetworkLinkControl')
            control.add('expires', kml.dateTime(self.expires))
            kmlfile.adopt(control)

        #control.add('expires', kml.dateTime(expires))
        kmlfile.toFile(self.KMLOutFile)
        unlink(tmpfile)

        self.kmlfile = kmlfile
        self.processed = True


    def makeKMZ(self, kmzfile):
        if not self.processed:
            self.process()
        kmz = kml.KMZFile(self.kmlfile, self.files_created)
        kmz.write(kmzfile)
        
    
def station2placemark(station, name, fig_url=None):
    if fig_url:
        description = "<img src=%s />\n"%figUrl
    else:
        description = ''
    return kml.Placemark(name, kml.Placemark.Point(station.x, station.y), description)

def station_figs_2_kml(stations, figures, outputfilename, url='', use_colors=True):
    kmldoc = kml.Document()

    if use_colors:
        blue = kml.make_marker_color(kmldoc, (255, 127, 255, 255))
        red = kml.make_marker_color(kmldoc, (255, 0, 200, 0))
        green = kml.make_marker_color(kmldoc, (255, 0, 200, 0))
        colormap = {'rural' : green, 'suburban' : blue, 'urban' : red}
        
        
    for station, filename in zip(stations, figures):
        description = '<img src="%s%s" />\n' % (url, filename)
        #placemk = kml.Placemark(station.code,
        placemk = kml.Placemark('',
                                kml.Placemark.Point(station.x, station.y),
                                description)
        if use_colors and station.area in colormap:
            placemk.append_simple(styleUrl=colormap[station.area])
        kmldoc.append(placemk)

    kmlfile = kml.File()
    kmlfile.append(kmldoc)
    if outputfilename.endswith('kmz'):
        kmz = kml.KMZFile(kmlfile, figures)
        kmz.write(outputfilename)
    else:
        kmldoc.tofile(outputfilename)
