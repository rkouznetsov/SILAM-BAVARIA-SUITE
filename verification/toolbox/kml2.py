from xml.etree import ElementTree as ET
import os, zipfile
from support import xmlpp

def formatdate(when):
    fmt = '%Y-%m-%dT%H:%M:%SZ'
    return when.strftime(fmt)

class KMLObject:
    def __init__(self, KMLName=None, name=None):
        """ KMLObject(KMLname=None, name=None) -> KMLObject
        
        Create a generic KML object corresponding to tag 'name', or
        name of the class if name == None. If given, the object will
        include element 'name' with value 'KMLName'.

        """
        if name is None:
            myname = self.__class__.__name__
            
        else:
            myname = name
        #print 'created ', self.__class__.__name__, 'as', myname
        self._elem = ET.Element(myname)
        
        if KMLName is not None:
            myKMLName = ET.SubElement(self._elem, 'name')
            myKMLName.text = KMLName
            
        self.attrib = self._elem.attrib
            
    def append_simple(self, **kwargs):
        for tag, text in kwargs.items():
            child = ET.SubElement(self._elem, tag)
            child.text = text
            
    def element(self):
        return self._elem
            
    def append(self, what):
        if what is None:
            return
        self._elem.append(what._elem)

    def write_as_tree(self, filename, pretty=True):
        tree = ET.ElementTree(self._elem)
        if pretty:
            try:
                xmlpp.pprint(tree.tostring(), output=filename)
            except AttributeError:
                import StringIO
                string = StringIO.StringIO()
                tree.write(string)
                xmlpp.pprint(string.getvalue(), output=filename)
            #outputfile.close()
        else:
            try:
                tree.write(filename, xml_declaration=True, method='text')
            except TypeError:
                # Older python versions
                tree.write(file, encoding='UTF-8')

class File(KMLObject):
    def __init__(self):
        KMLObject.__init__(self, name='kml')
        self.attrib['xmlns'] = 'http://earth.google.com/kml/2.2'
        
class Document(KMLObject):
    def tofile(self, filename):
        kmlfile = File()
        kmlfile.append(self)
        output = open(filename, 'w')
        kmlfile.write_as_tree(output)
        output.close()

class GroundOverlay(KMLObject):
    def __init__(self, name=None, box=None, tspan=None, icon=None):
        KMLObject.__init__(self, name)
        for x in (box, tspan, icon):
            if x is not None:
                self.append(x)
                
    class LatLonBox(KMLObject):
        def __init__(self, x0, x1, y0, y1):
            KMLObject.__init__(self)
            self.append_simple(north=y1, south=y0, east=x1, west=x0)

class Icon(KMLObject):
        def __init__(self, href):
            KMLObject.__init__(self)
            self.append_simple(href=href)

class TimeSpan(KMLObject):
    def __init__(self, begin, end):
        KMLObject.__init__(self)
        fmt = '%Y-%m-%dT%H:%M:%SZ'
        self.append_simple(begin=begin.strftime(fmt), end=end.strftime(fmt))

class ScreenOverlay(KMLObject):
    def __init__(self, name=None, icon=None, tspan=None):
        KMLObject.__init__(self, name)
        for x in (icon, tspan):
            self.append(x)

    def makeBottomLeft(self):
        overXY = KMLObject(name='overlayXY')
        screenXY = KMLObject(name='screenXY')
        for xy in (overXY, screenXY):
            xy.attrib['x'] = '0'
            xy.attrib['y'] = '1'
            xy.attrib['xunits'] = 'fraction'
            xy.attrib['yunits'] = 'fraction'
        self.append(overXY)
        self.append(screenXY)
        
class Placemark(KMLObject):
    def __init__(self, name, point, description=None):
        KMLObject.__init__(self, name)
        self.append(point)
        self.append_simple(description=description)
    
    class Point(KMLObject):
        def __init__(self, x, y):
            KMLObject.__init__(self)
            self.append_simple(coordinates="%3.6f,%3.6f" % (x,y))

class StyleMap(KMLObject):
    class Pair(KMLObject):
        def __init__(self, key, url):
            KMLObject.__init__(self)
            self.append_simple(key=key)
            self.append_simple(styleUrl=url)
            
    def __init__(self, id):
        KMLObject.__init__(self)
        self.attrib['id'] = id
        
    def append_pair(self, key, styleUrl):
        self.append(StyleMap.Pair(key, styleUrl))

        
def make_marker_color(document, abgr):
    # abgr = (alpha, blue, gree, red), 0...255
    style, style_hl = KMLObject(name='Style'), KMLObject(name='Style')
    hex = '%02x%02x%02x%02x' % abgr
    url_normal = 'cs_%s_norm' % hex
    style.attrib['id'] = url_normal
    iconstyle = KMLObject(name='IconStyle')
    iconstyle.append_simple(color=hex)
    style.append(iconstyle)

    url_hilite = 'cs_%s_hilit' % hex
    style_hl.attrib['id'] = url_hilite
    iconstyle = KMLObject(name='IconStyle')
    iconstyle.append_simple(color=hex)
    iconstyle.append_simple(scale='1.2')
    style_hl.append(iconstyle)

    url_map = 'sm_mrkr_%s' % hex
    stylemap = StyleMap(url_map)
    # The # means that the style url refers to the same file.
    stylemap.append_pair('normal', '#' + url_normal)
    stylemap.append_pair('highlight', '#' + url_hilite)

    for elem in (style, style_hl, stylemap):
        document.append(elem)

    return url_map

class KMZFile:
    def __init__(self, kmlfile, files=[]):
        self.kml = kmlfile
        self.files = files

    def putFiles(self, morefiles):
        self.files.extend(moreFiles)

    def _make_relative_hrefs(self, kmlobj, prefix):
        if kmlobj.tag == 'href':
            url = kmlobj.text
            basename = url.split('/')[-1]
            newurl = '/'.join((prefix, basename))
            kmlobj.text = newurl
            print url, newurl
        if kmlobj.tag == 'description':
            desc = kmlobj.text
            left, url, right = desc.split('"')
            if url.startswith('http') or os.path.sep == '/':
                sep = '/'
            else:
                sep = '\\'
            basename = url.split(sep)[-1]
            newurl = sep.join((prefix, basename))
            kmlobj.text = '"'.join((left, newurl, right))
            print url, newurl
            
        for child in kmlobj.getchildren():
            self._make_relative_hrefs(child, prefix)
        
    def write(self, kmz_filename):
        zf = zipfile.ZipFile(kmz_filename, 'w')
        if kmz_filename.endswith('.kmz'):
            kml_filename = kmz_filename[:-4] + '.kml'
        else:
            kml_filename = kmz_filaneme + '.kml'
        
        if isinstance(self.kml, basestring):
            zf.write(self.kml, os.path.basename(self.kml))
        else:
            # Some magic happens.
            import copy, StringIO
            kml_copy = copy.deepcopy(self.kml)
            self._make_relative_hrefs(kml_copy.element(), 'files')
            iostring = StringIO.StringIO()
            kml_copy.write_as_tree(iostring)
            zf.writestr(os.path.basename(kml_filename), iostring.getvalue())
            iostring.close()
            
        for datafile in self.files:
            datafile_basename = os.path.basename(datafile)
            zf.write(datafile, os.path.join('files', datafile_basename))

        zf.close()
            
    
        
        
        
        
        
