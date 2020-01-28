import tables
import os, datetime as dt
from os import path

class Tickdata(tables.IsDescription):
    timestamp = tables.StringCol(19)
    value = tables.Float32Col()
    duration_days = tables.Int32Col()
    duration_seconds = tables.Int32Col()


class Pointdata(tables.IsDescription):
    timestamp = tables.StringCol(19)
    value = tables.Float32Col()
    duration_days = tables.Int32Col()
    duration_seconds = tables.Int32Col()
    #station = tables.StringCol(32)


class PointdataIndex(tables.IsDescription):
    station = tables.StringCol(32)
    ind_start = tables.Int64Col()
    len = tables.Int32Col()
    #ind_end = tables.Int64Col()
    
def series2dict_by_quantity(series_list):
    seriesdict = {}
    for series in series_list:
        #print 'series2dict: ', series.station
        if not series.quantity in seriesdict:
            seriesdict[series.quantity] = [series]
        else:
            seriesdict[series.quantity].append(series)
    return seriesdict

def _series2hdf(series, h5file, group):
    if not series.station:
        code = '_undefined'
    else:
        code = '_' + series.station.code
    #print series.station
    table = h5file.createTable(group, code, Tickdata, 'Timeseries at %s' % code)
    times = series.times()
    recordptr = table.row
    for time in times:
        recordptr['timestamp'] = time.isoformat()
        recordptr['value'] = series[time]
        duration = series.duration(time)
        recordptr['duration_days'] = duration.days
        recordptr['duration_seconds'] = duration.seconds
        recordptr.append()
    table.flush()

def _series2hdf_indexed(series_list, h5file, group):
    series_dict = timeseries.series2dict(series_list)
    stations = sorted(series_dict.keys(), key=lambda x: x.code)
    table = h5file.createTable(group, 'data', Pointdata, 'data')
    index = h5file.createTable(group, 'index', PointdataIndex, 'index')
    row = table.row
    indexrow = index.row
    ind = 0
    for station in stations:
        indexrow['station'] = station.code
        indexrow['ind_start'] = ind
        series = series_dict[station]
        times = series.times()
        for time in times:
            row['timestamp'] = time.isoformat()
            row['value'] = series[time]
            duration = series.duration(time)
            row['duration_days'] = duration.days
            row['duration_seconds'] = duration.seconds
            #row['station'] = station.code
            row.append()
            ind += 1
        indexrow['len'] = ind - indexrow['ind_start'] + 1
        indexrow.append()
    table.flush()
    index.flush()
    
def series2hdf_indexed(series_in, filename):
    try:
        series_in.keys
        seriesdict = series_in
    except AttributeError:
        seriesdict = series2dict_by_quantity(series_in)

    if path.exists(filename):
        os.unlink(filename)
    h5file = tables.openFile(filename, mode='w', title='tsdata')
    for quantity in seriesdict:
        if quantity is None:
            quant = 'unknown'
        else:
            quant = quantity
        group = h5file.createGroup('/', quant, quant) 
        _series2hdf_indexed(seriesdict[quantity], h5file, group)
        #for series in seriesdict[quantity]:
        #    _series2hdf(series, h5file, group)
    h5file.close()
        
    
def series2hdf(series_in, filename):
    try:
        series_in.keys
        seriesdict = series_in
    except AttributeError:
        seriesdict = series2dict_by_quantity(series_in)

    if path.exists(filename):
        os.unlink(filename)
    h5file = tables.openFile(filename, mode='w', title='tsdata')
    for quantity in seriesdict:
        if quantity is None:
            quant = 'unknown'
        else:
            quant = quantity
        group = h5file.createGroup('/', quant, quant) 
        for series in seriesdict[quantity]:
            _series2hdf(series, h5file, group)
    h5file.close()

def _series4quantity2(h5file, stations, quantity):
    group = h5file.getNode(h5file.root, quantity)
    series_list = []
    for station in stations:
        code = '_' + station.code
        if not code in group:
            continue
        table = group._f_getChild(code)
        nrows = table.nrows
        series = timeseries.Timeseries([], [], [], station, quantity)
        for row in table.iterrows():
            timestamp = row['timestamp']
            time = dt.datetime(int(timestamp[0:4]), int(timestamp[5:7]), int(timestamp[8:10]),
                               int(timestamp[11:13]), int(timestamp[14:16]), int(timestamp[17:19]))
            duration = dt.timedelta(days=row['duration_days'], seconds=row['duration_seconds'])
            series._insert(time, row['value'], duration)
        if len(series) > 0:
            series_list.append(series)
    return series_list
    
def _series4quantity(h5file, stations, quantity):
    group = h5file.getNode(h5file.root, quantity)
    series_list = []
    for station in stations:
        code = '_' + station.code
        if not code in group:
            continue
        table = group._f_getChild(code)
        nrows = table.nrows
        times, durations, values = [], [], []
        for row in table.iterrows():
            timestamp = row['timestamp']
            time = dt.datetime(int(timestamp[0:4]), int(timestamp[5:7]), int(timestamp[8:10]),
                               int(timestamp[11:13]), int(timestamp[14:16]), int(timestamp[17:19]))
            times.append(time)
            durations.append(dt.timedelta(days=row['duration_days'], seconds=row['duration_seconds']))
            values.append(row['value'])
        series = timeseries.Timeseries(values, times, durations, station, quantity)
        series_list.append(series)
    return series_list

def _series4times(h5file, stations, quantity, time_start=None, time_end=None):
    if time_start and time_end:
        condition = '(timestamp >= "%s") & (timestamp <= "%s")' % (time_start.isoformat(),
                                                                   time_end.isoformat())
    elif time_start:
        condition = '(timestamp >= "%s")' % time_start.isoformat()
    elif time_end:
        condition = '(timestamp <= "%s")' % time_end.isoformat()
    else:
        condition = None
        
    group = h5file.getNode(h5file.root, quantity)
    series_list = []
    
    for station in stations:
        code = '_' + station.code
        if not code in group:
            continue
        table = group._f_getChild(code)
        nrows = table.nrows
        times, durations, values = [], [], []
        if condition:
            iterator = table.where(condition)
        else:
            iterator = table.iterrows()
        for row in iterator:
            timestamp = row['timestamp']
            
            time = dt.datetime(int(timestamp[0:4]), int(timestamp[5:7]), int(timestamp[8:10]),
                               int(timestamp[11:13]), int(timestamp[14:16]), int(timestamp[17:19]))
            times.append(time)
            durations.append(dt.timedelta(days=row['duration_days'], seconds=row['duration_seconds']))
            values.append(row['value'])
        series = timeseries.Timeseries(values, times, durations, station, quantity)
        #if code == '_1304':
        #    print times[0], values[0]
        #    print series.times()[0]
        
        series_list.append(series)
    return series_list


def _quantities(h5file):
    return list(h5file.iterNodes(h5file.root))

def fromhdf(stations, filename, quantity=None, time_start=None, time_end=None):
    h5file = tables.openFile(filename, 'r')

    if not quantity:
        quantities = _quantities(h5file)
        series_list = []
        for quantity in quantities:
            series_list.extend(_series4times(h5file, stations, quantity, time_start, time_end))
    else:
        series_list = _series4times(h5file, stations, quantity, time_start, time_end)
    h5file.close()
    return series_list

def _series4times_indexed(h5file, stations, quantity, time_start=None, time_end=None):
    if time_start and time_end:
        condition = '(timestamp >= "%s") & (timestamp <= "%s")' % (time_start.isoformat(),
                                                                   time_end.isoformat())
    elif time_start:
        condition = '(timestamp >= "%s")' % time_start.isoformat()
    elif time_end:
        condition = '(timestamp <= "%s")' % time_end.isoformat()
    else:
        condition = None
        
    group = h5file.getNode(h5file.root, quantity)
    series_list = []

    # load the index
    index = group._f_getChild('index')
    index_dict = {}
    for row in index:
        index_dict[row['station']] = (row['ind_start'], row['ind_start'] + row['len'])
    stations = sorted(stations, key=lambda x: x.code)

    table = group._f_getChild('data')
    ind_last = len(table)
    for station in stations:
        if not station.code in index_dict:
            #print '%s not in index' % station
            #raise timeseries.TimeseriesError()
            continue
        ind_start, ind_end = index_dict[station.code]
        times, durations, values = [], [], []
        if condition:
            iterator = table.where(condition, start=ind_start, stop=ind_end)
        else:
            iterator = table.iterrows(start=ind_start, stop=ind_end)
        for row in iterator:
            #if row['station'] != station.code:
            #    print row['station'], station.code, ind_start, ind_end
            timestamp = row['timestamp']
            time = dt.datetime(int(timestamp[0:4]), int(timestamp[5:7]), int(timestamp[8:10]),
                               int(timestamp[11:13]), int(timestamp[14:16]), int(timestamp[17:19]))
            times.append(time)
            durations.append(dt.timedelta(days=row['duration_days'], seconds=row['duration_seconds']))
            values.append(row['value'])
        series = timeseries.Timeseries(values, times, durations, station, quantity)
        series_list.append(series)
    return series_list

def _series4times_indexed2(h5file, stations, quantity, time_start=None, time_end=None, validator=None):
    if time_start and time_end:
        condition = '(timestamp >= "%s") & (timestamp <= "%s")' % (time_start.isoformat(),
                                                                   time_end.isoformat())
    elif time_start:
        condition = '(timestamp >= "%s")' % time_start.isoformat()
    elif time_end:
        condition = '(timestamp <= "%s")' % time_end.isoformat()
    else:
        condition = None
        
    group = h5file.getNode(h5file.root, quantity)
    series_list = []

    # load the index
    index = group._f_getChild('index')
    index_dict = {}
    for row in index:
        index_dict[row['station']] = (row['ind_start'], row['ind_start'] + row['len'])
    stations = sorted(stations, key=lambda x: x.code)

    table = group._f_getChild('data')
    ind_last = len(table)
    for station in stations:
        if not station.code in index_dict:
            #print '%s not in index' % station
            #raise timeseries.TimeseriesError()
            continue
        ind_start, ind_end = index_dict[station.code]
        # print station.code, ind_start, ind_end
        #times, durations, values = [], [], []
        values = {}
        durations = {}
        if condition:
            iterator = table.where(condition, start=ind_start, stop=ind_end)
        else:
            iterator = table.iterrows(start=ind_start, stop=ind_end)
        if validator:
            for row in iterator:
                #if row['station'] != station.code:
                #    print row['station'], station.code, ind_start, ind_end
                value = row['value']
                if not validator(value):
                    continue
                timestamp = row['timestamp']
                time = dt.datetime(int(timestamp[0:4]), int(timestamp[5:7]), int(timestamp[8:10]),
                                   int(timestamp[11:13]), int(timestamp[14:16]), int(timestamp[17:19]))

                durations[time] = dt.timedelta(days=row['duration_days'], seconds=row['duration_seconds'])
                values[time] = value
        else:
            for row in iterator:
                timestamp = row['timestamp']
                time = dt.datetime(int(timestamp[0:4]), int(timestamp[5:7]), int(timestamp[8:10]),
                                   int(timestamp[11:13]), int(timestamp[14:16]), int(timestamp[17:19]))

                durations[time] = dt.timedelta(days=row['duration_days'], seconds=row['duration_seconds'])
                values[time] = row['value']

        series = timeseries.Timeseries.fromdict(values, durations, station, quantity)
        series_list.append(series)
    return series_list



def fromhdf_indexed(stations, filename, quantity=None, time_start=None, time_end=None, validator=None):
    h5file = tables.openFile(filename, 'r')

    if not quantity:
        quantities = _quantities(h5file)
        series_list = []
        for quantity in quantities:
            series_list.extend(_series4times_indexed(h5file, stations, quantity,
                                                     time_start, time_end, validator))
    else:
        series_list = _series4times_indexed2(h5file, stations, quantity,
                                             time_start, time_end, validator)
    h5file.close()
    return series_list
    

if __name__ == '__main__':
    from toolbox import timeseries
    import time, datetime as dt, sys
    #stationfile = '/home/vira/silam/measurements/airbase/2006_Q1/stations_austria.dat'
    #seriesfile = '/home/vira/silam/measurements/airbase/2006_Q1/cnc_so2_austria.dat'

    variable = 'cnc_pm10'
    stationfile = '/home/vira/tmp/stations_%s_all.dat' % variable
    seriesfile = '/home/vira/tmp/%s_all.dat' % variable
    
    t1 = dt.datetime(2009, 2, 1)
    t2 = dt.datetime(2009, 2, 27)
    t1=t2=None
    stations = timeseries.readStations(stationfile)
    seconds = time.clock()
    series_list = timeseries.fromFile(stations.values(), seriesfile, variable, t1=t1, t2=t2)
    print 'Reading from ascii: %f' % (time.clock() - seconds)
    
    print 'Done reading'
    seconds = time.clock()
    series2hdf(series_list, 'test.hdf')
    print 'Writing hdf: %f' % (time.clock() - seconds)
    print 'Done writing'
    sys.exit()
    seconds = time.clock()
    series2hdf_indexed(series_list, 'test_i.hdf')
    print 'Writing hdf with index: %f' % (time.clock() - seconds)

    seconds = time.clock()
    series_list_from_ind_hdf = fromhdf_indexed(stations.values(),
                                               'test_i.hdf', variable,
                                               time_start=t1, time_end=t2)
    print 'Reading hdf (indexed): %f' % (time.clock() - seconds)
    print len(series_list_from_ind_hdf), len(series_list_from_ind_hdf[0])

    
    seconds = time.clock()
    outfile = open('test.dat', 'w')
    for series in series_list:
        series.toFile(outfile)
    outfile.close()
    print 'Writing ascii: %f' % (time.clock() - seconds)
 
    seconds = time.clock()
    series_list_from_hdf = fromhdf(stations.values(), 'test.hdf', variable, time_start=t1, time_end=t2)
    print 'Reading hdf: %f' % (time.clock() - seconds)

    
    dictascii = timeseries.series2dict(series_list)
    dicthdf = timeseries.series2dict(series_list_from_hdf)
    dicthdf_ind = timeseries.series2dict(series_list_from_ind_hdf)
    for station in dictascii:
        if not station in dicthdf:
            print station
        
    for station in dictascii:
        if not station in dicthdf_ind:
            print station
