import cPickle, datetime as dt, time, itertools
import numpy as np
import codecs
from numpy import ma
import copy
import re
import os
import subprocess

try:
    import pylab
    import matplotlib.pyplot as plt
    import matplotlib as mpl
except ImportError:
    print 'Warning: importing pylab/pyplot failed'
    pylab = None

import locale
from toolbox import gradsfile as gf, util, ncreader
from support import pupynere as netcdf

DATE_SMALL = dt.datetime(dt.MINYEAR,1,1)
DATE_BIG = dt.datetime(dt.MAXYEAR,12,31)

"""

A module for handling timeseries.

A timeseries is a set of values indexed with datetime objects. Each
value has a specified duration (of averaging or validity).
Timestamp is always at the END of avraging/validity interval.


A timeseries can be attached to a station with known location, code and
long name. The values can be arbitrary objects. 

The value of timeseries ts at time t is returned by ts[t]. The
associated duration (timedelta) is ts.duration(t). Similar to
dictionaries, the default iterator is over times. The default iterator
yields the times in arbitrary order.

If the values or times are needed in sorted order, use methods data()
and times().

Given a set of stations, timeseries can be initialized from textfiles
(fromFile), or from gridded data files (fromGradsfile). Timeseries can
be exported to text files with similar form (toFile).

"""

class TimeseriesError(Exception):
    pass

class Station:
    def __init__(self, code, name, x, y, height=0.0, area_type='', dominant_source=''):
        self.x = float(x)
        self.y = float(y)
        self.name = name
        self.code = code        
        self.area = area_type.lower()
        self.source = dominant_source.lower()
        self.hgt = height

    def __str__(self):
        return '%s %.2f %.2f %s' % (self.code, self.y, self.x, self.name.encode('utf-8'))
        
    def __hash__(self): 
        return self.code.__hash__()
    
    def __eq__(self, other): 
        return self.code == other.code 

    def __ne__(self, other):
        return not(self == other)
    
    def toFile(self, file, separator=u' '):
        fields = [format%value for format, value in (('%s', self.code), 
                                                     ('%f', self.y),
                                                     ('%f', self.x),
                                                     ('%.1f', self.hgt),
                                                     ('%s', self.name), 
                                                     ('%s', self.area),
                                                     ('%s', self.source))]
        file.write(separator.join(fields).encode('utf-8') + '\n')

def _validIfNotNegative(x):
    return x >= 0.0

def validIfNotNegative(x):
    return x >= 0.0

def getValidIfNotMissing(missing_value, tolerance=1e-9):
    funct = lambda x: np.abs(x-missing_value) > tolerance
    return funct

def alwaysValid(x):
    return True

def _alwaysValid(x):
    return True

def _stringToTimedelta(text):
    number, unit = text.split()
    units = {'hr' : 'hours', 'day' : 'days', 'sec' : 'seconds', 'min' : 'minutes'}
    if not unit in units:
        raise TimeseriesError('Strange time unit: %s' % unit)
    keyword = units[unit]
    return dt.timedelta(**{keyword : float(number)})

class Timeseries:
    def duration(self, time):
        return dt.timedelta(seconds=self._dict[time][1])

    @classmethod
    def fromtuples(cls, times_durations, values, station=None, quantity=None):
        series = cls([], [], [], station, quantity)
        for time_duration, value in itertools.izip(times_durations, values):
            series._dict[dt.datetime(*time_duration[:5])] = (value, time_duration[5])
        return series

    @classmethod
    def fromdict(cls, tsdict, station, quantity):
        series = cls([], [], [], station, quantity)
        series._dict = tsdict
        return series
    
    def __init__(self, values, times, durations, station=None,
                 quantity=None, validator=_validIfNotNegative):
        """Create a timeseries object with the given values and times. The
        argument durations can have several forms:
        - string: '[value] [unit]' where value is float and unit is hr, day, sec or min
        - a timedelta. Like the previous form, this is sets the same duration for each value.
        - a sequence of timedeltas with the same length as values and times.

        The correct range of values can be checked on the fly with the
        validator argument. The validator is a function returning True
        for each valid value. Some common ones are defined above. For
        historical reasons, the default is to ignore negative values.

        When a value is discarded during validation, it is not stored
        and cannot be recovered later. If explicit missing values are
        required, use the alwaysValid validator to allow them in the series."""
        
        if not len(values) == len(times):
            raise TimeseriesError('Values and times are of different size')

        self._dict = {}
        self.const_duration = None

        # Some savings can be had when the duration is constant - check it. BUT this is
        # currently only implemented for some functions for collocating two timeseries.
        if isinstance(durations, basestring):
            self._duration = _stringToTimedelta(durations)
            _duration = _stringToTimedelta(durations)
            seconds = _duration.days*24*3600 + _duration.seconds
            self.const_duration = _duration
            for value, time in zip(values, times):
                if not validator(value):
                    continue
                if not AllowDups and time in self._dict(): 
                    raise TimeseriesError('Duplicate times')
                    
                self._dict[time] = (value, seconds)
        elif isinstance(durations, dt.timedelta):
            self.const_duration = durations
            seconds = durations.days*24*3600 + durations.seconds
            for value, time in zip(values, times):
                if not validator(value):
                    continue
                self._dict[time] = (value, seconds)
        else:
            if len(durations) > 0:
                assert(isinstance(durations[0], dt.timedelta))
                self.const_duration = durations[0]
            else:
                # creating an empty series to be filled by fromdict or fromtuples, 
                # or via insert(). Constant-duration property cannot be exploited.
                # (so far).
                self.const_duration = False
            for value, time, duration in zip(values, times, durations):
                if not validator(value):
                    continue
                self._dict[time] = (value, duration.days*24*3600 + duration.seconds)
                if duration != self.const_duration:
                    self.const_duration = None
        self.quantity = quantity
        self.station = station
        self.validator = validator
                
    def __contains__(self, item):
        return item in self._dict
    
    def __len__(self):
        return len(self._dict)
    
    def __iter__(self):
        return self._dict.__iter__()
    
    def __getitem__(self, key):
        return self._dict[key][0]

    def iteritems(self):
        for time, tpl in self._dict.iteritems():
            yield time,  tpl[0], tpl[1]
    
    def insert(self, time, value, duration):
        if isinstance(duration, int):
            seconds = duration
        else:
            seconds = duration.days*24*3600 + duration.seconds
        if self.validator(value):
            self._dict[time] = (value, seconds)   
        if self.const_duration is not None and self.const_duration != duration:
            self.const_duration = None

    def update(self, series_from):
        """
        Update values of this series with those of series_from, overwriting existing
        values if necessary. The series must have same station and quantity.
        """
        if series_from.quantity != self.quantity:
            raise ValueError('The updated series must have same quantity')
        if series_from.station != self.station:
            raise ValueError('The updated series must have same station')
        for time, value, duration in series_from.iteritems():
            self.insert(time, value, duration)
        
            
    def copy(self):
        return copy.deepcopy(self)

    def times(self, t1=DATE_SMALL, t2=DATE_BIG):
        """Return a sorted array of times of the series, optionally limited to
        a given range."""
        out = [t for t in self._dict if not t < t1 and not t > t2]
        out.sort()
        return out

    def durations(self, t1=DATE_SMALL, t2=DATE_BIG):
        """Return an array of durations, sorted by the time, optionally limited to
        a given range."""
        out = [dt.timedelta(seconds=self._dict[t][1]) for t in self.times(t1, t2)]
        return out

    def data(self, t1=DATE_SMALL, t2=DATE_BIG):
        """Return an array of values, sorted by the time, limited to a given
        range."""
        times = self.times(t1,t2)
        dat = [self[t] for t in times]
        return dat

    def values(self):
        """Like above, but the values are returned without sorting."""
        return [tp[0] for tp in self._dict.values()]
    
    def filled(self, fill_value=np.nan, t1=DATE_SMALL, t2=DATE_BIG):
        """filled_times, filled_values = series.filled()

        Return array of times and values. For each interval of time
        not covered by the series (ie. outside the intervals
        [time-duration, time]) is filled by a missing value and
        duration corresponding to the gap. This is useful for plotting
        the series.
        """
        times = self.times(t1, t2)
        filled_times = []
        values = []
        for ind, time in enumerate(times):
            duration = self.duration(time)
            if ind > 0 and time - duration > times[ind-1]:
                values.append(fill_value)
                filled_times.append(time-duration)
            values.append(self[time])
            filled_times.append(time)
        return filled_times, values
        
    def plot(self, *args, **kwargs): plotSeries3(self, *args, **kwargs)

    def pack(self):
        """Return a binary representation of the series using cPickle. """
        if self.station:
            station = self.station.code
        else:
            station = ''
        times_durations = []
        values = []
        items = self._dict.iteritems()
        for tt, tpl in items:
            times_durations.append((tt.year, tt.month, tt.day, tt.hour, tt.minute, tpl[1]))
            values.append(tpl[0])
        values = np.array(values, dtype=np.float32)
        #times_durations = [ for tt, tpl in items]
        #values = np.array([tpl[0] for tt, tpl in items])
        return cPickle.dumps((times_durations, values, self.quantity, station))
    
    def toFile(self, file, append=False, comment=None, 
               columns=None, quantityAs=None, separator=' ', stdDev=''):
        """Create a MMAS type text file."""
        
        if columns: 
            ind_code = columns['code']
            ind_year = columns['year']
            ind_month = columns['month']
            ind_day = columns['day']
            ind_hour = columns['hour']
            ind_duration = columns['duration']
            ind_value = columns['value']
            ind_quantity = columns['quantity']
            if 'minute' in columns:
                ind_minute = columns['minute']
            else:
                ind_minute = False
            ind_stdDev = ind_value+1
        else:
            ind_code = 0
            ind_quantity = 1
            ind_year = 2
            ind_month = 3
            ind_day = 4
            ind_hour = 5
            ind_duration = 6
            ind_value = 7
            ind_stdDev = 8
            ind_minute = False

        if isinstance(file, str):
            if not append:
                fh = open(file, 'w')
            else:
                fh = open(file, 'a')
        else:
            fh = file
            
        if self.station is not None: 
            code = self.station.code
        else:
            code = 'unknown'

        if columns:
            array = ['***' for i in range(max(columns.values())+2)]
        else:
            array = ['***' for i in range(9)]
	
        if quantityAs is None:
            quantity = self.quantity
        else:
            quantity = quantityAs
		
        times = self.times()
        if comment:
            fh.write(comment.rstrip()+"\n")

        if ind_minute:
            for t in times:
                dt = self.duration(t)
                values = (code, quantity, t.year, t.month, 
                          t.day, t.hour, t.minute, dt.days*24 + dt.seconds / 3600. ,
                          '%g' % self[t],
                          stdDev)
                indices = (ind_code, ind_quantity, ind_year, ind_month, ind_day, 
                           ind_hour, ind_minute, ind_duration, ind_value, ind_stdDev)
                for i in range(len(indices)):
                    array[indices[i]] = str(values[i])
                fh.write(separator.join(array) + "\n")
        else:
            for t in times:
                dt = self.duration(t)
                values = (code, quantity, t.year, t.month, 
                          t.day, t.hour, dt.days*24 + dt.seconds / 3600. ,
                          '%g' % self[t],
                          stdDev)
                indices = (ind_code, ind_quantity, ind_year, ind_month, ind_day, 
                           ind_hour, ind_duration, ind_value, ind_stdDev)
                for i in range(len(indices)):
                    array[indices[i]] = str(values[i])

                fh.write(separator.join(array) + "\n") 
        if isinstance(file, str): fh.close()
    
    def timeAverage(self, window, start=None, values_required=1):
        """
        Compute time averages inside windows defined by start +
        n*window.  More precisely, return series with valid times t(i)
        = start + i*window, i > 0, where the average is taken over
        values with valid times between t(i-1) and t(i), including the
        end of window, excluding the beginning.

        The durations of the series to average are ignored -
        effectively we assume that they are << window.
        """
        if isinstance(window, basestring):
            window = _stringToTimedelta(window)
            
        times = self.times()
        if start is None:
            start = times[0]
        
        aver_times = []
        aver_values = []

        if start > times[-1]:
            raise ValueError('Averaging start is after last timestep')

        # The end of window is inclusive with respect to self.times(), the beginning is
        # exclusive.

        until = start + window
        itertimes = iter(times)
        now = itertimes.next()

        # Initially: we may have until < times[0] - then increment until so that until >= times[0].
        # Or, if until - window >= times[0] - find t so that t > until-window and set now = t.
        while until < times[0]:
            until += window
        while until - window >= now:
            now = itertimes.next()

        # Now start averaging. 
        value = 0.0
        count = 0

        have_values = True
        
        while have_values:
            #print 'now:', now
            while now <= until:
                #print 'averaging:', now, 'until', until
                value += self[now]
                count += 1
                try:
                    now = itertimes.next()
                except StopIteration:
                    have_values = False
                    break
            if count >= values_required:
                #print 'fill', until
                aver_times.append(until)
                aver_values.append(value/count)
            value = 0.0
            count = 0
            until += window

        return Timeseries(aver_values, aver_times, window, self.station, self.quantity)

    def mergeInto(self, target_series, values_required=1):
        values = []
        times = []
        durations = []
        for time, duration, value_target, value_aver in iterMerge(target_series, self, values_required):
            values.append((value_target, value_aver))
            times.append(time)
            durations.append(duration)
        return Timeseries(values, times, durations, self.station, self.quantity, validator=_alwaysValid)

    def apply(self, function):
        """Return a series consisting of f(x) for x in self.values()"""
        dct = {}
        for tt in self:
            val, dur = self._dict[tt]
            dct[tt] = (function(val), dur)
        return Timeseries.fromdict(dct, self.station, self.quantity)

    def timefilter(self, function):
        dct = {}
        for tt in self:
            if not function(tt):
                continue
            dct[tt] = (self._dict[tt])
        return Timeseries.fromdict(dct, self.station, self.quantity)
    
#***************************************************************************************************

def iterMerge(ts_lead, ts_follow, values_required=1):
    """
    Iterate over times of ts_lead, averaging ts_follow into the
    covered windows as needed. At every step, the iteration yields a
    tuple (end_time, duration, value_lead, value_averaged). The rules
    for averaging are similar to the more specific averaging functions
    above.
    
    """

    # itermerge is somewhat expensive. To save time we cover two
    # special cases of synchronous hourly or daily data - then we can
    # use izip3 avoiding all averaging.
    if (ts_lead.const_duration == ts_follow.const_duration == util.one_hour
        and
        ts_lead._dict.iterkeys().next().minute == ts_follow._dict.iterkeys().next().minute == 0):
        # Hourly series with synchronous timesteps.
        for tt, value1, value2 in izip3_timeseries(ts_lead, ts_follow):
            yield tt, ts_lead.const_duration, value1, value2
        return
    elif (ts_lead.const_duration == ts_follow.const_duration == util.one_day
          and
          ts_lead._dict.iterkeys().next().hour == ts_follow._dict.iterkeys().next().hour == 0):

        # Daily series with synchronous timesteps.
        for tt, value1, value2 in izip3_timeseries(ts_lead, ts_follow):
            yield tt, ts_lead.const_duration, value1, value2
        return
    lead_times = ts_lead.times()
    durations = ts_lead.durations()
    follow_times = ts_follow.times()

    #until = lead_times[0]
    #window = durations[0]
    iter_follow_times = iter(follow_times)
    iter_windows = itertools.izip(lead_times, durations)
    now = iter_follow_times.next()

    until, window = iter_windows.next()
    
    while until < follow_times[0]:
        until, window = iter_windows.next()
    while until - window >= now:
        now = iter_follow_times.next()

    # now we should have until - window < now < until. Do we?
    if until - window > now or until < now:
        # no values to iterate over.
        raise StopIteration()

    value = 0.0
    count = 0

    have_values = True

    while have_values:
        while now <= until - window:
            # outside the target windows
            #print 'skipping:', now
            now = iter_follow_times.next()
        while now <= until:
            #print 'averaging:', now, 'until', until, window
            value += ts_follow[now]
            count += 1
            try:
                now = iter_follow_times.next()
            except StopIteration:
                # ran out of values to merge from.
                have_values = False
                break
        if count >= values_required:
            #print 'fill', until
            yield until, window, ts_lead[until], value/count
        value = 0.0
        count = 0
        try:
            old_until = until
            until, window = iter_windows.next()
            if until - window < old_until:
                # The algorithm would fail in this case because we alreay
                # consumed the values.
                print until, window, old_until
                raise ValueError('The windows of ts_lead overlap.')
        except StopIteration:
            # ran out of values to merge into.
            break
    
def izip3_timeseries(ts1, ts2, ignore_duration=False):
    """Iterate over pairs of common values in two timeseries. It is an error to zip two
    timeseries with different durations unless ignore_duration is set to true. The valid
    timestamp is returned among with the values.

    Examples:

    # Find the maximum difference of ts1 and ts2:
    maximum = max(abs(a-b) for time, a, b in izip3_timeseries(ts1, ts2))
    # Use numpy.mean() to calculate the mean difference of ts1, ts2:
    bias = numpy.mean(numpy.array([a-b for time, a, b in izip3_timeseries(ts1, ts2)]))
    
    """
    if (ignore_duration 
        or ts1.const_duration is not None and ts1.const_duration == ts2.const_duration):
        for time in ts1.times():
            if not time in ts2:
                continue
            yield time, ts1[time], ts2[time]
    else:
        for time in ts1.times():
            if not time in ts2:
                continue
            value1, value2, duration = ts1[time], ts2[time], ts1.duration(time)
            if not duration == ts2.duration(time):
                raise TimeseriesError('Identical durations are requested for tszip')
            yield time, value1, value2
            
def izip_timeseries(ts1, ts2, ignore_duration=False):
    """Iterate over pairs of common values in two timeseries. It is an
    error to zip two timeseries with different durations unless ignore_duration is set to true.

    Examples:

    # Find the maximum difference of ts1 and ts2:
    maximum = max(abs(a-b) for a, b in izip_timeseries(ts1, ts2))
    # Use numpy.mean() to calculate the mean difference of ts1, ts2:
    bias = numpy.mean(numpy.array([a-b for a, b in izip_timeseries(ts1, ts2)]))
    
    """
    if ignore_duration or ts1.const_duration is not None and ts1.const_duration == ts2.const_duration:
        for time in ts1:
            if not time in ts2:
                continue
            yield ts1[time], ts2[time]
    else:
        for time in ts1:
            if not time in ts2:
                continue
            value1, value2, duration = ts1[time], ts2[time], ts1.duration(time)
            if not duration == ts2.duration(time):
                raise TimeseriesError('Identical durations are requested for tszip')
            yield value1, value2

izip = izip_timeseries
izip3 = izip3_timeseries

def zip_timeseries(ts1, ts2, ignore_duration=False):
    """
    Create a new Timeseries with the common values of ts1 and
    ts2. That is, zip_timeseries(ts1, ts2)[tt] = (ts1[tt], ts2[tt])
    for each time tt in intersection of ts1.times() and
    ts2.times8). Unless ignore_duration is set to True, every common
    point of ts1 and ts2 must have the same duration.

    """
    zipdict = {}
    if ts1.const_duration is not None and ts1.const_duration == ts2.const_duration:
        ignore_duration = True
    for time in ts1:
        if not time in ts2:
            continue
        value1 = ts1[time]
        value2 = ts2[time]
        duration = ts1.duration(time)
        if not ignore_duration and duration != ts2.duration(time):
            raise TimeseriesError('Identical durations are requested for tszip')
        zipdict[time] = ((value1, value2), duration)

    if ts1.quantity == ts2.quantity:
        quantity = ts1.quantity
    else:
        quantity = None
    if ts1.station == ts2.station:
        station = ts1.station
    else:
        station = None

    zipped = Timeseries.fromdict(zipdict, station, quantity)
    return zipped

#***************************************************************************************************
#
# Routines for sets of timeseries
#
#***************************************************************************************************

def operation(ts1, ts2, op, ignore_duration=False):
    """
    Perform a binary operation on the common values of ts1, ts2, and
    return as a timeseries. The common points must have the same
    duration.

    Arguments:

    ts1, ts2 : Timeseries
    op : two argument function
    ignore_duration (False) : boolean

    If op(x, y) is None, the result is ignored.

    Example:
    
    # compute the difference of ts1 and ts2 as a fraction of their mean:
    def compute_fraction(x,y):
        xpy = x+y
        if abs(xpy) > 0:
            return 2*(x-y)/xpy
        else:
            return None
    fraction = operation(ts1, ts2, compute_fraction)
    """
    if ts1.station == ts2.station:
        station = ts1.station
    else:
        station = None
    if ts1.quantity == ts2.quantity:
        quantity = ts1.quantity
    else:
        quantity = None
    resdict = {}
    for time in ts1:
        if not time in ts2:
            continue
        value1 = ts1[time]
        value2 = ts2[time]
        duration = ts1.duration(time)
        if not duration == ts2.duration(time) and not ignore_duration:
            raise TimeseriesError('Identical durations are requested for timeseries.operation')
        result = op(value1, value2)
        if result is None:
            continue
        resdict[time] = (result, util.dt2sec(duration))

    ts_result = Timeseries.fromdict(resdict, station, quantity)
    return ts_result

#**********************************************************************************
#
# Reading & writing timeseries and stations
#
#**********************************************************************************
    
def fromGradsfile(stations, file, quantity=None, t1=None, t2=None, scale=1.0, 
                  level=0, verbose=False, remove_if_outside=False):
    # do we have many or just one station?
    try: 
        iter(stations)
    except TypeError:
        stations = [stations]
        
    fh = None

    # One or many quantities?
    quantities = ['']
    if quantity:
        if isinstance(quantity, str):
            quantities = [quantity]
        else: 
            quantities = quantity

    nquantities = len(quantities)

    # filename or a handle?
        
    if (isinstance(file, str)):
        desc = gf.GradsDescriptor(file)
        if quantity == None:
            quantity = desc.vars[0]
            quantities = [quantity]
        fh = gf.Gradsfile(desc, quantities, level=level)

    elif isinstance(file, gf.GriddedDataReader):
        if file.nvars() > 1:
            raise ValueError('The grads file object must read only one variable.')
        #if quantity is not None and quantity != file.vars[0]:
        #    raise Exception('The specified quantity must be equal to the grads variable.')
        fh = file
        fh.rewind()

    else:
        raise ValueError('File must be either path to a grads ctl file or a GriddedDataReader object')

    times = np.array(fh.t())
    if t1 is not None:
        times = times[np.logical_not(times < t1)]
    if t2 is not None:
        times = times[np.logical_not(times > t2)]

    if len(times) == 0:
        raise ValueError('Time interval is outside file limits')

    code = None
    ix = None

    nx=len(fh.x())
    ny=len(fh.y())
    
    # check the locations of stations:
    stations_rf = []
    crd = []
    for s in stations:
        stXind,stYind = fh.indices(s.x, s.y)
        if stXind >= nx or stXind < 0 or stYind >= ny or stYind < 0:
            if remove_if_outside:
                if verbose:
                  print 'Station %s is outside the domain x=%f,y=%f ; ix=%f,iy=%f,  '%(s.code,s.x,s.y,stXind,stYind)
                continue
            raise ValueError('Station ' + s.code + ' is outside the domain')
        stations_rf.append(s)
        crd.append((stXind,stYind))

    stations = stations_rf
    nstations = len(stations)
    # create an empty 3d array (or array of lists in fact).
    if quantities:
        d2d = [[[] for q in quantities] for c in stations]
    else:
        d2d = [[[]] for s in stations]
        
    # find the x and y indices (z assumed 0)
    if  len(fh.z()) > 1:
        print (fh.z())
        raise ValueError('The file needs to read 2d fields')

    # Is times big? Do we report the progress?
    if len(times) > 500:
        counter = 0
    else:
        counter = None

    # Create the validator for the timeseries
    if np.isnan(fh.undef()):
        check_missing = lambda x: not np.isnan(x)
    else:
        check_missing = getValidIfNotMissing(fh.undef()*scale, 1e-5)
    
    for t in times:
        if counter and counter%100 != 0 or verbose:
            print 'Reading: ' + str(t)
        fh.goto(t)
        field = fh.read(1)
        field_filled = ma.filled(field, fh.undef())
        if nquantities > 1:
            for i in range(nstations):
                for j in range(nquantities):
                    d2d[i][j].append(field_filled[crd[i] + (j,)] * scale)
        else:
            for i in range(nstations):
                d2d[i][0].append(field_filled[crd[i]] * scale)

    if isinstance(file, str):
        fh.close()
    return [Timeseries(d2d[i][j][:], times, fh.dt(), stations[i], quantity=quantities[j],
                       validator=check_missing)
            for i in range(nstations) for j in range(nquantities)]


def fromGridded(stations, reader, t1=None, t2=None, verbose=False, remove_if_outside=False):
    try:
        # Maybe stations is given as a dictionary?
        stations = stations.values()
    except AttributeError:
        pass
    
    times = np.array(reader.t())
    if t1 is not None:
        times = times[np.logical_not(times < t1)]
    if t2 is not None:
        times = times[np.logical_not(times > t2)]

    if len(times) == 0:
        raise ValueError('Time interval is outside file limits')
    
    x = reader.x()
    y = reader.y()
    
    # check the locations of stations:
    stations_rf = []
    for s in stations:
        if s.x > x[-1] or s.x < x[0] or s.y > y[-1] or s.y < y[0]:
            if remove_if_outside:
                continue
            raise ValueError('Station ' + s.code + ' is outside the domain')
        stations_rf.append(s)

    stations = stations_rf
    nstations = len(stations)
    d2d = [[] for s in stations]

    # find the x and y indices (z assumed 0)
    if  len(reader.z()) > 1:
        raise ValueError('The file needs to read 2d fields')

    crd = [reader.indices(s.x, s.y) for s in stations]

    # Is times big? Do we report the progress?
    if len(times) > 500:
        counter = 0
    else:
        counter = None

    # Create the validator for the timeseries
    check_missing = getValidIfNotMissing(reader.undef(), 1e-5)
        
    for t in times:
        if counter and counter%100 != 0 or verbose:
            print 'Reading: ' + str(t)
        reader.goto(t)
        field = reader.read(1)
        field_filled = ma.filled(field, reader.undef())
        for i in range(nstations):
            d2d[i].append(field_filled[crd[i]])

    return [Timeseries(d2d[i][:], times, reader.dt(), stations[i], validator=check_missing)
            for i in range(nstations)]



#***************************************************************************************************

def fromFile_minutes(*args, **kwargs):
    columns = {'code' : 0, 'quantity' : 1, 
               'year' : 2, 'month' : 3, 
               'day' : 4, 'hour' : 5,
               'minute' : 6, 'duration' : 7,
               'value' : 8}
    kwargs['columns'] = columns
    return fromFile(*args, **kwargs)

def fromGenericFile(stations, file, quantity=None, t1=None, t2=None, validator=_alwaysValid,
             scale=1.0, timeshift=dt.timedelta(0), timePos='end', timefmt="%Y%m%d%H", 
             line_re=r"^(?P<code>\S+)\s+(?P<quantity>\S+)\s+(?P<time>\S+)\s+(?P<value>\S+)\s*$",
             dfl_duration=None, dfl_code=None):
    """
    Uses regexp and timefmt to match the data
    regexp must contain "time" and "value"
    Read timeseries from a text file. By default all values are
    allowed. The values are multiplied by scale, the argument
    timeshift is added to each valid time.

    timePos= (start|end)
    """
    if timePos not in ('start', 'end'):
        raise ValueError('Unknown timePos: %s' % timePos)
    addDuration = timePos == 'start'

    try: # is stations a dictionary?
        stations.keys
        station_dict = stations
    except AttributeError:
        station_dict = dict((station.code, station) for station in stations)

    linere=re.compile(line_re)
    groupdict = linere.groupindex
    if (dfl_code == None) !=  ("code" in groupdict): # aka XOR
      raise ValueError('code should be specified either in regexp or in kwarg')
    if (dfl_code != None):
      stcode = dfl_code

    if (dfl_duration == None) !=  ("duration" in groupdict):
      raise ValueError('duration should be specified either in regexp or in kwarg')
    if (dfl_duration != None):
      duration = dfl_duration

    quantity_select = None
    parse_quantity = False
    if "quantity" in groupdict:
        quantity_select = quantity  #Select the quantity from args
        parse_quantity = True
    elif quantity != None: 
        quantity_read = quantity # Force the quantity from args
    else:
      raise ValueError('Do not know where to take quantity...')

    try: 
        file.read
        reader = file
        close_file = False
    except AttributeError: # maybe file is a string with pathname?
        reader = codecs.open(file, 'r', encoding='utf-8')
        close_file = True
        
    seriesdict = {}
    if quantity:
        for code in station_dict:
            seriesdict[(code, quantity)] = [],[],[]
        
    for line in reader:
        if line.startswith('#') or line.isspace():
            continue

        r = linere.match(line)
        if r:
          if dfl_code == None: # Have the default code
            stcode = r.group("code")
          if not stcode in station_dict:
            continue

          if (parse_quantity):
            quantity_read = r.group("quantity")
            if quantity_select is not None and quantity_read != quantity:
              continue

          if (dfl_duration == None):
              duration = dt.timedelta(hours=float(r.group("duration")))

          time = dt.datetime.strptime(r.group("time"),timefmt) + timeshift

          if t1 and time < t1:
              continue
          if t2 and time > t2:
              continue

          val = float(r.group("value"))*scale

          if addDuration:
              time += duration

          dct_key = (stcode, quantity_read)
          try:
              times, values, durations = seriesdict[dct_key]
              times.append(time)
              values.append(val)
              durations.append(duration)
          except KeyError:
              seriesdict[dct_key] = [time],[val],[dur]

        else: # Failed to parse the line
          print "Trouble with line:",line

    if close_file:
        reader.close()

    return [Timeseries(values, times, durations, station_dict[code], quantity, validator)
            for (code, quantity), (times, values, durations) in seriesdict.iteritems() if values]

def fromFile(stations, file, quantity=None, t1=None, t2=None, columns=None, validator=_alwaysValid,
             scale=1.0, timeshift=dt.timedelta(0), timePos='end', separator=None):
    """
    Read timeseries from a text file. By default all values are
    allowed. The values are multiplied by scale, the argument
    timeshift is added to each valid time.

    timePos= (start|end)
    """
    if timePos not in ('start', 'end'):
        raise ValueError('Unknown timePos: %s' % timePos)
    addDuration = timePos == 'start'

    if columns: 
        ind_code = columns['code']
        ind_year = columns['year']
        ind_month = columns['month']
        ind_day = columns['day']
        ind_hour = columns['hour']
        if 'minute' in columns:
            ind_minute = columns['minute']
        else:
            ind_minute = None
        ind_duration = columns['duration']
        ind_value = columns['value']
        ind_quantity = columns['quantity']
    else:
        ind_code = 0
        ind_quantity = 1
        ind_year = 2
        ind_month = 3
        ind_day = 4
        ind_hour = 5
        ind_duration = 6
        ind_value = 7
	ind_minute = None

    try: # is stations a dictionary?
        stations.keys
        station_dict = stations
    except AttributeError:
        station_dict = dict((station.code, station) for station in stations)

    try: 
        file.read
        reader = file
        close_file = False
    except AttributeError: # maybe file is a string with pathname?
        if file.endswith(".bz2"):
           print "opening bzip2", file
           #reader = codecs.getreader('utf-8')(bz2.BZ2File(file), errors='replace')
#           reader = bz2.BZ2File(file) # Fails for some reason...
           reader= codecs.getreader('utf-8')(subprocess.Popen("bzcat %s"%(file), shell=True, bufsize=1000000, stdout=subprocess.PIPE).stdout)
        elif file.endswith(".gz"):
           print "opening gzip", file
           reader= codecs.getreader('utf-8')(subprocess.Popen("zcat %s"%(file), shell=True, bufsize=1000000, stdout=subprocess.PIPE).stdout)
        else:
          reader = codecs.open(file, 'r', encoding='utf-8')
        close_file = True
        
    #separator = ' '
    seriesdict = {}
    if quantity:
        for code in station_dict:
            seriesdict[(code, quantity)] = [],[],[]
        
    for line in reader.readlines():
        if line.startswith('#') or line.isspace():
            continue
        words = line.split(separator)
        #words = line.split()
        #print line,
        if len(words) == 9 and not columns:
            # Check for mmas files with the '4':
            del(words[2])
        quantity_read = words[ind_quantity]
        if quantity is not None and quantity_read != quantity:
            continue
        code = words[ind_code]
        if not code in station_dict:
            continue
        if ind_minute:
            minutes = int(words[ind_minute])
        else:
            minutes = 0
        #minutes = 0
        if int(words[ind_hour]) == 24:
            time = dt.datetime(int(words[ind_year]), int(words[ind_month]), int(words[ind_day]), 
                               0, minutes) + dt.timedelta(days=1) + timeshift
        else:
            time = dt.datetime(int(words[ind_year]), int(words[ind_month]), int(words[ind_day]), 
                               int(float(words[ind_hour])), minutes) + timeshift
        if t1 and time < t1:
            continue
        if t2 and time > t2:
            continue

        val = float(words[ind_value])*scale
        dur = dt.timedelta(hours=float(words[ind_duration]))
        if addDuration:
            time += dur

        dct_key = (code, quantity_read)
        
        #if dct_key in seriesdict:
        try:
            times, values, durations = seriesdict[dct_key]
            times.append(time)
            values.append(val)
            durations.append(dur)
            #seriesdict[dct_key].insert(time, val, dur)
        #else:
        except KeyError:
            #series = Timeseries((val,), (time,), (dur, ),
            #                    station_dict[code], quantity_read, validator)
            #seriesdict[dct_key] = series
            seriesdict[dct_key] = [time],[val],[dur]
            

    if close_file:
        reader.close()
    #return seriesdict.values()
    #return []
    return [Timeseries(values, times, durations, station_dict[code], quantity, validator)
            for (code, quantity), (times, values, durations) in seriesdict.iteritems() if values]


def series2dict(series):
    # Put a collection of series into a dictionary indexed by 
    # the station.
    dict = {}
    try: 
        iter(series)
    except TypeError:
        series = (series,)
    for s in series:
        if s.station in dict:
            raise TimeseriesError('Trying to make dictionary of multiple series with the same station')
        dict[s.station] = s
        
    return dict

#***************************************************************************************************

def readStations(file, columns=None, separator=None, have_full=False):
    if columns:
        ind_lat = columns['latitude']
        ind_lon = columns['longitude']
        ind_name = columns['name']
        ind_code = columns['code']
        if 'source' in columns:
            ind_source = columns['source']
        else:
            ind_source = None
        if 'area' in columns:
            ind_area = columns['area']
        else:
            ind_area = None
        if 'altitude' in columns:
            ind_alt = columns['altitude']
        else:
            ind_alt = None
    elif separator == ' ':
        ind_lat = 1
        ind_lon = 2
        ind_name = 4
        ind_code = 0
        if have_full:
            ind_alt = 3
            ind_area = -2
            ind_source = -1
        else:
            ind_area = None
            ind_source = None
            ind_alt = None
    else:
        ind_lat = 1
        ind_lon = 2
        ind_name = 4
        ind_code = 0
        ind_area = 5
        ind_source = 6
        ind_alt = 3

    stations = {}
    try:
        file.read
        fh = file
        Fname = None
    except AttributeError:
        #fh = codecs.open(file,'r', encoding='utf-8')  #Breaks for finnish stations
        fh = open(file,'r')
        Fname = file
    linenum = 0

    try:
        for l in fh:
            linenum += 1
            try:
                line=codecs.decode(l,"utf-8")
            except UnicodeDecodeError:
                print Fname,":Non-unicode line:", l
                line=codecs.decode(l,"latin1")
            if line.startswith('#') or line.isspace(): 
                continue
            dat = line.rstrip().split(separator)

            if len(dat) >= 4:
                area = dat[ind_area] if ind_area is not None else ''
                source = dat[ind_source] if ind_source is not None else ''
                alt = float(dat[ind_alt]) if ind_alt is not None else 0.0
                station = Station(dat[ind_code], dat[ind_name], 
                                  dat[ind_lon], dat[ind_lat], alt, area, source)
                stations[dat[ind_code]] = station
            elif (len(dat) == 3):
                # No name given
                stations[dat[ind_code]] = Station(dat[ind_code], '', dat[ind_lon], dat[ind_lat])
            else:
                raise ValueError('Invalid station record: ' + line)

    except Exception, err:
        print 'Parse error at %s:%i : %s' % (file, linenum, err.message)
        raise
        #raise ValueError()
    if fh is not file:
        fh.close()
    return stations
    
#**********************************************************************************
#
# Graphics with matplotlib
#
#**********************************************************************************

def plotSeries3(series, t1=None, t2=None, fill=True, ticks='auto', title=False, scale='linear',
                plotargs={}, plotter=pylab):
    """
    Plot the timeseries using Matplotlib.

    Arguments:

    t1, t2 (None) : restrict into this time range

    fill (True) : fill the series (see the method above) - results in disconnected graph
    where a data point is missing.

    ticks [daily|monthly|auto] (auto) : how to set xticks

    title : generate a title

    plotargs : keyword args for pylab.plot()
    
    """
    if not t1:
        t1 = DATE_SMALL
    if not t2:
        t2 = DATE_BIG

    if fill:
        times, values = series.filled(np.nan, t1, t2)
    else:
        times = series.times(t1, t2)
        values = [series[t] for t in times]
    if len(times) == 0:
        return
    try:
        locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')
    except:
        pass
    if not plotargs:
        plotargs = {'fmt':'-'}
    pylab.plot_date(times, values, **plotargs)
    if ticks == 'daily':
        one_day = dt.timedelta(days=1)
        if times[0].hour == 0:
            first = times[0]
        else:
            first = times[1].replace(hour=0, minute=0)
        last = times[-1].replace(hour=0, minute=0)
        ticks = []
        now = first
        while now <= last:
            ticks.append(now)
            now += one_day
        labels = [time.strftime('%Y/%m/%d') for time in ticks]
    if ticks == 'monthly':
        ticks = []
        year = times[0].year
        while year <= times[-1].year:
            for month in range(12):
                first_day = dt.datetime(year, month, 1)
                if first_day >= times[0] and first_day <= times[-1]:
                    ticks.append(first_day)
        labels = [time.strftime('%Y/%m/%d') for time in ticks]
    if ticks != 'auto':
        pylab.xticks(ticks, labels)
    pylab.grid(1)
    if title:
        pylab.title('%s at %s (%5.2f$^{\circ}$ N %5.2f$^{\circ}$ E)' % 
                    (series.quantity, series.station.name, 
                     series.station.y, series.station.x))
    font = mpl.font_manager.FontProperties(size='small')
    pylab.legend(prop=font)
    ax = mpl.pyplot.gca()
    if scale == 'log':
        try:
            ax.set_yscale('log')
        except:
            pass
    lb = ax.get_xticklabels()
    mpl.pyplot.setp(lb, rotation=30,fontsize=12)
    return times, values

def stations2map(stations, map, marker='bo',markersize=1, label='', text=None):
    x = np.zeros(len(stations), float)
    y = np.zeros(len(stations), float)
    for i in range(len(x)):
        x[i] = stations[i].x
        y[i] = stations[i].y
    X,Y = map(x,y)
    map.plot(X,Y, marker,ms=markersize, label=label)
    if text == 'code':
        for s,xpt,ypt in zip(stations, X, Y):
            plt.text(xpt, ypt, str(s.code))
    elif text == 'name':
        for s,xpt,ypt in zip(stations, X, Y):
            plt.text(xpt, ypt, str(s.name))


def stationScatter(stations, map, data, cmap):
    x = np.zeros(len(stations), float)
    y = np.zeros(len(stations), float)
    for i in range(len(x)-1):
        x[i] = stations[i].x
        y[i] = stations[i].y
    sc = gammaScaler(data,16)
  
    collection = map.scatter(x,y,c=data,marker='o',cmap=cm)
    plt.colorbar(collection)


sec_in_day = 24*3600
max_sec_i32 = 2**31

class NCIO:
    """
    Netcdf-based timeseries IO.

    The timeseries are written or read from a netcdf file some rudimentary cf-like
    conventions. Station code, lon and lat are included for compliance, but no further
    metadata. Currently there is no read support for the station information.

    The series can be of arbitrary shape and extend, but only one quantity is allowed per
    file.

    Writing the series is done in a buffered fashion: the series are added to the NCIO
    object, and written when close() is called.

    For examples see timeseries_tests.

    """

    def __init__(self, file_path, mode, ref_time=None, units=None):
        self.mode = mode
        if mode == 'w':
            self.file_path=file_path
            self.file_path_tmp="%s.pid%d"%(file_path,os.getpid())
            self.ncfile = netcdf.netcdf_file(self.file_path_tmp, 'w')
            self.series_to_write = {}
            self.const_duration = None
            self.closed = False
            self._num_records = 0
            if ref_time:
                self.ref_time = ref_time
            else:
                self.ref_time = dt.datetime(1980,1,1)
            self.quantity = None
            self.units = units
        elif mode == 'r':
            self.series_to_write = None
            self.ncfile = netcdf.netcdf_file(file_path, 'r')
            self.closed = False
            self.ref_time = ncreader.parse_epoch(self.ncfile.variables['time'].units)
        else:
            raise ValueError('Mode must be r or w')
        
    def _init_nc(self):
        ncfile = self.ncfile
        ncfile.createDimension('obs', self._num_records)
        ncfile.createDimension('station', len(self.series_to_write))
        #print 'timelen:', self.__class__._NC_TIME_FMT_LEN
        ncfile.createDimension('station_string', self.__class__._NC_STAT_STR_LEN)

        nclon = ncfile.createVariable('lon', 'f', ('station',))
        nclon.units = 'degrees_east'
        nclon.standard_name = 'longitude'

        nclat = ncfile.createVariable('lat', 'f', ('station',))
        nclat.units = 'degrees_north'
        nclat.standard_name = 'latitude'
        
        nccode = ncfile.createVariable('station_code', 'c', ('station', 'station_string'))
        nccode.cf_role = 'timeseries_id'
        nctime = ncfile.createVariable('time', 'i', ('obs',))
        nctime.units = ncreader.make_epoch(self.ref_time, 'seconds')
        nctime.long_name = 'time of measurement'
        nctime.standard_name = 'time'
        ncval = ncfile.createVariable('value', 'f', ('obs',))
        ncval.standard_name = self.quantity
        ncval.coordinates = 'time lat lon'
        if self.units:
            ncval.units = self.units
        self.ncfile.createVariable('start_index', 'i', ('station',))
        series_size = self.ncfile.createVariable('series_size', 'i', ('station',))
        series_size.long_name = 'number of observations for the stations'
        series_size.sample_dimension = 'obs'
        ncfile.quantity = self.quantity
        ncfile.featureType = 'timeSeries'
        

    def iter_series(self, stations=None, t1=None, t2=None):
        """
        Return an iterator over the series for the given set of stations. If provided,
        only data between t1 and t2 is read.
        """
        
        if self.mode != 'r':
            raise ValueError('Cannot iterate over file in write mode')

        if stations is None:
            stations = self._stations_from_nc()
        try: # is stations a dictionary?
            stations.keys
            station_dict = stations
        except AttributeError:
            station_dict = dict((station.code, station) for station in stations)

        strlen = self.__class__._NC_STAT_STR_LEN
        ncfile = self.ncfile
        ncstat = ncfile.variables['station_code']
        ncval = ncfile.variables['value']
        nctime = ncfile.variables['time']
        start_index = ncfile.variables['start_index']
        series_size = ncfile.variables['series_size']
        #end_index = ncfile.variables['end_index']
        if ncfile.const_duration < 0:
            ncdur = ncfile.variables['duration']
        for statcode, ind_start, length in itertools.izip(ncstat, start_index, series_size):
            statcode_str = statcode.tostring().rstrip()
            if  statcode_str not in station_dict:
                continue
            ind_end = ind_start + length - 1 # inclusive
            nctimes = nctime[ind_start:ind_end+1]
            times_all = [self._time_from_nc(tt) for tt in nctimes]
            ncslice = slice(ind_start, ind_end+1)

            if not t1 and not t2:
                times = times_all
                values = ncval[ncslice]
                if ncfile.const_duration >= 0:
                    duration = dt.timedelta(seconds=int(ncfile.const_duration))
                else:
                    duration = [dt.timedelta(seconds=int(sec)) for sec in ncdur[ncslice]]
                
                yield Timeseries(values, times, duration, station_dict[statcode_str],
                                 ncfile.quantity, alwaysValid)
                
            else:
                mask = np.zeros(ind_end-ind_start+1, dtype=bool)
                times = []
                if t1 and t2:
                    for ind_time, time in enumerate(times_all):
                        if t1 <= time <= t2:
                            mask[ind_time] = True
                            times.append(time)
                elif t1:
                    for ind_time, time in enumerate(times_all):
                        if t1 <= time:
                            mask[ind_time] = True
                            times.append(time)
                elif t2:
                    for ind_time, time in enumerate(times_all):
                        if time <= t2:
                            mask[ind_time] = True
                            times.append(time)

                duration = '1 hr'
                values = ncval[ncslice][mask]
                if ncfile.const_duration >= 0:
                    duration = dt.timedelta(seconds=int(ncfile.const_duration))
                else:
                    duration = [dt.timedelta(seconds=int(sec)) for sec in ncdur[ncslice][mask]]
                yield Timeseries(values, times, duration, station_dict[statcode_str],
                                 ncfile.quantity, alwaysValid)
                    

        
    def add(self, series):
        """
        Schedule a series for writing.
        """
        if self.mode != 'w':
            raise ValueError('NCIO not in write mode')
        if self.closed:
            raise ValueError('File is closed')
        if not self.quantity:
            self.quantity = series.quantity
        elif self.quantity != series.quantity:
            raise ValueError('Only one quantity per file')
        if not (len(series) > 0):
             raise ValueError('Empty series are not allowed')
        self._num_records += len(series)
        self.series_to_write[series.station] = series
        if self.const_duration is None:
            if series.const_duration is not None:
                self.const_duration = series.const_duration
            else:
                self.const_duration = False
        elif self.const_duration != series.const_duration:
            self.const_duration = False
            
    _NC_TIME_FMT = '%Y%m%d%H%M'
    _NC_TIME_FMT_LEN = 4 + 2 + 2 + 2 + 2
    _NC_STAT_STR_LEN = len(_NC_TIME_FMT)

    
    def _nctime(self, time):
        #return '%04i%02i%02i%02i%02i' % (time.year, time.month, time.day, time.hour, time.minute)
        #return cPickle.dumps((time.year, time.month, time.day, time.hour, time.minute))[:12]
        #return '%12i' % (time)
        delta = time - self.ref_time
        seconds = delta.days*sec_in_day + delta.seconds
        if not -max_sec_i32 < seconds < max_sec_i32:
            raise ValueError('Seconds overflow')
        return seconds
        #return '((time.year, time.month, time.day, time.hour, time.minute)))
        #return time.strftime(self._NC_TIME_FMT)

    def _time_from_nc(self, time_int):
        return self.ref_time + dt.timedelta(seconds=int(time_int))
        #string = string.tostring()
        #year, month, day, hour = (int(string[0:4]), int(string[4:6]), int(string[6:8]),
        #                          int(string[8:10]))
        #return dt.datetime(year, month, day, hour)

    def _stations_from_nc(self):
        nclon, nclat = self.ncfile.variables['lon'], self.ncfile.variables['lat']
        nccode = self.ncfile.variables['station_code']
        stations = {}
        for lon, lat, code in itertools.izip(nclon, nclat, nccode):
            codestr = code.tostring().rstrip()
            stations[codestr] = Station(codestr, '', lon, lat)
        return stations
    
    def close(self):
        """
        Close the netcdf file. If the file was opened for writing, the series included
        with the add() method are written to the file.
        """
        
        if self.mode == 'r':
            self.ncfile.close()
            return
        
        self._init_nc()
        self.closed = True
        ncfile = self.ncfile
        start_index = ncfile.variables['start_index']
        #end_index = ncfile.variables['end_index']
        length = ncfile.variables['series_size']
        have_const_dur = self.const_duration is not False
        
        if have_const_dur:
            self.ncfile.const_duration = int(util.dt2sec(self.const_duration))
        else:
            self.ncfile.const_duration = -1

        ind_series = 0
        ind_start = 0

        ncval = ncfile.variables['value']
        nctm = ncfile.variables['time']
        if not have_const_dur:
            ncdur = ncfile.createVariable('duration', 'i', ('obs',))
            ncdur.long_name = 'duration of the observation'
            ncdur.units = 's'
            
        station_code = ncfile.variables['station_code']
        lon, lat = ncfile.variables['lon'], ncfile.variables['lat']
        ind_rec = 0
        for ind_series, (station, series) in enumerate(self.series_to_write.iteritems()):
            len_series = len(series)
            station_code[ind_series] = station.code.encode("utf-8")
            lon[ind_series] = station.x
            lat[ind_series] = station.y
            ind_start = ind_rec
            start_index[ind_series] = ind_start
            length[ind_series] = len_series
            if len_series == 0:
                continue
            if have_const_dur:
                #for ind_time, time in enumerate(series.times()):
                for ind_time, time in enumerate(series):
                    ncval[ind_rec] = series[time]
                    nctm[ind_rec] = self._nctime(time)
                    ind_rec += 1
            else:
                #for ind_time, time in enumerate(series.times()):
                for ind_time, time in enumerate(series):
                    ind_rec = ind_start + ind_time
                    ncval[ind_rec] = series[time]
                    nctm[ind_rec] = self._nctime(time)
                    ncdur[ind_rec] = int(util.dt2sec(series.duration(time)))
                    ind_rec += 1
            #end_index[ind_series] = ind_rec-1

        self.ncfile.close()
        os.rename(self.file_path_tmp, self.file_path)
