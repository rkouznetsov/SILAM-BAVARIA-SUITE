from toolbox import statistician

import pypar

root = 0

class DistributedStatistician(statistician.DatabaseStatistician):
    def __init__(*args, **kwargs):
        self = args[0]
        statistician.DatabaseStatistician.__init__(*args, **kwargs)
        self.is_master = pypar.rank() == 0
        self.stations_local = None
        #assert pypar.size() > 1
        self.work_locally = not self.is_master
        
    def _get_statistic(self, source, statistic):
        if self.is_master:
            self.work_locally = True
            statistic_by_station = statistician.DatabaseStatistician._get_statistic(self, source, statistic)
            self.work_locally = False
            for sender in range(1, pypar.size()):
                statistic_by_station.update(pypar.receive(sender))
        else:
            if len(self.stations_local) > 0:
                statistic_by_station = statistician.DatabaseStatistician._get_statistic(self, source, statistic)
            else:
                statistic_by_station = {}
            pypar.send(statistic_by_station, root)
        return statistic_by_station
    
    def say(self, what):
        print 'PROC %i:'%pypar.rank(), what
        
    def report(self, model, filename=None, slaves_too=False):
        if not (self.is_master or slaves_too):
            return
        statistician.DatabaseStatistician.report(self, model, filename)

    def _get_stations_local(self):
        all_stations = self._get_all_stations()
        num_stations = len(all_stations)
        mpisize = pypar.size()
        stations_per_proc = num_stations / mpisize
        leftover = num_stations % mpisize
        _station_distr = []
        ind_start = 0
        ind_end = 0
        rank = pypar.rank()
        for ind_proc in range(mpisize):
            ind_start = ind_end
            count = stations_per_proc
            if ind_proc < leftover:
                count += 1
            ind_end = ind_start + count
            _station_distr.append((ind_start, ind_end))
        start_local, end_local = _station_distr[rank]
        
        print pypar.rank(), stations_per_proc, leftover, start_local, end_local
        if start_local >= num_stations:
            self.stations_local = set()
            self.say('No enough stations for this process')
        else:
            self.say('Station range: %i -- %i (%i)' % (start_local, end_local, num_stations))
            self.stations_local = set(all_stations[start_local:end_local])
        assert self.stations_local or not self.is_master
        
    def _get_all_stations(self):
        return statistician.DatabaseStatistician._filtered_obs_stations(self)
        
    def _filtered_obs_stations(self):
        if self.is_master and not self.work_locally:
            return self._get_all_stations()
        if self.stations_local is None:
            self._get_stations_local()
        return self.stations_local
    
    
