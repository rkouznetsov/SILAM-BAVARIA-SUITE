from mpi4py import MPI
import sys
import collections
import time

errs_set = False

def set_errs_fatal():
    if errs_set:
        return
    sys_excepthook = sys.excepthook
    def mpi_excepthook(type, value, traceback):
	# Make MPI abort if one process exits with error
	sys_excepthook(type, value, traceback)
	MPI.COMM_WORLD.Abort(1)
    sys.excepthook = mpi_excepthook

class Error(Exception):
    pass
    
class Dispatcher:
    job_null = None
    msg_ready = 'ready'
    msg_problem = 'problem'

    def _say(self, *args):
        if self.verbose:
            print ' '.join(str(arg) for arg in args)
    
    def __init__(self, work, comm=MPI.COMM_WORLD, chunksize=1, rank_master=0,
                 fail_on_exception=True, verbose=False, sleeptime=0.1):
        self.size = comm.Get_size()
        if self.size < 2:
            raise ValueError('Need at least 2 processes')
        if self.size > len(work) + 1:
            raise ValueError('Too many workers for %i tasks' % len(work))
        self.rank = comm.Get_rank()
        self.rank_master = rank_master
        self.master = self.rank == rank_master
        self.work = list(work)
        self.chunksize = chunksize
        self.comm = comm
        if fail_on_exception:
            set_errs_fatal()
        self.slaves = [rank for rank in range(self.size) if rank != rank_master]
        self.requests = {}
        self.verbose = verbose
        self.sleeptime = sleeptime
        self.status = {rank:'unknown' for rank in self.slaves}
        
    def _send_job(self, slave):
        job = []
        if not self.work:
            job = self.job_null
        else:
            while len(job) < self.chunksize and self.work:
                job.append(self.work.pop())
        self.comm.send(job, slave)
        if job is self.job_null:
            self._say('M: All done, sent null to %i' % (slave))
        else:
            self._say('M: Sent %i tasks to %i, %i left' % (len(job), slave, len(self.work)))
        self.status[slave] = 'busy'

    def _check_ready(self, slave):
        if not slave in self.requests:
            self.requests[slave] = self.comm.irecv(source=slave)
        ready, msg = self.requests[slave].test()
        if ready:
            self.status[slave] = msg
            del(self.requests[slave])
            
    def _finish(self):
        for slave in self.slaves:
            self._send_job(slave) # null job
            self.status[slave] = 'free'
            
    def _dispatch(self):
        while True:
            for slave in self.slaves:
                status = self.status[slave]
                if status == 'ready':
                    if self.work:
                        self._send_job(slave)
                elif status == 'busy' or status == 'unknown':
                    self._check_ready(slave)
                elif status == 'problem':
                    self._say('M: process %i had a problem' % (slave))
                    raise Error()
                else:
                    raise ValueError('Invalid status "%s" during _dispatch()' % status)
            time.sleep(self.sleeptime)
            if not self.work:
                self._finish()
                break
                    
    def _post_ready(self):
        self._say('S:%i, ready' % self.rank)
        self.comm.isend(self.msg_ready, dest=self.rank_master)
                    
    def _get_job(self):
        job = self.comm.recv(source=self.rank_master)
        self._say('S:%i, job received' % (self.rank), job)
        return job
    
    def share(self):
        if self.master:
            self._say('Rank %i, start' % self.rank)
            self._dispatch()
        else:
            self._post_ready()
            job = self._get_job()
            while job is not self.job_null:
                for task in job:
                    yield task
                self._post_ready()
                job = self._get_job()
        self._say('S:%i, done' % self.rank)
        
        # Barrier in the end. Several simultaneous share()'s will lead to weirdness. This
        # might be avoided using tags but why bother.
        self.comm.Barrier()
        
    def problem(self):
        self.comm.isend(self.msg_problem, dest=self.rank_master)
            
        
def demo2(problem=False):
    import random
    import time
    count_tasks = 25
    random.seed('abc')
    sleep_times = [5 * random.random() for task in range(count_tasks)]
    if problem:
        print 'Problem task:', sleep_times[8]
    disp = Dispatcher(sleep_times, chunksize=2, verbose=True, sleeptime=0.01)
    time_start = time.time()
    for task in disp.share():
        print 'S:%i, sleeping %f' % (disp.rank, task)
        if problem and abs(task - sleep_times[8]) < 1e-3:
            print 'S:%i, problem!' % disp.rank
            disp.problem()
        time.sleep(task)
    print 'done', disp.rank
    disp.comm.Barrier()
    time_end = time.time()
    if disp.master:
        print 'Total waiting time:%f' % sum(sleep_times)
    print 'Task %i, took %f seconds' % (disp.rank, time_end-time_start)

def demo1():
    words = 'abc def ghi jkl mno pqr stu vxy'.split()
    for word in share(words, verbose=True):
        print word
    
    
def share(tasks, chunksize=1, verbose=False):
    disp = Dispatcher(tasks, chunksize=chunksize, verbose=verbose)
    for item in disp.share():
        yield item
        
    
if __name__ == '__main__':
    demo1()
    demo2()
    #demo2(problem=True)
