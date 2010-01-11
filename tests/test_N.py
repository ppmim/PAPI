import os
import sys
import time
import threading


class ReduceThread(threading.Thread):

    """
    Main threaded class to execute the call to the remote object for reduction
    """

    #def __init__(self, resource_id, source_frame, dark_frame, flat_frame, out_frame ):
    #    threading.Thread.__init__(self)

    def __init__(self, id):
        threading.Thread.__init__(self)
        self.id  = id

    def run(self):
        print "Soy el hilo", self.id
        command='/disk-a/caha/panic/DEVELOP/PIPELINE/PANIC/tests/run_simple_%d.sh' %self.id
        os.system(command)


if __name__=="__main__":

    print "COMIENZA !!!!"

    start=time.time()
    threads=[]
    for i in range(1, 5):
        threads.append( ReduceThread(i) )
        threads[i-1].start()

    for t in threads:
        t.join()

    print "Finalizaron todas las HEBRAS en %f secs !!!" %(time.time()-start)
