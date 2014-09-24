from scipy.stats import norm, poisson
from collections import OrderedDict
from pylab import plt
import pymultinest
import json
import os
import numpy as np
from contextlib import contextmanager
import threading
import sys
import time
import pandas as pd


stdout_lock = threading.Lock()

@contextmanager
def set_stdout_parent(parent):
    """a context manager for setting a particular parent for sys.stdout

    the parent determines the destination cell of output
    """
    save_parent = sys.stdout.parent_header
    with stdout_lock:
        sys.stdout.parent_header = parent
        try:
            yield
        finally:
            # the flush is important, because that's when the parent_header actually has its effect
            sys.stdout.flush()
            sys.stdout.parent_header = save_parent

class ProgressPrinter(pymultinest.ProgressWatcher):
    """
        Continuously writes out the number of live and rejected points.
    """
    def run(self):
        import time
        thread_parent = sys.stdout.parent_header
        while self.running:
            time.sleep(self.interval_ms / 1000.)
            if not self.running:
                break
            try:
                with set_stdout_parent(thread_parent):
                    print(('rejected points: ', len(open(self.rejected, 'r').readlines())))
                    print(('alive points: ', len(open(self.live, 'r').readlines())))
            except Exception as e:
                print(e)



        
