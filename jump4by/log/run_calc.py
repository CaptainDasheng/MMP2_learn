#!/usr/bin/env python 
#================================================================
# performing calulcation accoding to the basic setting.
#================================================================
from ..utils.io import * 
import sys

def auto_submit():

    pool = load_pool(pool_path)

    options = load_options()
    
    PreProcess(options)



def calc(pool_path, task_path, overwrite=False):

    pool = load_pool(pool_path)

    pool['functional'].calc(task_path)


if __name__ == '__main__':

    #tpath = sys.argv[1]
    #ppath = sys.argv[2]
    
    calc(submit=True)

    auto_submit()

