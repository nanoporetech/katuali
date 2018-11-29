#!/usr/bin/env python
import os
import logging

import gpustat

def pick_gpu():
    """
    Get GPU to use from environmental variable SGE_HGR_gpu or find best GPU from gpustat
    """
    logging.basicConfig(format='[%(asctime)s - %(name)s] %(message)s', datefmt='%H:%M:%S', level=logging.INFO)
    logger = logging.getLogger('pick_gpu')
    gpu = os.getenv('SGE_HGR_gpu')
    if gpu is None:
        stats = gpustat.GPUStatCollection.new_query()
        sorter = lambda s: (s.memory_used, s.utilization, s.temperature)
        gpu = sorted(stats.gpus, key=sorter)[0].index
        logger.info('SGE_HGR_gpu was not set, setting GPU to {} based on memory and utilization'.format(gpu))
    else:
        gpu = gpu.replace('cuda', '')
        logger.info('Using gpu {} from SGE_HGR_gpu'.format(gpu))
    return gpu


if __name__ == '__main__':

    print(pick_gpu())
