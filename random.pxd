cimport cython
from libc.stdlib cimport rand, srand
from libc.math cimport log, sqrt, cos, round

cimport numpy as np


cdef extern from "stdlib.h":
    cdef int RAND_MAX # Need to this to normalize appropriately to uniform dist

cdef seed_random(unsigned long long seed)
cdef double exponential_rv(double Lambda)
cdef double uniform_rv()
cdef double normal_rv(double mean, double std)
cdef double gamma_rv(double k, double theta)
cdef double erlang_rv(double k , double theta)
cdef int sample_discrete(int choices, double* data, double Lambda)
cdef double array_sum(double* data, int length)
cdef unsigned binom_rnd(unsigned n, double p)
cdef unsigned binom_rnd_f(double N, double p)
cdef unsigned approx_binom_rnd(unsigned n, double p)
cdef unsigned approx_binom_rnd_f(double n, double p)




