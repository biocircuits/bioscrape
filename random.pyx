# cython: boundscheck=False
# cython: cdivision=True
# cython: wraparound=False

cimport cython
from libc.stdlib cimport rand, srand
from libc.math cimport log, sqrt, cos, round

cimport numpy as np
import time

cdef extern from "stdlib.h":
    cdef int RAND_MAX # Need to this to normalize appropriately to uniform dist

# MT Stuff

cdef unsigned NN = 312
cdef unsigned MM = 156
cdef unsigned long long MATRIX_A = 0xB5026F5AA96619E9ULL
cdef unsigned long long UM = 0xFFFFFFFF80000000ULL
cdef unsigned long long LM = 0x7FFFFFFFULL
cdef unsigned long long mt[312]
cdef unsigned mti = NN + 1
cdef unsigned long long mag01[2]

cdef mt_seed(unsigned long long seed):
    global mt
    global mti
    global mag01
    global NN
    global MATRIX_A
    mt[0] = seed
    for mti in range(1,NN):
        mt[mti] = (6364136223846793005ULL * (mt[mti-1] ^ (mt[mti-1] >> 62)) + mti)

    mag01[0] = 0ULL
    mag01[1] = MATRIX_A
    mti = NN


cdef unsigned long long genrand64():
    cdef int i
    cdef unsigned long long x
    global mag01
    global mti
    global mt
    global NN
    global MM
    global UM
    global LM

    if mti >= NN:
        for i in range(NN-MM):
            x = (mt[i]&UM) | (mt[i+1]&LM)
            mt[i] = mt[i+MM] ^ (x>>1) ^ mag01[int(x&1ULL)]

        for i in range(NN-MM, NN-1):
            x = (mt[i]&UM)|(mt[i+1]&LM)
            mt[i] = mt[i+(MM-NN)] ^ (x>>1) ^ mag01[int(x&1ULL)]

        x = (mt[NN-1]&UM)|(mt[0]&LM)
        mt[NN-1] = mt[MM-1] ^ (x>>1) ^ mag01[int(x&1ULL)]
        mti = 0

    x = mt[mti]
    mti += 1
    x ^= (x >> 29) & 0x5555555555555555ULL
    x ^= (x << 17) & 0x71D67FFFEDA60000ULL
    x ^= (x << 37) & 0xFFF7EEE000000000ULL
    x ^= (x >> 43);

    return x

def py_rand_int():
    return genrand64()

# Functions

# Seed the random number generator
cdef seed_random(unsigned long long seed):
    """
    Seed the C random number generator with the current system time.
    :return: none
    """
    if seed == 0:
        mt_seed(time.time())
    else:
        mt_seed(seed)

def py_seed_random(unsigned long long seed = 0):
    seed_random(seed)



cdef double exponential_rv(double Lambda):
    """
    Generate an exponentially distributed random number
    :param Lambda: (double) the rate parameter of the distribution
    :return: (double) a randomly distributed number
    """
    return -1.0/Lambda* log( uniform_rv() )

def py_exponential_rv(double Lambda):
    return exponential_rv(Lambda)



cdef double uniform_rv():
    """
    Generate a uniform random variable in [0,1]
    :return: (double) a random uniform number in [0,1]
    """
    return (genrand64() >> 11) * (1.0/9007199254740991.0)

def py_uniform_rv():
    return uniform_rv()





cdef double normal_rv(double mean, double std):
    """
    Generate a normal random variable
    :param mean: (double) the mean
    :param std: (double) the standard deviation
    :return:
    """
    cdef double u, v, R, theta
    u = uniform_rv()
    v = uniform_rv()
    R = sqrt(-2*log(u))
    theta = 2*3.141592653589793238462643383279502884*v
    return R*cos(theta)*std + mean

def py_normal_rv(double mean, double std):
    return normal_rv(mean,std)




cdef double gamma_rv(double k, double theta):
    """
    Generate a gamma random variable with parameters (k, theta). The mean is k*theta. Variance k*theta^2. Can think of
     this as adding together k independent exponentials each with mean theta.
    :param k: (double) the first parameter of the gamma distribution. if k is an integer, it is the number of
            independent exponentially distributed steps.
    :param theta: (double) the second parameter. theta is the mean of each exponential step.
    :return: (double) a randomly generated gamma random variable.
    """
    cdef double d, c, v, x, UNI
    d = k - 1.0/3
    c = 1/sqrt(9.0*d)

    while True:
        x = normal_rv(0,1)
        v = (1+c*x)**3
        UNI = uniform_rv()
        if v > 0 and log(UNI) < 0.5*x**2+d-d*v+d*log(v):
            return d*v*theta


def py_gamma_rv(double k, double theta):
    return gamma_rv(k,theta)




cdef double erlang_rv(double k , double theta):
    """
    Generate an Erlang random variable with k independent exponential steps each with mean theta. The mean is k*theta
     and the variance is k*theta^2
    :param k: (double) number of exponentially distributed steps. will be rounded to the nearest integer.
    :param theta: (double) mean of each step.
    :return: (double) the erlang random variable.
    """
    cdef double answer = 0.0
    cdef unsigned i = 0
    cdef unsigned num_iterations = int(k+0.5)
    for i in range(num_iterations):
        answer += -1.0 * theta * log( uniform_rv() )
    return answer

def py_erlang_rv(double k, double theta):
    return erlang_rv(k,theta)





cdef int sample_discrete(int choices, double* data, double Lambda):
    """
    Sample from a discrete set of options according to some un-normalized probability weights.
    :param choices: (int) the number of possible choices. should be positive.
    :param data: (double *) pointer to an array containing an un-normalized weight associated with each choice. (must
            have length = choices)
    :param Lambda: (double) the sum of all the probability weights in data. Used as the normalization factor.
    :return: (int) A non-negative index randomly sampled according to the weights. 0 <= return < choices
    """
    cdef double q = uniform_rv()*Lambda
    cdef int i = 0
    cdef double p_sum = 0.0
    while p_sum < q and i < choices:
        p_sum += data[i]
        i += 1
    return i - 1

def py_sample_discrete(int choices, np.ndarray[np.double_t,ndim=1] data, double Lambda):
    return sample_discrete(choices, <double*> (data.data), Lambda)





cdef double array_sum(double* data, int length):
    """
    Sum an array of floating point numbers
    :param data: (double *) pointer to the array of numbers (must have len(data) >= length)
    :param length: (int) the length of the array.
    :return: (double) the sum of all the numbers
    """
    cdef double answer = 0.0
    cdef int i = 0
    for i in range(length):
        answer += data[i]
    return answer


def py_array_sum(np.ndarray[np.double_t,ndim=1] data, int length):
    return array_sum(<double*> data.data, length)



cdef unsigned binom_rnd(unsigned n, double p):
    """
    Generate a binomial random number selected from binom(n,p)
    :param n: (unsigned) the population size, n>= 1
    :param p: (double) 0 <= p <= 1, the probability of success.
    :return: (unsigned) binomial random number
    """
    cdef unsigned answer = 0
    cdef unsigned i = 0
    for i in range(n):
        if uniform_rv() < p:
            answer += 1

    return answer

def py_binom_rnd(unsigned n, double p):
    return binom_rnd(n,p)


cdef unsigned binom_rnd_f(double N, double p):
    """
    Generate a binomial random number selected from binom(N,p), N >= 1
    :param N: (double) will be rounded to nearest integer
    :param p: (double) 0 <= p <= 1, the probability of success.
    :return: (unsigned) binomial random numbe
    """
    cdef unsigned answer = 0
    cdef unsigned n = int(N+0.5)
    cdef unsigned i = 0
    for i in range(n):
        if uniform_rv() < p:
            answer += 1

    return answer

def py_binom_rnd_f(double N, double p):
    return  binom_rnd_f( N,  p)



cdef unsigned approx_binom_rnd(unsigned n, double p):
    """
    Generate an approximately binomial random variable using the normal approximation.
    :param n: (unsigned) The population size n >= 1, but should be at least 20 for good approximation.
    :param p: (double) 0 <= p <= 1. p closer to 0.5 gives a better approximation
    :return: (unsigned) random number sample.
    """
    return int( normal_rv( n*p , sqrt(n*p*(1-p))  ) + 0.5 )

def py_approx_binom_rnd(unsigned n, double p):
    return approx_binom_rnd(n,p)




cdef unsigned approx_binom_rnd_f(double n, double p):
    """
    Generate an approximately binomial random variable using the normal approximation.
    :param n: (double) The population size n >= 1, but should be at least 20 for good approximation.
    :param p: (double) 0 <= p <= 1. p closer to 0.5 gives a better approximation
    :return: (unsigned) random number sample.
    """
    return int( normal_rv( n*p , sqrt(n*p*(1-p))  ) + 0.5 )

def py_approx_binom_rnd_f(double n, double p):
    return approx_binom_rnd_f(n,p)
