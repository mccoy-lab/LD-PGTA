#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 
POPCOUNTS

A module that provides an efficent function to count non-zero bits in positive integer. 

Daniel Ariad (daniel@ariad.org)
Dec 23, 2021
"""

class get_popcount():
    """ This is a generalization to n bits of popcount64b from 
        http://en.wikipedia.org/wiki/Hamming_weight . """
        
        
    def __init__(self, bits = 64):
        """ Initialization function """
        self.bits, self.bits_in_result = self.adjust_bits(bits)
        self.mask = (1 << self.bits_in_result) - 1
        self.M = tuple(self.constants().items()) 
        
    def adjust_bits(self, bits):
        """ Finds the smallest power of two that is greater than or equal to `bits`,
            denoted as `i`. Then, finds the smallest power of two, denoted as
            `j`, for which 2 to the power of `j` is greater than or equal to `i`. """
        i = 8; j = 8;
        while(i < bits): i <<= 1;
        while((1<<j) < i): j <<= 1;
        return i, j
    
    def constants(self):
        """ Constants used in the functions `__call__`. """
        M = {}
        i = 1
        j = self.bits
        while(i!=self.bits_in_result):
            j >>= 1
            M[i] = int(('0' * i + '1' * i ) * j, 2)
            i <<= 1
        return M
    
    def __call__(self,x):
        """ Counts non-zero bits in positive integer. """
        bt = x.bit_length()
        if x<0: 
            return -1
        if x==0:
            return 0
        if bt > self.bits:
            self.__init__(x.bit_length())
        for i,m in self.M:
            x = (x & m) + ((x >> i) & m)
        while(i<bt):
            i <<= 1
            x += x >> i
        return x & self.mask; 
    
if __name__ == "__main__":
    print('The module POPCOUNTS was invoked directly.')
    import random, sys
    
    popcount_naive = lambda x: bin(x).count('1')
    random.seed(a=2022, version=2)
    naive_result = [popcount_naive(random.getrandbits(1024)) for i in range(1000000)]
    
    popcount = get_popcount()
    random.seed(a=2022, version=2)
    result = [popcount(random.getrandbits(1024)) for i in range(1000000)]
    
    print('Algorithm validation: %s' % ('passed' if naive_result==result else 'failed'))
    sys.exit(0)
    
else:
    print("The module POPCOUNTS was imported.")