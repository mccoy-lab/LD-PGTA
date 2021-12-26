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
        self.half_of_bits = self.bits >> 1
        self.mask = (1 << self.bits_in_result) - 1
        self.M = tuple(self.constants().items()) 
        
    def adjust_bits(self, bits):
        """ Finds the smallest power of two that is larger than `bits`,
            denoted as `i`. Then, finds the smallest power of two, denoted as
            `j`, for which 2 to the power of `j` is larger than `i`. """
        i = 8; j = 8;
        while(i < bits): i <<= 1;
        while((1<<j) < i): j <<= 1;
        return i, j
    
    def constants(self):
        """ Constants used in the functions `__call__`. """
        i = 1
        M = {}
        while(i!=self.bits_in_result):
            M[i] = int(('0' * i + '1' * i ) * ((self.bits >> 1) // i), 2)
            i <<= 1
        return M
    
    def __call__(self,x):
        """ Counts non-zero bits in positive integer. """
        bt = x.bit_length()
        if x<0: 
            return -1
        if x==0:
            return 0
        if bt > self.bits: ### or (self.bits > 64 and x.bit_length() <= self.half_of_bits):
            self.__init__(x.bit_length())
        for i,m in self.M:
            x = (x & m) + ((x >> i) & m)
        while(i<bt):
            i <<= 1
            x += x >> i
        return x & self.mask; 