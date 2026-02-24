#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 21:53:46 2026

@author: cox
"""
from sympy import * 
init_printing()

y,x, z = symbols("y,x,z", cls=IndexedBase, real=True)
b= symbols("beta", real=True)
i,N=symbols("i, N", real=True, positive=True, integer=True)
g = Function("g")

class Phi(Function):
    ## define what the derivative of Phi is
    def fdiff(self, argindex=1): 
        # argindex indexes the args, starting at 1
        return phi(self.args[0])
    
class phi(Function):
    ## define what the derivative of phi is
    def fdiff(self, argindex=1):
        # argindex indexes the args, starting at 1
        return phi(self.args[0])*(-self.args[0])
    
class G(Function):
    ## We'll call is G instead of Lambda, because
    ## lambda means something special in python 
    ## and can mess things up
    ## define the derivative
    def fdiff(self, argindex=1):
        # argindex indexes the args, starting at 1

        return G(self.args[0])*(1-G(self.args[0]))
    
Lp = Sum(log(Phi(z[i]*b)), (i,1,N))

sp = Lp.diff(b)
sp
hp = sp.diff(b).simplify()
hp


Ll= Sum(log(G(z[i]*b)), (i,1,N))

sl = Ll.diff(b)
sl
hl = sl.diff(b).simplify()
hl
