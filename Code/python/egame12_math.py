#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 14:19:03 2026

@author: cox
"""
from sympy import * 
init_printing()

x1, x2, x3,z = symbols("x_1, x_2, x_3, z", cls=IndexedBase, real=True)
a1, a2 = symbols("a_1, a_2", cls=IndexedBase, real=True, positive=True, integer=True)

b1, b2, b3, gamma =symbols("beta_1, beta_2, beta_3, gamma", real=True)

i,N =symbols("i, N", real=True, positive=True, integer=True)
g = Function("g")

class G(Function):
    ## We'll call is G instead of Lambda, because
    ## lambda means something special in python 
    ## and can mess things up
    ## define the derivative
    def fdiff(self, argindex=1):
        # argindex indexes the args, starting at 1

        return G(self.args[0])*(1-G(self.args[0]))

class g(Function):
    ## We'll call is G instead of Lambda, because
    ## lambda means something special in python 
    ## and can mess things up
    ## define the derivative
    def fdiff(self, argindex=1):
        # argindex indexes the args, starting at 1

        return g(self.args[0])*(1-2*G(self.args[0]))
        
    
u2 = z[i]*gamma
p2 = G(u2)
u1 = -x1[i] *b1 + p2*x3[i]*b3 +(1-p2)*x2[i]*b2
p1 = G(u1)
sq = (1-a1[i])
bd = a1[i]*(1-a2[i])
sf = a1[i]*a2[i]

L = Sum( sq*log(1-p1) + bd*(log(p1)+log(1-p2)) + sf*(log(p1)+log(p2)), (i,1,N))

P1, P2, U1,d1 =symbols("P1,P2, U1, d1")

L.diff(b1).replace(p1, P1).simplify()
L.diff(b2).replace(p1, P1).replace(p2,P2).simplify()
L.diff(b3).replace(p1, P1).replace(p2,P2).simplify()
L.diff(gamma).replace(p1, P1).replace(u1.diff(gamma), d1).replace(p2, P2).collect(z[i])
