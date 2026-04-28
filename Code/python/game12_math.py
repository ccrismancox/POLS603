#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 15:53:21 2026

@author: cox
"""
from sympy import * 
init_printing()

y,x1,x2,x3, z = symbols("y,x_1, x_2, x_3,z", cls=IndexedBase, real=True)
b1,b2,b3,g0= symbols("beta_1, beta_2, beta_3, gamma", real=True)
i,N=symbols("i, N", real=True, positive=True, integer=True)
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
        
    
pb = G(z[i]*g0)
Usq = x1[i]*b1
Ubd = x2[i]*b2
Usf = x3[i]*b3
EU = -Usq + Ubd*(1-pb)+Usf*pb
p1 = G(EU)
L2 = Sum(y[i]*log(p1)+(1-y[i])*log(1-p1), (i,1,N))


P1, P2 =symbols("p1,p2", real=True)

L.diff(b1).replace(p1, P1).simplify()
L.diff(b2).replace(p1, P1).replace(p2,P2).simplify()
L.diff(b3).replace(p1, P1).replace(p2,P2).simplify()
L.diff(gamma).replace(p1, P1).replace(p2, P2).collect(z[i])



## For two step
L2gb1 = L2.diff(g0).diff(b1)
L2gb1.replace(p1, P1).replace(pb,P2)


L2gb2 = L2.diff(g0).diff(b2)
L2gb2.replace(p1, P1).replace(pb,P2).replace(P2*(1-P2), d2)

L2gb3 = L2.diff(g0).diff(b3)
L2gb3.replace(p1, P1).replace(pb,P2)
