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

class g(Function):
    ## We'll call is G instead of Lambda, because
    ## lambda means something special in python 
    ## and can mess things up
    ## define the derivative
    def fdiff(self, argindex=1):
        # argindex indexes the args, starting at 1

        return g(self.args[0])*(1-2*G(self.args[0]))
        
    
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






### marginal
### here we're looking at ame(x1) with x2 being all other x's
b1, b2 = symbols("beta_1, beta_2", real=True)
x1, x2,m = symbols("x_1, x_2,mu", cls=IndexedBase, real=True)

ame_p = 1/N * Sum(b1 * phi(x1[i]*b1+ x2[i]*b2), (i,1,N))
Db1 = ame_p.diff(b1)
Db1 = Db1.replace((x1[i]*b1+ x2[i]*b2), m[i])
Db1 = Db1.replace(-(x1[i]*b1+ x2[i]*b2), -m[i])

Db = ame_p.diff(b2)
Db = Db.replace((x1[i]*b1+ x2[i]*b2), m[i])
Db = Db.replace(-(x1[i]*b1+ x2[i]*b2), -m[i])


Db1.expand()
Db



ame_l = 1/N * Sum(b1 * g(x1[i]*b1+ x2[i]*b2), (i,1,N))
Db1 = ame_l.diff(b1)
Db1 = Db1.replace((x1[i]*b1+ x2[i]*b2), m[i])

Db = ame_l.diff(b2)
Db = Db.replace((x1[i]*b1+ x2[i]*b2), m[i])


Db1
Db
