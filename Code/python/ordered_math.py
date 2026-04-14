#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 14:08:06 2026

@author: cox
"""
from sympy import * 
init_printing()

x = symbols("x", cls=IndexedBase, real=True)
y0, y1, y2, y3, y4 = symbols("y_0, y_1, y_2, y_3, y_4", cls=IndexedBase, real=True, positive=True, integer=True)

b, t0= symbols("beta, tau_0", real=True)
d1, d2, d3 = symbols("delta_1, delta_2, delta_3", real=True)

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
        
    
xb = x[i]*b
t1 = t0+exp(d1)
t2 = t0 + exp(d1)+exp(d2)
t3 = t0 + exp(d1)+exp(d2)+exp(d3)
L = Sum(y0[i]*log(G(t0-xb)-G(-oo-xb)) + y1[i]*log(G(t1-xb)-G(t0-xb))+ 
        y2[i]*log(G(t2-xb)-G(t1-xb))+ y3[i]*log(G(t3-xb)-G(t2-xb)), (i,1,N))


tau1, tau2, tau3 = symbols("tau1, tau2, tau3", real=True)
Db=L.diff(b)
Db=Db.replace(G(-oo),0).replace(G(oo),1)
Db=Db.simplify()

Db=Db.replace(G(-xb+t1), G(tau1-xb)).replace(G(-xb+t2), G(tau2-xb)).replace(G(-xb+t3), G(tau3-xb))
Db


Dt=L.diff(t0)
Dt=Dt.replace(G(-oo),0).replace(G(oo),1)
Dt=Dt.simplify()

Dt=Dt.replace(G(-xb+t1), G(tau1-xb)).replace(G(-xb+t2), G(tau2-xb)).replace(G(-xb+t3), G(tau3-xb))
Dt

Dd1=L.diff(d1)
Dd1=Dd1.replace(G(-oo),0).replace(G(oo),1)

Dd1=Dd1.replace(G(-xb+t1), G(tau1-xb)).replace(G(-xb+t2), G(tau2-xb)).replace(G(-xb+t3), G(tau3-xb))
Dd1

Dd2=L.diff(d2)
Dd2=Dd2.replace(G(-oo),0).replace(G(oo),1)
Dd2=Dd2.replace(G(-xb+t1), G(tau1-xb)).replace(G(-xb+t2), G(tau2-xb)).replace(G(-xb+t3), G(tau3-xb))
Dd2








xb = x[i]*b
t1 = t0+exp(d1)
t2 = t0 + exp(d1)+exp(d2)
t3 = t0 + exp(d1)+exp(d2)+exp(d3)
L = Sum(y0[i]*log(Phi(t0-xb)-Phi(-oo-xb)) + y1[i]*log(Phi(t1-xb)-Phi(t0-xb))+ 
        y2[i]*log(Phi(t2-xb)-Phi(t1-xb))+ y3[i]*log(Phi(t3-xb)-Phi(t2-xb)), (i,1,N))


tau1, tau2, tau3 = symbols("tau1, tau2, tau3", real=True)
Db=L.diff(b)
Db=Db.replace(Phi(-oo),0).replace(Phi(oo),1)
Db=Db.simplify()

Db=Db.replace(Phi(-xb+t1), Phi(tau1-xb)).replace(Phi(-xb+t2), Phi(tau2-xb)).replace(Phi(-xb+t3), Phi(tau3-xb))
Db=Db.replace(phi(-xb+t1), phi(tau1-xb)).replace(phi(-xb+t2), phi(tau2-xb)).replace(phi(-xb+t3), phi(tau3-xb))
Db

Dt=L.diff(t0)
Dt=Dt.replace(Phi(-oo),0).replace(Phi(oo),1)
Dt=Dt.simplify()

Dt=Dt.replace(Phi(-xb+t1), Phi(tau1-xb)).replace(Phi(-xb+t2), Phi(tau2-xb)).replace(Phi(-xb+t3), Phi(tau3-xb))
Dt=Dt.replace(phi(-xb+t1), phi(tau1-xb)).replace(phi(-xb+t2), phi(tau2-xb)).replace(phi(-xb+t3), phi(tau3-xb))
Dt

Dd1=L.diff(d1)
Dd1=Dd1.replace(Phi(-oo),0).replace(Phi(oo),1)

Dd1=Dd1.replace(Phi(-xb+t1), Phi(tau1-xb)).replace(Phi(-xb+t2), Phi(tau2-xb)).replace(Phi(-xb+t3), Phi(tau3-xb))
Dd1

Dd2=L.diff(d2)
Dd2=Dd2.replace(Phi(-oo),0).replace(Phi(oo),1)
Dd2=Dd2.replace(Phi(-xb+t1), Phi(tau1-xb)).replace(Phi(-xb+t2), Phi(tau2-xb)).replace(Phi(-xb+t3), Phi(tau3-xb))
Dd2
