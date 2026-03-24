#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 17:02:58 2026

@author: cox
"""
from sympy import * 
init_printing()

y,x,z,s = symbols("y,x,z,s", cls=IndexedBase, real=True)
b,b2, g= symbols("beta, beta_2,gamma", real=True)
u,t = symbols("eta, tau", real=True)
i,N=symbols("i, N", real=True, positive=True, integer=True)
r, s2= symbols("rho,sigma^2")

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
    

L = s*log(exp(-(1/2) * (y[i]-x[i]*b)**2 / exp(t))/sqrt(exp(t)*2*pi)*Phi((tanh(u)*(y[i]-x[i]*b)/sqrt(exp(t)) + z[i]*g)/sqrt(1-tanh(u)**2))) + (1-s)*log(Phi(-z[i]*g))
L = expand_log(L)
L.replace((tanh(u)*(y[i]-x[i]*b)/sqrt(exp(t)) + z[i]*g)/sqrt(1-tanh(u)**2), mu)

Db = L.diff(b)
m = symbols("mu", real=True)
Db = Db.replace((tanh(u)*(y[i]-x[i]*b)/sqrt(exp(t)) + z[i]*g)/sqrt(1-tanh(u)**2), mu)
Db = Db.replace(exp(-t), 1/s2).replace(exp(-t/2), 1/sqrt(s2)).replace(tanh(u), r)
Db.collect(x[i])


Dt = L.diff(t)
Dt = Dt.replace((tanh(u)*(y[i]-x[i]*b)/sqrt(exp(t)) + z[i]*g)/sqrt(1-tanh(u)**2), mu)
Dt = Dt.replace(exp(-t), 1/s2).replace(exp(-t/2), 1/sqrt(s2)).replace(tanh(u), r)
Dt


De = L.diff(u)
De = De.replace((tanh(u)*(y[i]-x[i]*b)/sqrt(exp(t)) + z[i]*g)/sqrt(1-tanh(u)**2), mu)
De = De.replace(exp(-t), 1/s2).replace(exp(-t/2), 1/sqrt(s2)).replace(tanh(u), r)
De.simplify()

Dg = L.diff(g)
Dg = Dg.replace((tanh(u)*(y[i]-x[i]*b)/sqrt(exp(t)) + z[i]*g)/sqrt(1-tanh(u)**2), mu)
Dg = Dg.replace(exp(-t), 1/s2).replace(exp(-t/2), 1/sqrt(s2)).replace(tanh(u), r)
Dg.collect(z[i])



