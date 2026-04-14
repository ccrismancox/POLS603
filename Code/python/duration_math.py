#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 14:25:46 2026

@author: cox
"""
from sympy import * 
init_printing()

x = symbols("x", cls=IndexedBase, real=True)
y,s = symbols("y,s", cls=IndexedBase, real=True, positive=True)

a =symbols("alpha", real=True, positive=True)
b = symbols("beta", real=True)
i,N=symbols("i, N", real=True, positive=True, integer=True)
A=symbols("a", real=True, integer=True)


## exponential 
### set A=1 for PH and -1 for AFT
l_exp =  Sum(s[i]*A*x[i]*b-exp(A*x[i]*b)*y[i], (i,1,N))
l_exp.diff(b)


## weibull aft
gamma = exp(x[i]*b)
lam_a = gamma**(-a)
lhaz_aft = log(a)+log(lam_a)+(a-1)*log(y[i])
lS_aft = -lam_a*(y[i]**a)

l_wei_aft =  Sum(s[i]*lhaz+ lS, (i,1,N))
l_wei_aft.diff(b)
l_wei_aft.subs(a, exp(a)).diff(a).subs(exp(a),a)


## weibull ph
lam= exp(x[i]*b)
lhaz_ph = log(a)+log(lam)+(a-1)*log(y[i])
lS_ph= -lam*(y[i]**a)

l_wei_ph =  Sum(s[i]*lhaz_ph+ lS_ph, (i,1,N))
l_wei_ph.diff(b)
l_wei_ph.subs(a, exp(a)).diff(a).subs(exp(a),a)


## expected values 

# exp
x = symbols("x", real=True)
y,s = symbols("y,s", real=True, positive=True)
fe = simplify(exp(A*x*b)*exp(-exp(A*x*b)*y))
fe_ph=fe.subs(A,1)
integrate(fe_ph*y, (y,0,oo)) ## hazard
fe_aft=fe.subs(A,-1)
Ee=integrate(fe_aft*y, (y,0,oo)) ## Duration

## marginal
Ee.diff(x)




# weibull aft
gamma = exp(x*b)
lam_a = gamma**(-a)
haz_aft = a*lam_a*y**(a-1)
S_aft = exp(-lam_a*(y**a))
fw_aft = simplify(haz_aft*S_aft)
integrate(fw_aft*y, (y,0,oo)).simplify()
1-S_aft


# weibull ph
lam= exp(x*b)
haz_ph = a*lam*y**(a-1)
S_ph = exp(-lam*(y**a))
fw_ph = simplify(haz_ph*S_ph)
integrate(fw_ph*y, (y,0,oo)).simplify()
1-S_ph


haz_ph.diff(x)
