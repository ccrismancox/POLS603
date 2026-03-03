#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 21:53:46 2026

@author: cox
"""
from sympy import * 
init_printing()

y,x = symbols("y,x", cls=IndexedBase, real=True)
b= symbols("beta", real=True)
i,N=symbols("i, N", real=True, positive=True, integer=True)

    
L = Sum(y[i]*x[i]*b -exp(x[i]*b), (i,1,N))

s = L.diff(b).simplify()
H = s.diff(b).simplify()






### marginal
### here we're looking at ame(x1) with x2 being all other x's
b1, b2 = symbols("beta_1, beta_2", real=True)
x1, x2,l = symbols("x_1, x_2,lambda", cls=IndexedBase, real=True)

ame = 1/N * Sum(b1 * exp(x1[i]*b1+ x2[i]*b2), (i,1,N))
Db1 = ame.diff(b1)
Db1 = Db1.replace((x1[i]*b1+ x2[i]*b2), l[i])

Db = ame.diff(b2)
Db = Db.replace((x1[i]*b1+ x2[i]*b2), l[i])


Db1.expand()
Db



## bias?
## Here we consider a different (but identical) version of the score function
## for a model with just the constant
l, lhat = symbols("ell, \\hat{\\lambda}")
score0 = Sum(y[i]-exp(b), (i,1,N)).expand().doit()
bhat = solve(score0, b)[0].replace(Sum(y[i],(i,1,N)),N*lhat).simplify()
bhat

## derives for a 2nd order Taylor
D1 = bhat.diff(lhat).simplify().subs(lhat,l)
D2 = D1.diff(l).simplify()

tay = bhat.subs(lhat,l) + D1 *(lhat-l) +   D2*(lhat-l)**2 / 2
Etay = tay.subs([((lhat-l)**2, l/N),(lhat-l,0) ])
Etay
