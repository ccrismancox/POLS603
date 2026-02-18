#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 18 08:11:51 2026

@author: cox
"""
## Load the package
from sympy import * 

##optional but makes the printed math nicer
init_printing(forecolor="White") 

## define parameters. Do 2 endogenous variables to match example
g1, g2,b,t, r1, r2 =symbols("gamma_1, gamma_2 beta, tau, rho_1, rho_2", real=True)
s = symbols("sigma^2", positive=True)
i, N= symbols("i, N", positive=True, integer=True)

## define data
y = IndexedBase('y', real=True, positive=True)
x = IndexedBase('x', real=True)
z = IndexedBase('z', real=True)
d = IndexedBase('d', real=True)

## some shortcuts 
v1 = d[i]-z[i]*g1
v2 = d[i]-z[i]*g2
xb = x[i]*b + d[i]*t+ v1*r1+v2*r2

## likelihoods 
L1a = -N*log(2*pi)/2 -  N*log(s)/2 - Sum((d[i]-z[i]*g1)**2, (i,1,N))/ (2*s)
L1b = -N*log(2*pi)/2 -  N*log(s)/2 - Sum((d[i]-z[i]*g2)**2, (i,1,N))/ (2*s)
L2 = Sum(y[i]*xb -exp(xb),(i,1,N) )


## For H2, note that we are fixing the new inputs so the derive wrt to beta
## will give us the pattern we need for coding the full H2
H2 = L2.diff(b).diff(b)
## -X' (Iexp(XB)) X like the notes

## more important for us is C

## C_OPG needs the Jacobian wrt to gamma
DL2_Dg = L2.diff(g1)
## Z' rho * (lambda-y)


## C_H uses cross derivatives
C_bg = -L2.diff(b).diff(g1)
C_tg = -L2.diff(t).diff(g1)
C_r1g = -L2.diff(r1).diff(g1)
C_r2g = -L2.diff(r2).diff(g1)


C_bg
## rho * (X' (Iexp(XB)) Z)

C_r1g
## rho * (v' *lambda) +(lambda - y)'Z

C_r2g
## rho * (v' *lambda)
