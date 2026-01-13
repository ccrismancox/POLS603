#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 24 12:30:39 2025

@author: cox
"""
## Load the package
from sympy import * 

##optional but makes the printed math nicer
init_printing(forecolor="White") 

## We first need to declare what our variables are
m = symbols("mu", real=True)
s = symbols("sigma^2", positive=True)
i, N = symbols("i, N", positive=True, integer=True)


## This will be used to sum over
y = IndexedBase('y', real=True)

## Now we can define the likelihood
L = -N*log(2*pi)/2 -  N*log(s)/2 - Sum((y[i]-m)**2, (i,1,N))/ (2*s)


L

## First derivatives
FOC_mu = L.diff(m)
FOC_s = L.diff(s)

## Take a look
FOC_mu
FOC_s

## different ways to simplify and solve for mu
solve(FOC_mu, m)
simplify(FOC_mu)
simplify(FOC_mu).doit()
simplify(FOC_mu).expand()
simplify(FOC_mu).expand().doit()
solve(simplify(FOC_mu).expand().doit(),m)



## Finding the MLEs
FOC_mu = -N*m/s + Sum(y[i], (i,1,N))/s
mu_hat = solve(FOC_mu, m)
s2_hat = solve(FOC_s, s)

mu_hat, s2_hat




## Second order conditions
D2_mu = FOC_mu.diff(m)
D2_mu

D2_s = FOC_s.diff(s)
e = symbols("e", real=True)
D2_s = D2_s.subs([(s,e*e/N),(m,mu_hat[0])])
D2_s = D2_s.subs([(Sum((y[i]-mu_hat[0])**2, (i,1,N)), e*e)])
D2_s


