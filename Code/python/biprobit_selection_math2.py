#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  9 14:14:00 2025

@author: cox
"""
x, z, d,g,b, r, y, t,c, w= symbols("x, z, d, gamma, beta, rho, y, tau, psi, w", real=True)
s = symbols("sigma", positive=True)
Phi = Function("Phi", real=True)
phi = Function("phi")
Phi2 = Function("Phi_2")
phi2 = Function("phi_2")

class Phi(Function):
    def fdiff(self, argindex=1):

        # argindex indexes the args, starting at 1

        return phi(self.args[0])
    
class phi(Function):
    def fdiff(self, argindex=1):

        # argindex indexes the args, starting at 1

        return phi(self.args[0])*(-self.args[0])

class Phi2(Function):

    def fdiff(self, argindex):
        if argindex == 1:
            return phi(self.args[0])*Phi( (self.args[1] - self.args[2]*self.args[0])/sqrt(1-self.args[2]**2))
        elif argindex == 2:
            return phi(self.args[1])*Phi( (self.args[0] - self.args[2]*self.args[1])/sqrt(1-self.args[2]**2))
        else:
            return phi2(self.args[0], self.args[1], self.args[2])
    
    
class phi2(Function):
    def fdiff(self, argindex):
        if argindex == 3:
          return phi2(self.args[0], self.args[1], self.args[2]) *( (-self.args[2]*(-2*self.args[2]*self.args[0]*self.args[1]+self.args[1]**2+self.args[0]**2))/((1-self.args[2]**2)**2) + self.args[2]/(1-self.args[2]**2) + self.args[1]*self.args[0]/(1-self.args[2]**2))
    
        
    
E  = Phi2(z*g +w*c, x*b +w*t, r)/ Phi(z*g +w*c) 
ame1 = E.diff(x)
ame1

m1, m2 = symbols("mu_1,mu_2", )

ame2 = E.diff(w).subs(z*g+w*c, m2).subs(x*b+w*t, m1)
ame2
