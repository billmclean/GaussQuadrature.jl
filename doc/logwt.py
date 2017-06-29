from sympy import *

x = symbols('x')

def p(l, x):
    C = binomial(2*l,l)
    return (legendre(l,2*x-1) / C).expand()

def w(r, x):
    return (x**r) * log(1/x)

def nu(r, l):
    return integrate(p(l-1,x)*w(r,x), (x, 0, 1))

r = 0
q0 = 1

alpha1 = Rational(1,2) + nu(r,2) / nu(r,1)
q1 = x - alpha1

for l in range(1,6):
    sigma = integrate(p(l-1,x)*q1*w(r,x), (x,0,1))
    print("sigma_{},2 = ".format(l), sigma)

sigma12 = integrate(p(0,x)*q1, (x,0,1))
sigma22 = integrate(p(1,x)*q1, (x,0,1))
sigma32 = integrate(p(2,x)*q1, (x,0,1))
