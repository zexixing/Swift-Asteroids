import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
from scipy.special import erf
from scipy import interpolate

def gauss(x, A, mu, sigma, normalize=True):
    mi = -np.power((x-mu)/sigma,2)/2
    if normalize:
        norm = 1./(sigma*np.sqrt(2*np.pi))
    else:
        norm = 1.
    return A*np.power(np.e,mi)*norm

def Fun1(p,y,x):
   '''compute the residuals for gaussian fit in get_components '''
   a0, x0, sig0 = p
   return y - gauss(x,a0,x0,sig0)

def gauss_motion(x, A, mu, sigma, motion, normalize=True):
    t1, t2 = -motion/2, motion/2
    m1, m2 = (t1-(x-mu))/(np.sqrt(2)*sigma), (t2-(x-mu))/(np.sqrt(2)*sigma)
    n1, n2 = m1*m1, m2*m2
    fifth = -(np.exp(-n2)-np.exp(-n1))
    sixth = np.sqrt(np.pi/2)*(x-mu)/sigma*(erf(m2)-erf(m1))
    forth = fifth + sixth
    third = np.exp(np.power((x-mu)/sigma,2)/2)*2*np.power(sigma,2)*forth
    secnd = -1/(2*np.power(sigma,2))*third
    def first_f(t):
        return np.exp(-np.power(t/sigma,2)/2+t*(x-mu)/np.power(sigma,2))
    first = first_f(t2)-first_f(t1)
    zeroth = np.power(sigma,2)/(x-mu)*(first - secnd)
    if normalize == True:
        norm = 1./(sigma*np.sqrt(2*np.pi))
    else:
        norm = 1.
    result = A*norm/motion*np.exp(-np.power((x-mu)/sigma,2)/2)*zeroth
    return result

def Fun2(p,y,x):
   '''compute the residuals for gaussian fit in get_components '''
   a0, x0, sig0, motion0 = p
   return y - gauss_motion(x,a0,x0,sig0,motion0)

def frac(sigma, motion):
    A = 50
    mu = 100
    x = np.arange(0.5,200.5)
    y = gauss_motion(x, A, mu, sigma, motion, normalize=True)
    f = interpolate.interp1d(x,y,fill_value="extrapolate")
    x1 = mu-(sigma+motion/4)*2.5
    x2 = mu+(sigma+motion/4)*2.5
    step = 0.1
    x_aster = np.arange(x1,x2+step,step)
    x_aster = (x_aster[:-1]+x_aster[1:])/2
    y_aster = f(x_aster)
    total = np.sum(y)*(x[1]-x[0])
    aster = np.sum(y_aster)*step
    return aster/total
print(frac(3.2,0.001))

def GaussianHalfIntegralFraction(x):
   ''' 
   Computes the normalised integrated gaussian from 
   -x to +x. For x=inf, the result equals 1.
   
   x is in units of sigma
   
   Abramowitz & Stegun, par. 26.2 
   '''
   import numpy as np    
   from past.utils import old_div
   Z_x = old_div(np.exp( -x*x/2.), np.sqrt(2.* np.pi))
   p = .33267
   a1 =  .4361836
   a2 = -.1201676
   a3 =  .9372980
   t = old_div(1.,(1.+p*x))
   P_x = 1. - Z_x * t* (a1 + t*(a2 + a3*t) ) 
   A_x = 2.*P_x - 1
   return  A_x  

print(GaussianHalfIntegralFraction(2.5))

x = np.arange(0,200)
A = 50
sigma=3.2
'''
y = np.zeros(x.shape)
for mu in np.arange(95,105)+0.5:
    y += gauss(x, A, mu, sigma, normalize=True)
plt.plot(x,y,'r-',label='motion=10pix')
(A0,mu0,sigma0), ier = leastsq(Fun1, (50,100,3.2), args=(y,x) )
plt.plot(x,gauss(x, A0, mu0, sigma0),'r--',alpha=0.3)

y = np.zeros(x.shape)
for mu in np.arange(92,108)+0.5:
    y += gauss(x, A, mu, sigma, normalize=True)
plt.plot(x,y,'b-',label='motion=16pix')
(A0,mu0,sigma0), ier = leastsq(Fun1, (50,100,3.2), args=(y,x) )
plt.plot(x,gauss(x, A0, mu0, sigma0),'b--',alpha=0.3)
'''
y = np.zeros(x.shape)
for mu in np.arange(89,111)+0.5:
    y += gauss(x, A, mu, sigma, normalize=True)
plt.plot(x,y,'g-',label='superposed')
#(A0,mu0,sigma0), ier = leastsq(Fun1, (50,100,3.2), args=(y,x) )
#plt.plot(x,gauss(x, A0, mu0, sigma0),'g--',alpha=0.3)
y = gauss_motion(x, A*22, 100, sigma, 22, normalize=True)
plt.plot(x,y,'k--',alpha=0.5,label='integration')
plt.legend()
plt.title('motion=22pixel')
plt.show()
'''

mu = 100
x = np.delete(x,np.where(x==mu))
motion=20
A = 50
y = gauss_motion(x, A, mu, sigma, motion, normalize=True)
plt.plot(x,y,'k-',label='motion=20pix')
plt.vlines(mu-(sigma+motion/4)*2.5,0,np.max(y),color='k')
plt.vlines(mu+(sigma+motion/4)*2.5,0,np.max(y),color='k')
#(A0,mu0,sigma0,motion0), ier = leastsq(Fun2, (50,100,3.2,20), args=(y,x) )
#plt.plot(x,gauss_motion(x, A0, mu0, sigma0, motion0),'r--')
print((sigma+motion/4)*2.5)


plt.legend()
plt.xlim(50,150)
plt.show()
'''