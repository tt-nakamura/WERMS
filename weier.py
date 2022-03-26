from sympy import integrate,I,re,symbols,exp
from sympy.abc import u,v,w

a,b = symbols('a,b', real=True)

def weier(f,g,w=w,u=u,v=v,polar=False):
    """ Weierstrass-Enneper representation of minimal surfaces
    f,g = analytic functions to generate surface
    w = complex independent variable of f,g
    u,v = parameters on surface
    if polar==False, u+iv is substituted for w
    else, exp(u+iv) is substituted for w
    return x,y,z as functions of (u,v)
    reference: A. Gray
      "Modern Differential Geometry of Curves and Surfaces"
        3rd edition, section 22.5
    """
    p = (integrate(f*(1-g**2), w),
         integrate(f*(1+g**2), w)*I,
         integrate(f*g, w)*2)
    if polar:
        p = (x.subs(w, exp(w)) for x in p)
    p = (re(x.subs(w, a+I*b)) for x in p)
    p = (x.subs(((a,u),(b,v))) for x in p)
    return (x.simplify() for x in p)
