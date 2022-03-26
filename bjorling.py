from sympy import integrate,diff,I,re,symbols,exp,rootof
from sympy.abc import u,v,t,s

a,b = symbols('a,b', real=True)

def bjorling(x,y,t=t,u=u,v=v,polar=False):
    """ Bjorling's formula for minimal surfaces that
    contain given curve (x(t),y(t),0) as geodesic
    x,y = curve on xy-plane (functions of t)
    t = parameter on curve
    u,v = parameters on surface
    if polar==False, u+iv is substituted for t
    else, exp(u+iv) is substituted for t
    return x,y,z as functions of (u,v)
    reference: A. Gray
      "Modern Differential Geometry of Curves and Surfaces"
        section 22.6
    """
    l = diff(x,t)**2 + diff(y,t)**2
    l = l.subs(t,2*t).simplify()
    l = rootof(s**2 - l, s, 1).subs(t,t/2)
    p = (x, y, I*integrate(l,t))
    if polar:
        p = (x.subs(t, exp(t)) for x in p)
    p = (re(x.subs(t, a+I*b)) for x in p)
    p = (x.subs(((a,u),(b,v))) for x in p)
    return (x.simplify() for x in p)
