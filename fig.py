from sympy import pi,I,log,cos,sin,exp,arg
from sympy.abc import u,v,w,t
from sympy.plotting import plot3d_parametric_surface
from sympy.plotting import plot3d_parametric_line
from sympy.plotting import plot3d
from weier import weier
from bjorling import bjorling

# catenoid
x,y,z = weier(1/(2*w**2), w, polar=True)
print(x, ',', y, ',', z)
plot3d_parametric_surface(
    x,y,z,(u,-2,2),(v,-pi,pi),
    xlabel='x', ylabel='y', title='catenoid')

# helicoid
x,y,z = weier(I/(2*w**2), w, polar=True)
z = z.subs(arg(exp(I*v)), v)
print(x, ',', y, ',', z)
plot3d_parametric_surface(
    x,y,z,(u,-2,2),(v,-pi,pi),
    xlabel='x', ylabel='y', title='helicoid')

# Enneper
x,y,z = weier(1,w)
print(x,y,z)
plot3d_parametric_surface(
    x,y,z,(u,-1,1),(v,-1,1),
    xlabel='x', ylabel='y', title='Enneper')

# Scherk
x,y,z = weier(2/(1-w**4), w)
print(x, ',', y, ',', z)
plot3d(log(cos(v)/cos(u)), (u,-1.56,1.56),(v,-1.56,1.56),
       xlabel='x', ylabel='y', title='Scherk')

# Richmond
x,y,z = weier(w**2, 1/w**2)
print(x, ',', y, ',', z)
plot3d_parametric_surface(
    (x,y,z,(u,-0.5,0.5),(v,0.1,1)),
    (x,y,z,(u,-0.5,0.5),(v,-1,-0.1)),
    (x,y,z,(u,-0.5,-0.1),(v,-1,1)),
    (x,y,z,(u,0.1,0.5),(v,-1,1)),
    xlabel='x', ylabel='y', title='Richmond')

# bat
x,y,z = weier(I*w,w**3)
print(x, ',', y, ',', z)
plot3d_parametric_surface(
    x,y,z,(u,-1,1),(v,-1,1),
    xlabel='x', ylabel='y', title='bat')

# Catalan
p1 = plot3d_parametric_line(
    t-sin(t), 1-cos(t), t*1e-6, (t,-2*pi,2*pi), show=False)
x,y,z = bjorling(t-sin(t), 1-cos(t))
print(x, ',', y, ',', z)
p2 = plot3d_parametric_surface(
    x,y,z,(u,-2*pi,2*pi),(v,0,1),show=False,
    xlabel='x', ylabel='y', title='Catalan')
p2.extend(p1)
p2.show()

# asteroid
p1 = plot3d_parametric_line(
    cos(t)**3, sin(t)**3, t*1e-6, (t,0,2*pi), show=False)
x,y,z = bjorling(cos(t)**3, sin(t)**3)
print(x, ',', y, ',', z)
p2 = plot3d_parametric_surface(
    x,y,z,(u,0,2*pi),(v,-0.5,0),show=False,
    xlabel='x', ylabel='y', title='asteroid')
p2.extend(p1)
p2.show()
