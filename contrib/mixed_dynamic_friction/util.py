from getfem import *
from numpy import *
from scipy import *
from pylab import * # see http://matplotlib.sourceforge.net/



with_graphics=True
try:
    import getfem_tvtk
except:
    print "\n** Could NOT import getfem_tvtk -- graphical output disabled **\n"
    import time
    time.sleep(2)
    with_graphics=False

print 'Some tests with python',


# m=Mesh('load', '../../tests/meshes/disc_P2_h4.mesh')
m=Mesh('import', 'gmsh', '/media/disk/rectangularQ3.msh')


if with_graphics:
    fig = getfem_tvtk.Figure()
    fig.show_mesh(m, faces=0, edges=1)
    print "Press Q to continue.."
    fig.set_colormap('tripod')
    fig.loop()



# convergence graphics
params = {'backend': 'ps',
          'axes.labelsize': 20,
          'text.fontsize': 20,
          'legend.fontsize': 20,
          'xtick.labelsize': 20,
          'ytick.labelsize': 20,
          'text.usetex': True,
          'figure.figsize': (9,6)}
rcParams.update(params)


# elastodyn case (P2)
h = array([0.5, 1., 2., 4.]);
el2 = array([0.416, 0.891, 3.236, 8.909]); al2=stats.linregress(log(h),log(el2));
eh1 = array([0.619, 1.095, 3.388, 9.098]); ah1=stats.linregress(log(h),log(eh1));
# ms = taille des "points", mfc = couleurs des "points"
# alpha = transparence des "points", lw = linewidth 
# fig = figure(figsize=(9,6))
grid(True)


# ax = fig.add_subplot(111)
loglog(h, el2, 'o-', ms=15, lw=2, alpha=0.9, mfc='orange')
loglog(h, eh1, 'd-', ms=15, lw=2, alpha=0.9, mfc='red')
loglog(h, exp(al2[0]*log(h)+al2[1]), 'k:')
loglog(h, exp(ah1[0]*log(h)+ah1[1]), 'k:')
xlim((0.4, 10.0))
ylim((0.2, 17.0))
# lines marker : [ ‘+’ | ‘*’ | ‘,’ | ‘.’ | ‘1’ | ‘2’ | ‘3’ | ‘4’ | ‘<’ | ‘>’ | ‘D’ | ‘H’ | ‘^’ | ‘_’ | ‘d’ | ‘h’ | ‘o’ | ‘p’ | ‘s’ | ‘v’ | ‘x’ | ‘|’ | TICKUP | TICKDOWN | TICKLEFT | TICKRIGHT | ‘None’ | ‘ ‘ | ‘’ ]
# line style [ ‘-‘ | ‘_’ | ‘-.’ | ‘:’ | ‘None’ | ‘ ‘ | ‘’ ] 

# title('the title', fontsize=20)
xlabel(r"mesh size $h$", fontsize=20)
ylabel(r'$H_1$ and $L_2$ error', fontsize=20)
legend((r'$L_2$ error (rate %f)' % al2[0], r'$H_1$ error (rate %f)' % ah1[0]), 'upper left', shadow=True)
show()

# scalar case (P2)
h = array([0.0125, 0.025, 0.05, 0.1, 0.2]);
el2 = array([0.001547, 0.003366, 0.005068, 0.007915, 0.007619]); al2=stats.linregress(log(h),log(el2));
eh1 = array([0.04558,  0.06532, 0.06480, 0.07417, 0.0751705]); ah1=stats.linregress(log(h),log(eh1));
loglog(h, el2, 'o-', ms=15, lw=2, alpha=0.9, mfc='orange')
loglog(h, eh1, 'd-', ms=15, lw=2, alpha=0.9, mfc='red')
loglog(h, exp(al2[0]*log(h)+al2[1]), 'k:')
loglog(h, exp(ah1[0]*log(h)+ah1[1]), 'k:')
xlim((0.4, 10.0))
ylim((0.2, 17.0))
xlabel(r"mesh size $h$", fontsize=20)
ylabel(r'$H_1$ and $L_2$ error', fontsize=20)
legend((r'$L_2$ error (rate %f)' % al2[0], r'$H_1$ error (rate %f)' % ah1[0]), 'upper left', shadow=True)
show()



if 0:

    x = array([0,1,2,3,4]); y = array([1,5,3,9,5])
    a=stats.linregress(x,y)
    z = [amin(x), amax(x)]
    
    plot(x, y, z, a[0]*z+a[1])
    
    
    title(r'$\ddot{o}\acute{e}\grave{e}\hat{O}\breve{i}\bar{A}\int_a^b\tilde{n}\vec{q}$', fontsize=20)
    
    xlabel(r"""$\"o\ddot o \'e\`e\~n\.x\^y$""", fontsize=20)
    
    ylabel(r'toto$\"u$', fontsize=20)


    show()
