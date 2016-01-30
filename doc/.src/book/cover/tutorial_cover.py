"""Comic strip for illustrating Euler's method for ODEs."""

from pysketcher import *
import numpy as np
xkcd = True
#xkcd = False

xmin = 0
drawing_tool.set_coordinate_system(xmin=xmin, xmax=4,
                                   ymin=0, ymax=2.5,
                                   axis=True, xkcd=xkcd)
drawing_tool.set_linecolor('blue')
# Turn off tick marks at the top and the right axes
drawing_tool.mpl.tick_params(
    axis='x',
    top='off',
    labeltop='off',
    )
drawing_tool.mpl.tick_params(
    axis='y',
    right='off',
    labelright='off',
    )

def ForwardEuler(I, a, T, dt):
    u = [I]
    t = [0]
    while t[-1] <= T:
        u_new = u[-1] - a*dt*u[-1]
        u.append(u_new)
        t.append(t[-1] + dt)
    return np.array(u), np.array(t)

def make_fig(dt=0.5, heading=''):
    I = 2
    a = 0.5
    T_e = 3
    T_FE = 1
    t_fine = np.linspace(0, T_e, 101)
    u_e = I*np.exp(-a*t_fine)

    u, t = ForwardEuler(I, a, T_FE, dt)

    # y = slope*(x - x0) + y0
    # u_extrapolated = -a*u[-1]*(t - t[-1]) + u[-1]
    t_future = t[-1] + 1.5   # let the line be longer than one step
    line = Line((t[-1], u[-1]), (t_future, -a*u[-1]*(t_future - t[-1]) + u[-1]))

    circles = {
        i: Circle((t[i], u[i]), 0.05).set_linecolor('red').set_filled_curves('red')
        for i in range(1, len(u))}
    # Add next predicted point
    t_next = t[-1] + dt
    u_next = -a*u[-1]*(t_next - t[-1]) + u[-1]
    circles[0] = Circle((t_next, u_next), 0.05).\
                 set_linecolor('red').set_filled_curves('red')
    circles = Composition(circles)

    curves = Composition(dict(
        exact=Curve(t_fine, u_e).set_linestyle('dashed'),
        numerical=Curve(t, u),
        extrapolation=line.set_linecolor('red').set_linewidth(3)))

    tpos = 2.5
    tpos = 2.7
    pos = (tpos, 1)
    pos = (tpos, 0.8)
    text_exact = Text_wArrow("exact solution", pos, (tpos, I*np.exp(-a*tpos)),
                             alignment='left')

    pos = (1.7, 1.7)
    pos = (1.7, 1.5)
    text_predict = Text_wArrow("Here we know the slope:\n   u' = f(u)\nLet the solution continue\nlinearly along that slope!",
                               pos, (t[-1], u[-1]),
                               alignment='left')
    text_next = Text_wArrow("New predicted point",
                            (1, 0.25), (t_next, u_next),
                            alignment='left')

    fig = Composition(dict(curves=curves,
                           circles=circles,
                           exact=text_exact,
                           predict=text_predict,
                           next=text_next))
    if heading:
        fig['comment'] = Text(heading, (0.3, 2.05), alignment='left')

    return fig

#fig = make_fig(dt=0.5, heading="How to solve u' = -au\nwith programming:")
fig = make_fig(dt=0.5, heading="")
fig.draw()
drawing_tool.display()

# Springer demands 300 dpi for cover image
drawing_tool.savefig('tmp1.png', dpi=300)

import os, commands
# Get file size
failure, output = commands.getstatusoutput('identify tmp1.png')
width, height = output.split()[2].split('x')
# Remove background color
os.system('convert tmp1.png -transparent white tmp1b.png')
# ImageMagick colors: http://www.imagemagick.org/script/color.php
gradient = 'navy-snow'
gradient = 'yellow-white'
gradient = 'lightyellow-white'
gradient = 'orange-white'
gradient = 'lime-lightyellow'
gradient = 'RoyalBlue1-SlateGray1'
gradient = 'SkyBlue1-white'
gradient = 'SeaGreen3-white'  # 2
gradient = 'DarkSeaGreen2-white'  # 3
gradient = 'PaleGreen3-white'  # 1
gradient = 'yellow2-LemonChiffon' # 1
gradient = 'PaleGreen3-LightYellow'  # FINAL!
#gradient = 'RoyalBlue1-SlateGray1'

# Combine plot and graded background
filename = 'tmp0-%s.png' % gradient
os.system('convert -size %sx%s gradient:%s tmp1b.png -composite %s' % (width, height, gradient, filename))

# Generate a web page that mimics the front cover
outfile = open('tmp.html', 'w')
outfile.write(u"""
<body bgcolor="black">
  <center>
    <h1 style="font-size:300%%; color:white; font-family:helvetica;">Finite Difference Computing with<br> Exponential Decay Models</h1>
    <img src="%s" width=700>
    </center>
</body>
""" % filename)


#input()
