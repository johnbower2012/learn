import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt   # For plotting graphs.
import numpy as np
import subprocess                 # For issuing commands to the OS.
import os
import sys                        # For determining the Python version.
import glob

data = np.loadtxt("../data/wave.dat")
columns = np.size(data,1)

style = {0: '--',
         1: '-',
         2: '-.',
         3: ':'}
color = {0: 'k',
         1: 'r',
         2: 'g',
         3: 'b'}
n = 0
print('Working...')
for col in range(0,(columns)/50+1):
    n = col*50
    plt.figure()
    ax = plt.subplot(111)
    ax.plot(data[:,0],data[:,n+2],label='wave',c=color[0],ls=style[1])
    ax.set_xlabel('x')
    ax.set_ylabel('amplitude')
    ax.set_ylim((-1.0,1.0))
#    str = 'Time {0:.1} ms'.format(.1*n)
#    plt.title(str)
    leg = ax.legend(loc='best', shadow=True, fancybox=True)
    fig = "../frames/frame_%05d"%(n)
    plt.savefig(fig+'.png',format='png')
    plt.clf()

command = ('mencoder',
           'mf://*png',
           '-mf',
           'type=png:w=800:h=600:fps=5',
           '-ovc',
           'lavc',
           '-lavcopts',
           'vcodec=mpeg4',
           '-oac',
           'copy',
           '-o',
           'profile_LEM.avi')
print "\n\nabout to execute:\n%s\n\n" % ' '.join(command)
subprocess.check_call(command)
command = ('convert',
           'frame*png',
           '+dither',
           '-layers',
           'Optimize',
           '-colors',
           '256',
           '-depth',
           '8',
           'profile_LEM.gif')

subprocess.check_call(command)
print "\n\n The movie was written to 'output.avi'"
