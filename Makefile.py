import sys
import os
import numpy
from numpy.distutils import fcompiler
from numpy import f2py
#distutils: define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
os.environ['CFLAGS']='-Ofast -funroll-all-loops -ffast-math -march=native -pipe'

source = open('src/FortranFunctions.f90').read()
modulename = 'FortranFunctions'

fcompileropts = ["--opt=-Ofast -funroll-all-loops -ffast-math -march=native -pipe"]
check = f2py.compile(source=source,modulename=modulename,extra_args=fcompileropts,verbose=True,extension='.f90',)

if check>0:
    print(f'Error during build')
    sys.exit()

files = os.listdir()
for f in files:
    fname = f.split('.')
    if (fname[0] == 'FortranFunctions')&(fname[-1]=='so'):
        libname = f
        if check==0:
            print('Build completed succesfully!')
        break

if not os.path.exists('libs'):
    os.mkdir('libs')

os.rename(libname,'libs/'+libname)
files = os.listdir('libs')
for f in files:
    fname = f.split('.')
    if (fname[0] == 'FortranFunctions')&(fname[-1]=='so'):
        print(f'Moved {f:s} to libs/{f:s}')
