import os
import glob

debug = ARGUMENTS.get('debug',0)
compiler = ARGUMENTS.get('compiler','clang++')

linkflags = ''
if compiler=='g++':
    cflags = ' -Wfatal-errors '
    if int(debug):
        cflags += ' -Og -g -pg '
        linkflags = ' -g -pg '
    else:
        cflags += ' -O3 '
elif compiler=='clang++':
    cflags = '-Weverything -Werror -ferror-limit=1 -Wno-error=padded'
    if int(debug):
        cflags += ' -O0 -g -pg'
        linkflags = ' -g -pg'
    else:
        cflags += ' -O3 -march=native'
elif compiler=='mpiCC':
    cflags = ' '
    if int(debug):
        cflags += ' -O0 -g -pg -Wfatal-errors -DRICH_MPI '
        linkflags = ' -g -pg '
    else:
        cflags += ' -O3 -Wfatal-errors -DRICH_MPI '
else:
    raise NameError('unsupported compiler')

if int(debug):
    f90flags = ' -g -Og '
else:
    f90flags = ' -O3 '

build_dir = 'build/'+compiler
if int(debug):
    build_dir+='/debug'
else:
    build_dir+='/release'

VariantDir(build_dir,'source')
env = Environment(ENV = os.environ,
                  CXX=compiler,
                  CPPPATH=[os.environ['RICH_ROOT']+'/source',
                           os.environ['RICH_ROOT']],
                  LIBPATH=[os.environ['RICH_ROOT']+'/'+build_dir,
                           os.environ['HDF5_LIB_PATH']],
                  LIBS=['rich','hdf5','hdf5_cpp','gfortran'],
                  LINKFLAGS=linkflags,
                  F90FLAGS=f90flags,
                  CXXFLAGS=cflags)
                  
env.Program(build_dir+'/rich',
            Glob(build_dir+'/*.cpp')+Glob(build_dir+'/*.f90'))
