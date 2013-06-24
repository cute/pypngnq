from distutils.core import setup, Extension

include_dirs = ['/usr/include', '/usr/local/include']
library_dirs = ['/usr/local/lib', ]
libraries = ['png', 'z', ]
runtime_library_dirs = []
extra_objects = []
define_macros = []
sources = ['src/pngnq.c',
          'src/colorspace.c',
          'src/neuquant32.c',
          'src/rwpng.c']

setup(name = "pypngnq",
      version = "0.1.1",
      author = "Guangming Li",
      author_email = "leeful@gmail.com",
      license = "Apache License",
      url = "http://liguangming.com/pngnq",
      packages = ["pypngnq"],
      ext_package = "pypngnq",
      ext_modules = [Extension(name = "pngnq",
                                sources = sources,
                                include_dirs = include_dirs,
                                library_dirs = library_dirs,
                                runtime_library_dirs = runtime_library_dirs,
                                libraries = libraries,
                                extra_objects = extra_objects,
                                define_macros = define_macros,
                                extra_compile_args = ['-std=c99', ],
                              )],
      )

