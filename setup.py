from distutils.core import setup, Extension

include_dirs = ['/usr/include', '/usr/local/include']
library_dirs = ['/usr/local/lib', ]
libraries = ['png', 'z', ]
runtime_library_dirs = []
extra_objects = []
define_macros = []
sources = ['src/errors.h',
          'src/pngnq.c',
          'src/colorspace.h',
          'src/colorspace.c',
          'src/neuquant32.h',
          'src/neuquant32.c',
          'src/nqcvector.h',
          'src/pngnqhelp.h',
          'src/rwpng.h',
          'src/rwpng.c']

setup(name = "pypngnq",
      version = "0.1.2",
      author = "Guangming Li",
      author_email = "leeful@gmail.com",
      license = "Apache License",
      description = 'a utility for optimizing PNG files',
      url = "https://github.com/cute/pypngnq",
      keywords = ['TinyPNG', 'PngNQ', 'ImageOptim', 'PngOptim'],
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

