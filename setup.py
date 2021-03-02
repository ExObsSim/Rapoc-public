import codecs
import os

from setuptools import setup, find_packages

packages = find_packages(exclude=('tests', 'docs'))

provides = ['rapoc', ]

install_requires = ['astropy',
                    'matplotlib',
                    'numpy',
                    'scipy',
                    'h5py',
                    'molmass'
                    ]

classifiers = [
    'Development Status :: 4 - Beta',
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: BSD License',
    'Operating System :: POSIX',
    'Operating System :: POSIX :: Linux',
    'Operating System :: Unix',
    'Operating System :: OS Independent',
    'Programming Language :: Python :: 3',
    'Topic :: Scientific/Engineering',
    'Topic :: Software Development :: Libraries',
]


def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()


def get_info(rel_path, info):
    for line in read(rel_path).splitlines():
        if line.startswith('__%s__' % info):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find %s string." % info)


with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(name='rapoc',
      provides=provides,
      version=get_info("rapoc/__version__.py", 'version'),
      description=get_info("rapoc/__about__.py", 'summary'),
      url=get_info("rapoc/__about__.py", 'url'),
      author=get_info("rapoc/__about__.py", 'author'),
      author_email=get_info("rapoc/__about__.py", 'email'),
      license=get_info("rapoc/__about__.py", 'license'),
      long_description=long_description,
      long_description_content_type="text/markdown",
      packages=packages,
      classifiers=classifiers,
      install_requires=install_requires,
      include_package_data=True,
      python_requires='>=3.8',
      zip_safe=False)
