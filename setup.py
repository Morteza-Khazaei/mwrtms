from setuptools import setup, find_packages, Extension

try:
    from Cython.Build import cythonize
    import numpy as np
except ImportError:  # pragma: no cover - optional speedup
    cythonize = None
    np = None

def read_requirements(file):
    with open(file) as f:
        return f.read().splitlines()

def read_file(file):
   with open(file) as f:
        return f.read()
    
long_description = read_file('README.md')
version = read_file('VERSION')
requirements = read_requirements('requirements.txt')

ext_modules = []
if cythonize is not None and np is not None:
    extensions = [
        Extension(
            name="ssrt.surface._i2em_cy",
            sources=["src/ssrt/surface/_i2em_cy.pyx"],
            include_dirs=[np.get_include()],
        )
    ]
    ext_modules = cythonize(extensions, compiler_directives={"language_level": "3"})

setup(
    name = 'pySSRT',
    version = version,
    author = 'Morteza Khazaei',
    author_email = 'morteza.khazaei@usherbrooke.ca',
    url = 'https://github.com/Morteza-Khazaei/SSRT',
    description = 'Single Scattering Radiative Transfer model (SSRT) in Python',
    long_description_content_type = 'text/markdown',  # If this causes a warning, upgrade your setuptools package
    long_description = long_description,
    license = 'MIT license',
    package_dir={'': 'src'},
    packages = find_packages(
        where='src', 
        exclude=['dataset']
    ),  # Don't include test directory in binary distribution
    install_requires = requirements,
    ext_modules=ext_modules,
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ]  # Update these accordingly
)
