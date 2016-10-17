from distutils.core import setup

setup(
    name='hamopy',
    version='0.2.0',
    author='S. Rouchier',
    author_email='s.rouchier@gmail.com',
    packages=['hamopy', 'hamopy.materials', 'hamopy.benchmarks'],
    scripts=['bin/BM3_simul.py'],
    url='http://pypi.python.org/pypi/hamopy/',
    license='LICENSE.txt',
    description='Heat, air and moisture transfer modelling in python',
    long_description=open('README.txt').read(),
    #install_requires=["numpy >= 1.7.1","pandas >= 0.11.0",],
)