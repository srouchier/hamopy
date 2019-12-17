import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='hamopy',
    version='0.4.0',
    author='S. Rouchier',
    author_email='simon.rouchier@univ-smb.com',
    packages=setuptools.find_packages(),
    url='http://pypi.python.org/pypi/hamopy/',
    license='LICENSE.txt',
    description='Heat, air and moisture transfer modelling in python',
    long_description=long_description,
    #install_requires=["numpy >= 1.7.1","pandas >= 0.11.0",],
)