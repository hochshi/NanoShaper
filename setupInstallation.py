#!/usr/bin/env python

from distutils.core import setup
from distutils.sysconfig import get_python_lib
import shutil 

#searching version

file = open('../src/nanoshaper.i','r')
strs = file.read().split()
ind = 0
ver = ''
for l in strs:
	if (l=='VERSION'):
		ver = strs[ind+1]
		ll = list(ver)
		p = ll.index('"')
		del(ll[p])
		p = ll.index('"')
		del(ll[p])
		ver ="".join(ll)
		break
	ind = ind+1

setup(name='NanoShaper',
      version=ver,
      description='NanoShaper',
      author='Sergio Decherchi',
      author_email='sergio.decherchi@iit.it',
      url='http://www.electrostaticszone.eu',
      license='GPL',
      py_modules = ['NanoShaper'],
      packages=['NanoShaper'],
      package_data={'NanoShaper':['_NanoShaper.so']},
     )