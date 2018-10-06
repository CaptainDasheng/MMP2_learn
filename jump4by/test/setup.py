# coding: utf-8
# Copyright (c) JUMP2 Development Team.
# Distributed under the terms of the JLU License.


#=================================================================
# This file is part of JUMP2.
#
# Copyright (C) 2017 Jilin University
#
#  Jump2 is a platform for high throughput calculation. It aims to 
#  make simple to organize and run large numbers of tasks on the 
#  superclusters and post-process the calculated results.
#  
#  Jump2 is a useful packages integrated the interfaces for ab initio 
#  programs, such as, VASP, Guassian, QE, Abinit and 
#  comprehensive workflows for automatically calculating by using 
#  simple parameters. Lots of methods to organize the structures 
#  for high throughput calculation are provided, such as alloy,
#  heterostructures, etc.The large number of data are appended in
#  the MySQL databases for further analysis by using machine 
#  learning.
#
#  Jump2 is free software. You can redistribute it and/or modify it
#  under the terms of the GNU General Public License as published 
#  by the Free sofware Foundation, either version 3 of the License,
#  or (at your option) and later version.
# 
#  You should have recieved a copy of the GNU General Pulbic Lincense
#  along with Jump2. If not, see <https://www.gnu.org/licenses/>.
#=================================================================

import shutil
from setuptools import find_packages, setup
import os, sys


#required packages.

#setup (install_requires=['ipython >= 2.0'])
#setup (install_requires=['django  >= 1.7'])
#setup (install_requires=['python  >= 2.7.0'])
#setup (install_requires=['numpy   >= 1.6.0'])

# optional required packages.

#setup (install_requires=['numpy >= 1.6.0'])


# jump2 packages.
#JUMP2 = 'module'
#JUMP2 = 'src'
JUMP2 = "jump2"
PACKAGES = [JUMP2] + ['%s.%s' % (JUMP2, i) \
            for i in find_packages(JUMP2)]

setup(
      name=JUMP2,
      version='V1.0',
      author='Lijun Zhang, Xin-Gang Zhao, Yuhao Fu, Jian Lv, \
              Guangren Na, Shulin Luo, Bangyu Xing, Xin He, 
              Tianshu Li, Dongwen Yang, Yuanhui Sun, Yawen Li,
              Qiaoling Xu, Zhun Liu, Dianlong Zhao, Xinjiang Wang',
      author_email='lijun_zhang@jlu.edu.cn'
      maintainer='Xin-Gang Zhao, Yuhao Fu, Shulin Luo, Guangren Na',
      maintainer_email='xgzhao0201@gmail.com',
      license = 'Jilin University',
      description = 'Python Package for High Throughput Materials-design \
                     calculation.',
      platforms='Linux',
      install_requires=['ipython >= 2.0', 'django >= 1.7', 'numpy >= 1.6.0'],
      packages = PACKAGES)


#
#if not os.path.exists('/home/gordon/.local/jump2/bin/jump2'):

#    os.makedirs('/home/gordon/.local/jump2/bin')
#    shutil.copy('jump2.x', '/home/gordon/.local/jump2/bin/jump2')
#    os.system('export PATH=/home/gordon/.local/jump2/bin:$PATH')


