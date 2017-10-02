#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
print (sys.version)
import time

my_dir = os.path.abspath(os.path.dirname(__file__))

# mylib/faf
my_library = os.path.expanduser('~/.pylib')
sys.path.append(my_library)
# libraries:
# from mylib.FindAllFiles import *
# from mylib.moving_average import *
from mylib.cp import *
# from mylib.FindAllFiles import *
# from mylib.highway_check import *
# from mylib.moving_average import *
# from mylib.regex import reg_ex
# from mylib.run_command import run_command


dir_temp = sys.argv[1]
x1 = int(sys.argv[2])
x2 = int(sys.argv[3])


for i in range(x1,x2+1):
    print(i)
    print(dir_temp)
    new_dirname = dir_temp.split('_')[0] + '_' + str(i).zfill(3)
    print(new_dirname)
    # print my_dir
    cp_tree(my_dir,dir_temp,my_dir,new_dirname)
