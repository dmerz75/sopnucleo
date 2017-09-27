#!/usr/bin/env python
import os,sys,time

my_library = os.path.expanduser('~/.pylib')
sys.path.append(my_library)

# imports from my_library
from cp import *
from regex import *


my_dir = os.path.abspath(os.path.dirname(__file__))

dest_dir = os.path.join(my_dir,'sopnuc1')
if not os.path.exists(dest_dir): os.makedirs(dest_dir)
struct_dir = os.path.join(my_dir,'sopnuc1/Struct_data')
if not os.path.exists(struct_dir): os.makedirs(struct_dir)
coord_dir  = os.path.join(my_dir,'sopnuc1/Coord')
if not os.path.exists(coord_dir): os.makedirs(coord_dir)

# submitf  = os.path.join(my_library,'py_qsub/submit.py')
cp_file(my_library,'py_qsub/submit.py',dest_dir,'submit.py')
cp_file(my_library,'py_qsub/template_sopnucleo.pbs',\
        dest_dir,'template_sopnucleo.pbs')
# template = os.path.join(my_library,'py_qsub/template_sopnucleo.pbs')


cp_file(my_dir,'run_protein',dest_dir,'run_protein')
cp_file(my_dir,'pdbADP4B9Q.pdb',dest_dir,'pdbADP4B9Q.pdb')
cp_file(my_dir,'Contacts_intra_inter_b_ADP4B9Q_j29',dest_dir,'Contacts_intra_inter_b_ADP4B9Q_j29')
cp_file(my_dir,'def_param.h',dest_dir,'def_param.h')


