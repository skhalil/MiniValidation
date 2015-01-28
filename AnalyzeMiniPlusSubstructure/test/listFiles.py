#!/usr/bin/python

import sys
import os, commands
from optparse import OptionParser

parser = OptionParser()

parser.add_option('--path', metavar='F', type='string', action='store',
                  default = '/eos/uscms/store/user/skhi/Tpwb1200/MINIAOD',
                  dest='path',
                  help='Input path')

parser.add_option('--outputText', metavar='F', type='string', action='store',
                  default = "outputText",
                  dest='outputText',
                  help='output file')

# Parse and get arguments
(options, args) = parser.parse_args()

path = options.path
textName = options.outputText
f = open(textName+'.txt', 'w')

cmd = 'ls %s/' %(path)

for i in commands.getoutput(cmd).split('\n'):
    print path+"/"+i
    f.write(path+"/"+i+'\n')
    
