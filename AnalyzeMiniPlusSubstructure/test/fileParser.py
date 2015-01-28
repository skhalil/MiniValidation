#!/usr/bin/python

import re, math, os, sys
from optparse import OptionParser

def main():

  parser = OptionParser()

  parser.add_option('--input', metavar='F', type='string', action='store',
                    default = 'ttjetsPHYS14Files.txt',
                    dest='input',
                    help='Input text file')

  parser.add_option('--outputDir', metavar='F', type='string', action='store',
                    default = 'ttbar',
                    dest='outputDir',
                    help='output directory')

  parser.add_option('--dummyConfig', metavar='F', type='string', action='store',
                    default = "dummy_config.py",
                    dest='dummyConfig',
                    help='the template config')

  # Parse and get arguments
  (options, args) = parser.parse_args()
  inputFile = options.input
  outDir = options.outputDir
  tempConfig = options.dummyConfig
  # get total number of files stored in a text file
  with open(inputFile) as f:
    f_num = sum(1 for _ in f)
  print f_num
  # create a copy of parameter config for each file
  in_file = open(inputFile)
  for n_file in range (1, f_num+1):
     temp = open(tempConfig)
     o_file = open(outDir+'/'+outDir+'_'+str(n_file)+'.py', 'w')
     file_path = in_file.readline()
     print file_path, int(n_file)
     for line in temp:
        line = line.replace('dummy',outDir+str(n_file))
        line = line.replace('file_path',str(file_path))	
        o_file.writelines(line)
     temp.close()
     o_file.close()

if __name__ == '__main__':
    main()
