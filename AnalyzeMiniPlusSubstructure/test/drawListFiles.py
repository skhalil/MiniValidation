#!/usr/bin/env python

import subprocess
fPath = '/eos/uscms/store/user/skhi'
options = [
    [fPath+'/Tpwb1200/MINIAOD', 'TpWb_m_1200'],
    [fPath+'/Tpwb800/MINIAOD', 'TpWb_m_800'],
    ]
command = 'python listFiles.py --path={0:s} --outputText={1:s} '

for option in options :
    s = command.format(option[0], option[1])
    subprocess.call( ["echo --------------------------------------------------------------------------",""], shell=True)
    subprocess.call( ["echo %s"%s,""], shell=True)
    subprocess.call( ["echo --------------------------------------------------------------------------",""], shell=True)

    subprocess.call( [s, ""], shell=True )
