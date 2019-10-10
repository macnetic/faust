#!/usr/bin/env python

import re, sys

EXP_BLOCK_START="^(#|%) experimental block start"
EXP_BLOCK_END="(#|%) experimental block end"

if(len(sys.argv) >= 1):
    for script in sys.argv[1:]:
        in_exp_block = False
        for line in open(sys.argv[1]):
            if(not in_exp_block and re.match(EXP_BLOCK_START,
                                             line)):
                in_exp_block = True
            elif(in_exp_block and re.match(EXP_BLOCK_END, line)):
                in_exp_block = False
            elif(not in_exp_block):
                print(line,end='')

