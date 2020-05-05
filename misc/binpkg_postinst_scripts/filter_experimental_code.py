#!/usr/bin/env python

import re, sys

EXP_BLOCK_START="^(#|%) experimental block start"
EXP_BLOCK_END="^(#|%) experimental block end"

out_lines = []
print(sys.argv[0], "executing...")
# option to not write empty out file
if("--no-empty" in sys.argv[:]):
    noempty = True
    sys.argv.remove("--no-empty")
else:
    noempty = False

if(len(sys.argv) >= 2):
    for script in [sys.argv[1]]:
        in_exp_block = False
        if sys.version_info[0] < 3:
            # python2 doesn't provide encoding argument
            fd = open(script)
        else:
            fd = open(script, encoding='utf-8')
        print("processing", script)
        for line in fd:
            if(not in_exp_block and re.match(EXP_BLOCK_START,
                                             line)):
                in_exp_block = True
            elif(in_exp_block and re.match(EXP_BLOCK_END, line)):
                in_exp_block = False
            elif(not in_exp_block):
                out_lines += [ line ]


# print("out_lines:", ''.join(out_lines))
# print("matching empty=", re.match('^\s*$',''.join(out_lines)))

if(len(sys.argv) >= 3):
    if(not noempty or not re.match('^\s*$',''.join(out_lines))):
        if sys.version_info[0] < 3:
            outf = open(sys.argv[2], 'w')
        else:
            outf = open(sys.argv[2], 'w', encoding='utf-8')
        outf.writelines(out_lines)
        outf.close()
    else:
        print("deleting file:", sys.argv[2])
        # delete the file
        import os
        os.remove(sys.argv[2])
else:
    print(out_lines)

