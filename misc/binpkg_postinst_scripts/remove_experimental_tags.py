from sys import argv
import re

EXP_BLOCK_START="^(#|%) experimental block start"
EXP_BLOCK_END="^(#|%) experimental block end"

if __name__ == '__main__':
    if len(argv) > 1:
        fp = argv[1]
        with open(fp) as f:
            lines = f.readlines();
        with open(fp, 'w+') as f:
            for l in lines:
                if not re.match(EXP_BLOCK_START, l) and not re.match(EXP_BLOCK_END, l):
                    f.write(l)

