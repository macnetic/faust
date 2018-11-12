import re
import sys

if(len(sys.argv) < 2):
    sys.exit()

while(True):
    try:
        line = input();
        for ns in sys.argv[1:]:
            if(re.match(r'^.*@namespace\s+'+ns,line)):
                break
        else:
            print(line);
    except:
        break

