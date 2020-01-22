import shutil
import sys
import re

"""
This script intends to filter out the doxygen non-documented classes (set in EXCLUDE_SYMBOLS)
"""

# example to filter
# <tr class="memitem:"><td class="memItemLeft" align="right" valign="top">class &#160;</td><td class="memItemRight" valign="bottom"><b>ConstraintReal</b></td></tr>                                                 
# <tr class="memdesc:"><td class="mdescLeft">&#160;</td><td class="mdescRight">This class represents a real constraint on a matrix. <br /></td></tr>                                                                
# <tr class="separator:"><td class="memSeparator" colspan="2">&#160;</td></tr>

# example to keep
# <tr class="memitem:"><td class="memItemLeft" align="right" valign="top">class &#160;</td><td class="memItemRight" valign="bottom"><a class="el" href="classmatfaust_1_1factparams_1_1ParamsFact.html">ParamsFact</a></td></tr>
# <tr class="memdesc:"><td class="mdescLeft">&#160;</td><td class="mdescRight">The parent abstract class to represent the general factorization parameters.  <a href="classmatfaust_1_1factparams_1_1ParamsFact.html#details">More...</a><br /></td></tr>
# <tr class="separator:"><td class="memSeparator" colspan="2">&#160;</td></tr>


in_block_to_del=False
count_in_block_lines = 0
for file in sys.argv[1:]:
    print("filtering file:", file)
    tmp = open("tmp.html", mode='w')
    for line in open(file):
        if(in_block_to_del):
            count_in_block_lines -= 1
            if(count_in_block_lines <= 0):
                in_block_to_del = False
        elif(re.match('.*class &#160.*"bottom"><b>', line)):
            in_block_to_del=True
            count_in_block_lines = 2
        else:
            tmp.write(line)
    tmp.close()
    shutil.copyfile("tmp.html", file)
