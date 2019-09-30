from sys import argv
from glob import glob
from os.path import join


html_tags = [ '<br>', '<br/>', '<b>', '<b/>', '</b>', '<code>', '</code>' ]

if __name__ == "__main__":
    py_mods = glob(join(argv[1], '**.py'))
    for script2filter in py_mods:
        f = open(script2filter)
        lines = f.readlines()
        filtered_lines = []
        for line in lines:
            filtered_line = line
            for tag in html_tags:
                filtered_line = filtered_line.replace(tag, '')
            filtered_lines += [ filtered_line ]
        f.close()
        # overwrite with filtered lines
        outf = open(script2filter, "w")
        outf.writelines(filtered_lines)
        outf.close()

