# this script is temporary, TODO: delete it when
# gen_matlab_inline_doc_from_doxy_blocks will handle matlab function files in
# addition to class files
import sys, shutil, os.path, os
if(len(sys.argv) < 2):
    print("error: two arguments are needed, the input and ouput files.")
    exit(1)
IN_FILE=sys.argv[1]
OUT_FILE=sys.argv[2]
OUT_DIR=os.path.dirname(OUT_FILE)
if not os.path.exists(IN_FILE):
    exit(0)
    # raise Exception(IN_FILE+ " not found.") # not an error because some files
    # are totally experimental code and deleted
if not os.path.exists(OUT_DIR):
    os.mkdir(OUT_DIR)
shutil.copy(IN_FILE, OUT_FILE)
