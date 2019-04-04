#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function

import urllib
import zipfile
from os.path import join, isdir, isfile
import tempfile
from sys import argv, version_info
if version_info[0] == 3:
    import urllib.request
from os import listdir

ARCH_NAME = "@REMOTE_DATA_FILE@"
BASE_URL = "@REMOTE_DATA_URL@"

def download_uncompress(uncompress_dir=None):
    ARCH_URL = "/".join([BASE_URL, ARCH_NAME])

    TMP_DIR = tempfile.gettempdir()
    DEST_FILE = join(TMP_DIR, ARCH_NAME)

    def reporthook(x,y,z):
        print("\rDownloading FAµST data:", int((x*y)*100.0/z),'%', end='')

    if(uncompress_dir):
        if(not isdir(uncompress_dir)):
            raise Exception(uncompress_dir+" is not an existing "
                            "directory/folder.")
        loc_files = [f for f in listdir(uncompress_dir) if isfile(join(uncompress_dir, f))]
        if(len(loc_files) > 0):
            print("It seems FAµST data is already available locally. To renew"
                  " the download please empty the directory:", uncompress_dir)
            return

    if version_info[0] == 2:
        urllib.urlretrieve (ARCH_URL, DEST_FILE, reporthook=reporthook)
    else:
        urllib.request.urlretrieve (ARCH_URL, DEST_FILE, reporthook=reporthook)

    print()

    print("====== data downloaded:", DEST_FILE)

    if uncompress_dir and isdir(uncompress_dir):
        print("Uncompressing zip archive to", uncompress_dir)
        zip_ref = zipfile.ZipFile(DEST_FILE, 'r')
        zip_ref.extractall(uncompress_dir)
        zip_ref.close()



if __name__ == '__main__':
    uncompress_dir=None
    if len(argv) > 1 and isdir(argv[1]):
        uncompress_dir=argv[1]
    download_uncompress(uncompress_dir)
