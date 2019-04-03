#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function

import urllib
import zipfile
from os.path import join, isdir
import tempfile
from sys import argv, version_info
if version_info[0] == 3:
    import urllib.request


ARCH_NAME = "faust_data-2.4.2.zip"
BASE_URL = "https://gforge.inria.fr/frs/download.php/file/37960/"

def download_uncompress(uncompress_dir=None):
    ARCH_URL = join(BASE_URL, ARCH_NAME)

    TMP_DIR = tempfile.gettempdir()
    DEST_FILE = join(TMP_DIR, ARCH_NAME)

    def reporthook(x,y,z):
        print("\rDownloading FAµST data:", int((x*y)*100.0/z),'%', end='')

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
