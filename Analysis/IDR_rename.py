import os
import fnmatch

directory="/Users/pui/Documents/Lab_Vienna/ATAC/IDR/stage_20+ *"

for filename in os.listdir(directory):
    if fnmatch.fnmatch(filename, 'stage_2*_idr.dms'):
        print(filename)