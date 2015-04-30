__author__ = 'liangshengzhang'

import os
import fnmatch
import subprocess

count = None

count_file = "./parameters/count.txt"
if os.path.exists(count_file):
    f = open(count_file, 'r')
    line = f.readline()
    count = line.split()[0]

if count is not None:
    for file in os.listdir('./parameters'):
        if fnmatch.fnmatch(file, '*' + count + ".dat"):
            start = file.find(".dat")
            filename = file[:start]
            subprocess.call("mv ./parameters/" + file + " ./parameters/" + filename + ".temp")
        elif fnmatch.fnmatch(file, "*.dat"):
            subprocess.call("rm ./parameters/" + file)

    for file in os.listdir("./parameters"):
        if fnmatch.fnmatch(file, "*.temp"):
            start = file.find(".temp")
            filename = file[:start]
            subprocess.call("mv ./parameters/" + file + " ./parameters/" + filename + ".dat")