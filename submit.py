__author__ = 'liangshengzhang'

import file_process as fp
import subprocess

def adroit_submit(data, program_name, count):
    T = fp.RunTime()
    # Obtain estimated running time
    non_neg = False

    while non_neg is False:
        T.hour = raw_input("Running wall time -- hour: ")
        if T.hour.isdigit() is False:
            print fp.bcolors.FAIL + "Hour must be non-negative integer." + fp.bcolors.ENDC
        else:
            non_neg = True

    non_neg = False
    while non_neg is False:
        T.min = raw_input("Running wall time -- minutes: ")
        if T.min.isdigit() is False:
            print fp.bcolors.FAIL + "Minutes must be non-negative integer." + fp.bcolors.ENDC
        elif int(T.min) > 60:
            print fp.bcolors.FAIL + "Minutes must be within 60." + fp.bcolors.ENDC
        else:
            non_neg = True

    # Generate submitting files
    submit_program_name = fp.slurm_file_gen(data, T, count)

    if submit_program_name != program_name:
        print "Program names do not match:"
        print "Name in make: " + program_name
        print "Name in slurm: " + submit_program_name
        raise Exception("Name Error")

    sub = subprocess.Popen("sbatch " + program_name + ".run", stdout=subprocess.PIPE, shell=True)
    print(sub.communicate()[0])

def feynman_submit(data, program_name, count):
    T = fp.RunTime()
    # Obtain estimated running time
    non_neg = False

    while non_neg is False:
        T.hour = raw_input("Running cpu time -- hour: ")
        if T.hour.isdigit() is False:
            print fp.bcolors.FAIL + "Hour must be non-negative integer." + fp.bcolors.ENDC
        else:
            non_neg = True

    non_neg = False
    while non_neg is False:
        T.min = raw_input("Running cpu time -- minutes: ")
        if T.min.isdigit() is False:
            print fp.bcolors.FAIL + "Minutes must be non-negative integer." + fp.bcolors.ENDC
        elif int(T.min) > 60:
            print fp.bcolors.FAIL + "Minutes must be within 60." + fp.bcolors.ENDC
        else:
            non_neg = True

    # Obtain memeory
    mem_valid = False
    while mem_valid is False:
        T.mem = raw_input("Memory -- XXmb or XXgb: ")

        if T.mem[-2:] != 'mb' and T.mem[-2:] != 'gb' and T.mem[-2:] != 'MB' and T.mem[-2:] != 'GB':
            print fp.bcolors.FAIL + "Momory must end in mb or gb." + fp.bcolors.ENDC
        elif T.mem[:-2].isdigit() is False:
            print fp.bcolors.FAIL + "Memory must be a positive integer." + fp.bcolors.ENDC
        elif (T.mem[-2:].lower() == "gb") and (int(T.mem[:-2]) > 120):
            print fp.bcolors.FAIL + "Memory must be smaller than 120gb." + fp.bcolors.ENDC
        else:
            mem_valid = True

    # Generate submitting files
    submit_program_name = fp.qsub_file_gen(data, T, count)

    if submit_program_name != program_name:
        print "Program names do not match:"
        print "Name in make: " + program_name
        print "Name in qsub: " + submit_program_name
        raise Exception("Name Error")

    sub = subprocess.Popen("qsub " + program_name + ".run", stdout=subprocess.PIPE, shell=True)
    print(sub.communicate()[0])
