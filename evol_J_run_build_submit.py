__author__ = 'liangshengzhang'

import subprocess
import task_model as tm
import file_process as fp
import submit as sb
import os
import fnmatch

data = fp.DATA()
tasks_models = tm.TASK_MODEL()

modify_words = ["size", "num_realizations", "J", "init_func_name", "time_step", "log_time"]

success_run = False # Check whether the program has successfully run

# Obtain system
proc = subprocess.Popen("uname", stdout=subprocess.PIPE, shell=True)
system = proc.communicate()[0]

if system == "Linux\n":
    subprocess.call("module load intel", stdout=subprocess.PIPE, shell=True)

# Obtain count and update count file
count = None
count = fp.count_obtain()

print fp.bcolors.BOLD, "\nThis is " + str(count) + "th running.\n", fp.bcolors.ENDC

file_modify = {} # Check whether files have been modified so that new files have been created
if count > 1:
    match_name = "_" + str(count-1) + ".dat" # A suffix for files
else:
    match_name = ".dat"
for file in os.listdir('./parameters'):
    if fnmatch.fnmatch(file, "*" + match_name):
        start = file.find(match_name)
        file_modify[file[:start]] = False


# Obtain tasks and models
tm.task_model(tasks_models)

# Generate parameters
fp.para_gen("generic_para", tasks_models, data, count, modify_words)
file_modify["generic_para"] = True

# Generate Model data
fp.para_gen(data.model.lower(), tasks_models, data, count, modify_words)
file_modify[data.model.lower()] = True

# Generate Task data
fp.para_gen(data.task.lower(), tasks_models, data, count, modify_words)
file_modify[data.task.lower()] = True

print fp.bcolors.BOLD, "System size is " + str(data.size) + "\n", fp.bcolors.ENDC

# Update all files which are not modified to have correct file names
new_match_name = "_" + str(count) + ".dat"
for key in file_modify:
    if file_modify[key] is False:
        subprocess.call("cp ./parameters/" + key + match_name
                        + " ./parameters/" + key + new_match_name, shell=True )
        file_modify[key] = True

# Ask whether user want running details
run_detail = None
while run_detail is None:
    detail = raw_input("Show running details?")
    if (detail == '') or detail.lower().startswith('n'):
        run_detail = False
    elif detail.lower().startswith('y'):
        run_detail = True
        fp.display_file("generic_para", data, count)
        fp.display_file(data.model.lower(), data, count)
        fp.display_file(data.task.lower(), data, count)

if system == "Darwin\n":
    # It is a mac system
    cont = raw_input("Continue? ")
    if (cont == '') or cont.startswith('y') or cont.startswith('Y'):
        make_process = subprocess.Popen("make auto OUT=mbl_auto -j4",stderr=subprocess.STDOUT, shell=True)
        if make_process.wait() != 0:
            fp.file_clean(count, file_modify)
            raise Exception("Make Error")

        run = subprocess.Popen("./mbl_auto " + str(count), stdout=subprocess.PIPE, shell=True)
        stdout = []
        while True:
            line = run.stdout.readline()
            stdout.append(line)
            print line,
            if line == '' and run.poll() != None:
                break
        success_run = True

elif system == "Linux\n":
    # Make the program
    progname = data.model.lower() + "_" + data.size + "_" + data.task.lower()
    make_process = subprocess.Popen("make auto OUT=" + progname
                                    + " -j4",stderr=subprocess.STDOUT, shell=True)
    if make_process.wait() != 0:
        fp.file_clean(count, file_modify)
        raise Exception("Make Error")

    valid_choice = False
    while valid_choice is False:
        choice = raw_input("Submit (s), Direct Run (r), Valgrind (v) or Exit (e)?\n")

        if choice.startswith("r") or choice.startswith("R"):
            run = subprocess.Popen("./" + progname + " " + str(count), stdout=subprocess.PIPE, shell=True)
            stdout = []
            while True:
                line = run.stdout.readline()
                stdout.append(line)
                print line,
                if line == '' and run.poll() != None:
                    break
            valid_choice = True
            success_run = True

        elif choice.startswith("v") or choice.startswith("R"):
            if int(data.num_threads) > 1:
                print "Valgrind has problems handling multi-thread."
            else:
                val = subprocess.Popen("valgrind --leak-check=full ./" + progname + " " + str(count),
                                       stdout=subprocess.PIPE, shell=True)
                stdout = []
            while True:
                line = val.stdout.readline()
                stdout.append(line)
                print line,
                if line == '' and val.poll() != None:
                    break
                valid_choice = True
                success_run = True

        elif choice.startswith("e") or choice.startswith("E") or choice.startswith('q') or choice.startswith('Q'):
            # Exit
            valid_choice = True

        elif choice.startswith("s") or choice.startswith("S"):

            # Determine server type
            f = open("submit_template.run",'r')
            for line in f:
                if line.startswith("#Server"):
                    server = line.split()[1]
                    break
                elif line.startswith("# Server"):
                    server = line.split()[2]
                    break

            if (server == "Adroit") or (server == "adroit"):
                sb.adroit_submit(data, progname, count)

            elif (server == "Feynman") or (server == "feynman"):
                sb.feynman_submit(data, progname, count)

            valid_choice = True
            success_run = True

        else:
            print fp.bcolors.FAIL + "Invalid choice." + fp.bcolors.ENDC
else:
    print fp.bcolors.FAIL + "Unknown system type." + fp.bcolors.ENDC

if success_run:
    file = "./parameters/count.txt"
    file_new = "./parameters/count_temp.txt"
    if count > 1:
        subprocess.call("rm " + file, shell=True)
    subprocess.call("mv " + file_new + " " + file, shell=True)

else:
    fp.file_clean(count, file_modify)







