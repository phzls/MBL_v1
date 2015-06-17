__author__ = 'liangshengzhang'

import subprocess
import task_model as tm
import file_process as fp
import core as core
import os
import fnmatch

data = fp.DATA()
tasks_models = tm.TASK_MODEL()
exe_folder = "executables" # The folder name for executables

success_run = False # Check whether the program has successfully run

# Obtain system
proc = subprocess.Popen("uname", stdout=subprocess.PIPE, shell=True)
system = proc.communicate()[0]

if system == "Linux\n":
    subprocess.call("module load intel", stdout=subprocess.PIPE, shell=True)

    # Construct a folder for executables if it does not exist
    if os.path.isdir("./" + exe_folder) is False:
        subprocess.call("mkdir " + exe_folder, stdout=subprocess.PIPE, shell=True)

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

# In case any exception is thrown, the program can end gracefully
try:
    # Obtain tasks and models
    tm.task_model(tasks_models)

    # Generate generic parameters
    fp.para_gen("generic_para", tasks_models, data, count)

    file_modify["generic_para"] = True

    """
    print "Tasks:"
    for n in tasks_models.tasks:
        print n
    print "Models:"
    for n in tasks_models.models:
        print n
        m = n.lower()
        model_name = ' '.join(n.split())
        f = open("./parameters/" + m+".dat",'w')
        f.write("//" + model_name + "Parameters\n")
        f.write("//\n")
        f.write("J 0.6 // Disorder Parameter\n")
        f.close()
    """

    print fp.bcolors.BOLD, "System size is " + str(data.size) + "\n", fp.bcolors.ENDC

    # Ask whether need to change output data
    output_a = False
    while output_a is False:
        answer = raw_input("Modify output parameters?\n")
        if (answer == "Yes") or (answer == "yes") or (answer == "Y") or (answer == "y"):
            fp.para_gen("output_para", tasks_models, data, count)
            output_a = True
            file_modify["output_para"] = True
        elif (answer == 'No') or (answer == 'no') or (answer == 'N') or (answer == 'n'):
            output_a = True
        else:
            print fp.bcolors.FAIL + "Answer must be Yes(Y) or No(N)." + fp.bcolors.ENDC

    # Generate Model data
    fp.para_gen(data.model.lower(), tasks_models, data, count)
    file_modify[data.model.lower()] = True

    # Generate Task data
    fp.para_gen(data.task.lower(), tasks_models, data, count)
    file_modify[data.task.lower()] = True

    # Update all files which are not modified to have correct file names
    new_match_name = "_" + str(count) + ".dat"
    for key in file_modify:
        if file_modify[key] is False:
            subprocess.call("cp ./parameters/" + key + match_name
                            + " ./parameters/" + key + new_match_name, shell=True )
            file_modify[key] = True

    # Handle the executable
    temp_success_run = core.core(system=system, count=count, file_modify=file_modify,
                                 exe_folder=exe_folder, data=data)
    success_run = success_run and temp_success_run

except Exception:
    print "An error occurred."
    fp.file_clean(count, file_modify)
    raise

if success_run:
    file = "./parameters/count.txt"
    file_new = "./parameters/count_temp.txt"
    if count > 1:
        subprocess.call("rm " + file, shell=True)
    subprocess.call("mv " + file_new + " " + file, shell=True)

else:
    fp.file_clean(count, file_modify)





