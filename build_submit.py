__author__ = 'liangshengzhang'

import subprocess
import task_model as tm
import file_process as fp
import submit as sb

data = fp.DATA()
tasks_models = tm.TASK_MODEL()

# Obtain tasks and models
tm.task_model(tasks_models)

# Generate parameters
fp.para_gen("generic_para", tasks_models, data)
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


# Ask whether need to change output data
output_a = False
while output_a is False:
    answer = raw_input("Modify output parameters?\n")
    if (answer == "Yes") or (answer == "yes") or (answer == "Y") or (answer == "y"):
        fp.para_gen("output_para", tasks_models, data)
        output_a = True
    elif (answer == 'No') or (answer == 'no') or (answer == 'N') or (answer == 'n'):
        output_a = True
    else:
        print "Answer must be Yes(Y) or No(N)."

# Generate Model data
fp.para_gen(data.model.lower(), tasks_models, data)

# Generate Task data
fp.para_gen(data.task.lower(), tasks_models, data)

proc = subprocess.Popen("uname", stdout=subprocess.PIPE, shell=True)
system = proc.communicate()[0]

if system == "Darwin\n":
    # It is a mac system
    cont = raw_input("Continue? ")
    if (cont == '') or cont.startswith('y') or cont.startswith('Y'):
        make_process = subprocess.Popen("make auto OUT=mbl_auto -j4",stderr=subprocess.STDOUT, shell=True)
        if make_process.wait() != 0:
            print "Make Error"

        run = subprocess.Popen("./mbl_auto", stdout=subprocess.PIPE, shell=True)
        print run.communicate()[0]

elif system == "Linux\n":
    # Make the program
    progname = data.model.lower() + "_" + data.size + "_" + data.task.lower()
    make_process = subprocess.Popen("make auto OUT=" + progname
                                    + " -j4",stderr=subprocess.STDOUT, shell=True)
    if make_process.wait() != 0:
        print "Make Error"

    valid_choice = False
    while valid_choice is False:
        choice = raw_input("Submit (s), Direct Run (r), Valgrind (v) or Exit (e)?\n")

        if choice.startswith("r") or choice.startswith("R"):
            run = subprocess.Popen("./" + progname, stdout=subprocess.PIPE, shell=True)
            print run.communicate()[0]

        elif choice.startswith("v") or choice.startswith("R"):
            if int(data.num_threads) > 1:
                print "Valgrind has problems handling multi-thread."
            else:
                val = subprocess.Popen("valgrind --leak-check=full ./" + progname,
                                       stdout=subprocess.PIPE, shell=True)
                out = val.communicate()[0]

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
                sb.adroit_submit(data, progname)

            elif (server == "Feynman") or (server == "feynman"):
                sb.feynman_submit(data, progname)

            valid_choice = True

        else:
            print fp.bcolors.FAIL + "Invalid choice." << fp.bcolors.ENDC
else:
    print fp.bcolors.FAIL + "Unknown system type." << fp.bcolors.ENDC

sb.feynman_submit(data,"aa")






