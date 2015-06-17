__author__ = 'liangshengzhang'

import subprocess
import file_process as fp
import submit as sb

def core(system, count, file_modify, exe_folder, data):
    """
    This function executes different commands regarding an executable on different systems. So far it only
    accepts either OSX or adroit/feynman cluster system
    :param system: The type of the system
    :param count: The count number for file execution
    :param file_modify: A map recording which files have been modified
    :param exe_folder: The name of the folder where the executable lies
    :param data: The data that contains some information about the execution
    :return: A boolean specifies whether the run has been successful
    """
    success_run = False
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
        progname = data.model.lower() + "_" + data.size + data.model_para + "_" \
                   + data.task.lower() + data.task_para
        make_process = subprocess.Popen("make auto OUT=" + progname
                                        + " -j4",stderr=subprocess.STDOUT, shell=True)
        if make_process.wait() != 0:
            fp.file_clean(count, file_modify)
            raise Exception("Make Error")

        # Move the executable to the folder
        subprocess.call("mv " + progname + "./" + exe_folder, stdout=subprocess.PIPE, shell=True)

        valid_choice = False
        while valid_choice is False:
            choice = raw_input("Submit (s), Direct Run (r), Valgrind (v) or Exit (e)?\n")

            if choice.startswith("r") or choice.startswith("R"):
                run = subprocess.Popen("./" + exe_folder+ '/' + progname + " " + str(count),
                                       stdout=subprocess.PIPE, shell=True)
                stdout = []
                while True:
                    line = run.stdout.readline()
                    stdout.append(line)
                    print line,
                    if line == '' and run.poll() != None:
                        break
                valid_choice = True
                success_run = True

            elif choice.startswith("v") or choice.startswith("V"):
                if int(data.num_threads) > 1:
                    print "Valgrind has problems handling multi-thread."
                else:
                    val = subprocess.Popen("valgrind --leak-check=full ./" + exe_folder + '/' + progname
                                           + " " + str(count), stdout=subprocess.PIPE, shell=True)
                    stdout = []
                    while True:
                        line = val.stdout.readline()
                        stdout.append(line)
                        print line
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
                    sb.adroit_submit(data, progname, count, exe_folder)

                elif (server == "Feynman") or (server == "feynman"):
                    sb.feynman_submit(data, progname, count, exe_folder)

                valid_choice = True
                success_run = True

            else:
                print fp.bcolors.FAIL + "Invalid choice." + fp.bcolors.ENDC
    else:
        print fp.bcolors.FAIL + "Unknown system type." + fp.bcolors.ENDC

    return success_run
