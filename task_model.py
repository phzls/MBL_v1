__author__ = 'liangshengzhang'

import subprocess
import file_process as fp


class TASK_MODEL(object):
    """
    Record tasks, models and initial conditions of the program
    """
    def __init__(self):
        self.tasks = []
        self.models = []
        self.inits = []


def task_model(tasks_models):
    """
    Build and run task_model program, and record the outputted tasks and models. Here
    it is assumed that in the output only three lines have keywords "tasks", "models"
    or "Initial State Construction Functions", and below these lines the names of tasks,
    models and init_func are given as the second word.
    :param tasks_models: The data structure that holds tasks, models and init_func.
    :return: None
    """
    make_process = subprocess.Popen("make task_model -j4",stderr=subprocess.STDOUT, shell=True)
    if make_process.wait() != 0:
        fp.file_clean()
        raise Exception("Make Error")

    proc = subprocess.Popen("./task_model", stdout=subprocess.PIPE, shell=True)

    find_model = False
    find_task = False
    find_init = False

    while True:
        line = proc.stdout.readline()
        if line != '':
            # Not end of the output

            find_content = False
            if line.find("models") > -1:
                # The line with keyword "models"
                find_model = True
                find_task = False
                find_init = False
            elif line.find("tasks") > -1:
                # The line with keyword "tasks"
                find_model = False
                find_task = True
                find_init = False
            elif line.find("Initial State Construction Functions") > -1:
                # The line with keyword "Initial State Construction Functions"
                find_model = False
                find_task = False
                find_init = True
            elif line.find("Initial Density Construction Functions") > -1:
                find_model = False
                find_task = False
                find_init = True
            elif len(line) > 1:
                # Not empty line
                find_content = True

            if find_content:
                content = line.split()

                if find_model:
                    tasks_models.models.append(content[1])
                elif find_task:
                    tasks_models.tasks.append(content[1])
                elif find_init:
                    tasks_models.inits.append(content[0])
        else:
            break
    tasks_models.inits = list(set(tasks_models.inits))