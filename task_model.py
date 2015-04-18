__author__ = 'liangshengzhang'

import subprocess


class TASK_MODEL(object):
    """
    Record tasks and models of the program
    """
    def __init__(self):
        self.tasks = []
        self.models = []


def task_model(tasks_models):
    """
    Build and run task_model program, and record the outputted tasks and models. Here
    it is assumed that in the output only two lines have keywords "tasks" or "models",
    and below two lines the names of tasks or models are given as the second word.
    :param tasks_models: The data structure that holds tasks and models
    :return: None
    """
    make_process = subprocess.Popen("make task_model -j4",stderr=subprocess.STDOUT, shell=True)
    if make_process.wait() != 0:
        print "Make Error"

    proc = subprocess.Popen("./task_model", stdout=subprocess.PIPE, shell=True)

    find_model = False
    find_task = False

    while True:
        line = proc.stdout.readline()
        if line != '':
            # Not end of the output

            find_content = False
            if line.find("models") > -1:
                # The line with keyword "models"
                find_model = True
                find_task = False
            elif line.find("tasks") > -1:
                # The line with keyword "tasks"
                find_model = False
                find_task = True
            elif len(line) > 1:
                # Not empty line
                find_content = True

            if find_content:
                content = line.split()

                if find_model:
                    tasks_models.models.append(content[1])
                elif find_task:
                    tasks_models.tasks.append(content[1])
        else:
            break