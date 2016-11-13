import fnmatch
import os
import time
import platform
import subprocess
import re

import numpy as np

import matplotlib.pyplot as plt
plt.ion()

def get_processor_name():
    if platform.system() == "Windows":
        return platform.processor()
    elif platform.system() == "Darwin":
        os.environ['PATH'] = os.environ['PATH'] + os.pathsep + '/usr/sbin'
        command ="sysctl -n machdep.cpu.brand_string"
        return subprocess.check_output(command).strip()
    elif platform.system() == "Linux":
        command = "cat /proc/cpuinfo"
        all_info = subprocess.check_output(command, shell=True).strip()
        for line in all_info.split("\n"):
            if "model name" in line:
                return re.sub( ".*model name.*:\s*", "", line,1)
    return ""

def getlinestyle(method):
    color = None
    marker = None
    if method=='ooura':
        color = 'magenta'
        marker = 's'
    if method=='ipp':
        color = 'black'
        marker = '.'
    if method=='fftw3':
        color = 'blue'
        marker = 'o'
    if method=='djbfft':
        color = 'yellow'
    if method=='ffts':
        color = 'red'
    if method=='pffft':
        color = 'green'
    if method=='fftreal':
        color = 'cyan'

    return color, marker

if  __name__ == "__main__" :
    print('Plotting benchmark results...')

    methods = dict()
    for root, dirnames, filenames in os.walk('.'):
        for fresult in fnmatch.filter(filenames, 'benchmark.log'):
            fresult = os.path.join(root,fresult)
            print('Reading '+fresult+'...')
            with open(fresult) as f:
                for line in f:
                    res = line.split()
                    method = res[0]
                    if not method in methods:
                        methods[method] = dict()
                        methods[method]['Ns'] = list()
                        methods[method]['MFlops'] = list()
                        methods[method]['dur5'] = list()
                        methods[method]['dur50'] = list()
                        methods[method]['dur95'] = list()
                        methods[method]['acc5'] = list()
                        methods[method]['acc50'] = list()
                        methods[method]['acc95'] = list()

                    methods[method]['Ns'].append(int(res[2]))
                    methods[method]['MFlops'].append(float(res[3]))

                    methods[method]['dur5'].append(float(res[4]))
                    methods[method]['dur50'].append(float(res[5]))
                    methods[method]['dur95'].append(float(res[6]))

                    methods[method]['acc5'].append(float(res[7]))
                    methods[method]['acc50'].append(float(res[8]))
                    methods[method]['acc95'].append(float(res[9]))


    # Guess the precision
    usefloat32 = np.mean(methods[methods.keys()[0]]['acc50'])>1e-10
    if usefloat32: print('Assume single precision (float32)')
    else:          print('Assume double precision (float64)')

    f, axs = plt.subplots(2, 1, sharex=True, sharey=False)
    for m in methods:
        method = methods[m]
        Ns = np.array(method['Ns'])
        method['dur5'] = np.array(method['dur5'])
        method['dur50'] = np.array(method['dur50'])
        method['dur95'] = np.array(method['dur95'])
        method['acc5'] = np.array(method['acc5'])
        method['acc50'] = np.array(method['acc50'])
        method['acc95'] = np.array(method['acc95'])

        color, marker = getlinestyle(m)
#        axs[0].plot(np.log2(method['Ns']), method['MFlops'], label=m, color=color, marker=marker)

        MFlops = 0.5*5*Ns*np.log2(Ns)/(1e6*method['dur50'])
        MFlops_up = 0.5*5*Ns*np.log2(Ns)/(1e6*method['dur5'])
        MFlops_bt = 0.5*5*Ns*np.log2(Ns)/(1e6*method['dur95'])
        axs[0].fill_between(np.log2(Ns), MFlops_bt, MFlops_up, facecolor=color, alpha=0.5)
        axs[0].plot(np.log2(Ns), MFlops, label=m, color=color, marker=marker)

        if usefloat32:
            method['acc5'][method['acc5']==0.0] = 10**-8.0
        else:
            method['acc5'][method['acc5']==0.0] = 10**-17.0

        axs[1].fill_between(np.log2(Ns), np.log10(method['acc5']), np.log10(method['acc95']), facecolor=color, alpha=0.5)
        axs[1].plot(np.log2(Ns), np.log10(method['acc50']), label=m, color=color, marker=marker)

    if usefloat32: axs[1].set_ylim((-7.6, -6.6))
    else:          axs[1].set_ylim((-16.5, -15.25))

    axs[0].legend(loc='upper left')
    axs[0].grid()
    axs[1].grid()
    axs[0].set_ylabel('Speed [MFlops]')
    axs[1].set_ylabel('Accuracy [log10 RMS]')
    axs[1].set_xlabel('FFT size N')
    if usefloat32:  title = 'Single precision (float32)'
    else:           title = 'Double precision (float64)'
    axs[1].set_title('(ranges show 5% to 95% percentiles)\n')
    f.canvas.draw()

    labels = [item.get_text() for item in axs[0].get_xticklabels()]
    for li in range(len(labels)):
        labels[li] = str(2**int(float(labels[li])))
    axs[0].set_xticklabels(labels)
    axs[1].set_xticklabels(labels)

#    labels = [item.get_text() for item in axs[1].get_yticklabels()]
#    print(labels)
#    labels = [w.replace(u'\N{MINUS SIGN}', '-') for w in labels]
#    print(labels)
#    for li in range(len(labels)):
#        labels[li] = '$10^{'+labels[li]+'}$'
#    axs[1].set_yticklabels(labels)

    # Get CPU info
    title = get_processor_name()+'\n'+title+'\n'
    axs[0].set_title(title)

    from IPython.core.debugger import  Pdb; Pdb().set_trace()

