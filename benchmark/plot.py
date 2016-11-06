import fnmatch
import os
import time

import numpy as np

import matplotlib.pyplot as plt
plt.ion()

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

        method['acc5'][method['acc5']==0.0] = 10**-8.0
        axs[1].fill_between(np.log2(Ns), np.log10(method['acc5']), np.log10(method['acc95']), facecolor=color, alpha=0.5)
        axs[1].plot(np.log2(Ns), np.log10(method['acc50']), label=m, color=color, marker=marker)

    axs[1].set_ylim((-7.6, -6.6))
    axs[0].legend(loc='upper left')
    axs[0].grid()
    axs[1].grid()
    axs[0].set_ylabel('Speed [MFlops]')
    axs[1].set_ylabel('Accuracy [log10 RMS]')
    axs[1].set_xlabel('FFT size N')
    axs[0].set_title('(range shows 5% to 95% percentiles)')
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

    from IPython.core.debugger import  Pdb; Pdb().set_trace()

