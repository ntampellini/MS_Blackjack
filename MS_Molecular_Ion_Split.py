###################### MS MOLECULAR PEAK SPLIT PATTERN PREDICTOR ################################

formula = 'C14H9Cl5'     # Spaceless format, case-sensitive. Ex: 'C14H9Cl5'
deuterium = False        # Utilizes two isotopes for Hydrogen. EXPENSIVE!
plot = True              # Generate a plot that shows molecular ion splitting pattern
nanopeaks = False        # Plot and display minuscule peaks below 0.05 % relative abundance

#################################################################################################
from isotopic_abundances import pt
import matplotlib.pyplot as plt
import numpy as np
import time
import re

t0 = time.time()

if not deuterium:
    pt['H'] = {1:100}

flattened = []                                  #parser for chemical formula interpretation
def flatten(l):
    for i in l: 
        if type(i) == list: 
            flatten(i) 
        else: 
            flattened.append(i) 
    return flattened

parsed = re.findall('[A-Z][^A-Z]*', formula)
parsed = [re.split('(\d+)',i) for i in parsed]
parsed = flatten(parsed)
for i in parsed:
    if i == '':
        parsed.remove(i)
# print(parsed)

for i in range(len(parsed)):                    #adding '1' for single atoms like 'Cl' and 'O' in 'C2H3OCl'
    if i != len(parsed) - 1:
        if not parsed[i].isdigit() and not parsed[i+1].isdigit():
            parsed.insert(i+1, '1')
    else:
        if not parsed[i][0].isdigit():
            parsed.insert(i+1, '1')
# print(parsed)
# exit()

atoms=[]
for i in range(len(2*parsed)):                    #adding '1' for single atoms like 'Cl' and 'O' in 'C2H3OCl'
    if i != len(parsed) - 1:
        try:
            if not parsed[i].isdigit() and not parsed[i+1].isdigit():
                parsed.insert(i+1, '1')
        except:
            pass
    else:
        if not parsed[i].isdigit():
            try:
                parsed.insert(i+1, '1')
            except:
                pass

        for p in range(len(parsed)):
            if p % 2 == 0:
                atoms.append(parsed[p])
        abundance=[]
        for p in range(len(parsed)):
            if p % 2 != 0:
                abundance.append(int(parsed[p]))
        atoms_verbose = []
        for atom in range(len(atoms)):
            for _ in range(abundance[atom]):
                atoms_verbose.append(atoms[atom])

# print(atoms)
# print(abundance)
# print(atoms_verbose)

def cartesian_product(args, repeat=1):
    # cartesian_product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
    # cartesian_product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
    t1 = time.time()
    pools = [tuple(pool) for pool in args] * repeat
    result = [[]]
    for pool in pools:
        result = [x+[y] for x in result for y in pool]
    for prod in result:
        yield tuple(prod)
    t2 = time.time()
    print(f'-> Cartesian Function took {round(t2 - t1, 2)} s')

cartesian_list = [np.asarray(list(pt[atoms_verbose[i]].values()))/100 for i in range(len(atoms_verbose))]
# print(cartesian_list)

matrix = np.asarray(list(cartesian_product(cartesian_list)))
peaks = np.sum(np.asarray(list(cartesian_product([list(pt[atoms_verbose[i]].keys()) for i in range(len(atoms_verbose))]))), axis=1)
# print(matrix)
# print(peaks)

product_vector = np.asarray(np.prod(matrix, axis=1))
# print(product_vector)

set_peaks = np.sort(np.asarray(list(set(peaks))))
# print(set_peaks)

intensity_peaks = []
for i in range(len(set_peaks)):
    a = 0
    for n in range(len(peaks)):
        if set_peaks[i] == peaks[n]:
            a += product_vector[n]
    intensity_peaks.append(a)
intensity_peaks = np.asarray(intensity_peaks)
intensity_peaks *= 100/np.max(intensity_peaks)

intensity_peaks = list(intensity_peaks)
set_peaks = list(set_peaks)
# print(intensity_peaks)

if not nanopeaks:
    old_int_peaks = intensity_peaks[:]
    old_set_peaks = set_peaks[:]
    for p in range(len(old_int_peaks)):
        if old_int_peaks[p] < 0.05:
            intensity_peaks.remove(old_int_peaks[p])
            set_peaks.remove(old_set_peaks[p])

t3 = time.time()
print(f'\n-> Whole program took {round(t3 - t0, 2)} s')
print('------------------------------------------\n')
print('-> Peaks List:\n\nm/z    Rel. Ab.\n***************')
for i in range(len(set_peaks)):
    if intensity_peaks[i] > 0.05:
        print('%-6s %-5s' % (set_peaks[i], round(intensity_peaks[i], 2)))


if plot:
    fig = plt.figure()
    plot = plt.vlines(set_peaks, [0 for _ in range(len(set_peaks))], intensity_peaks)
    plt.xlabel('m/z')
    plt.ylabel('% Abundance')
    plt.title(formula)
    plt.xlim(min(set_peaks) - 1, max(set_peaks) + 1)
    plt.show()

