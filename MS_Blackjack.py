from bs4 import BeautifulSoup
from periodic_table import pt, pt_exact
import numpy as np
import requests
import time
import sys
import re

def pprint(_list):
    for l in _list:
        print(l)

def guide(colors=False):
    if colors:
        b = '\033[1m'  # BOLD code
        e = '\033[0m'  # END code
        y = '\033[93m' # YELLOW code
    else:
        b = ''
        e = ''
        y = ''
    return f'''**************************************************************************************************************
    {y}{b}Mass Fragments Blackjack{e} by Nicolò Tampellini and Alessandro Brusa.

    Insert mass, optional atoms present in the fragment and optional keyword for the search algorithm.
    Example - \'59 CNO2 mol extensive\'

    {b}Keywords: (if multiple, exactly in THIS order){e}
    {y}exact{e} - request exact mass peak (single isotope, best match reported with associate error)
    {y}mol{e} - requests only structures with integer saturation index (molecules or OE(.+) ions, not
        EE(+) ions) and with a C/H ratio compatible with a proper molecule. Be careful, may discard
        low-hydrogen molecules like HC3N or C4N2, which may be proper for human standards.
    {y}extensive{e} - broadens search criterias (required for low %C(m/m) molecules like CS2 and Freons)
    {y}wild{e} - explores all possible chemical space, no assumptions made. (SLOW, USE WITH CARE!)
    {y}only_*{e} - use only selected atoms to build structure. Example - \'only_CHClBr\'
    {y}nist{e} - for each found structure, look up on NIST if there is any matching mass spectra
    **************************************************************************************************************\n'''

def input_alg(colors=False):
    if colors:
        b = '\033[1m'  # BOLD code
        e = '\033[0m'  # END code
        y = '\033[93m' # YELLOW code
    else:
        b = ''
        e = ''
        y = ''
    inp = ''
    while len(inp) == 0:
        try: 
            inp = input(f'{y}m/z{e} [atoms] [keyword(s)] : ')
        except KeyboardInterrupt: 
            sys.exit()
    return inp

def blackjack_alg(inp):

    global everything_is_ok
    everything_is_ok = True

    global even_saturation_index
    even_saturation_index = False

    global exact_mass
    exact_mass = False

    global extensive
    extensive = False

    global wild
    wild = False

    global only_option
    only_option = False

    global nist
    nist = False

    global guided
    guided = False

    global print_list
    print_list = []

    global best_gap
    best_gap = 1e6

    global best_fit
    best_fit = []

    global pt

    def _weight(a):
        return pt[a]

    def _verbose_weight(atoms_verbose):
        return [pt[a] for a in atoms_verbose]

    def _CH_list_generator(mass):
        global everything_is_ok

        c_max = int(np.ceil((mass*1.2-2)/14)) # 1.2 coefficient is ok for as many insaturations as coronene, modeled on (C2H)n fragment
        z = 2                                 #rises to 6 for damned freons and CS2, extensive keyword

        if extensive:
            z = 6

        c_min = max(0, int(np.floor(((mass*2)/25 - 1) / z))) # higher z means lower carbon atoms. Freon require z = 6
        
        if guided and atoms_verbose.count('C') > c_min:  #if user enters more C than minimum, keep those as minimum
            c_min = atoms_verbose.count('C')

        if wild:
            c_min = 0
            c_max = mass // 12
        
        h_min = 0
        h_max = 2 * c_max + 2

        if wild:
            h_max = mass // 1

        if guided and atoms_verbose.count('H') > h_min:
            h_min = atoms_verbose.count('H')

        if guided:
            starting_mass = sum([pt[i] for i in atoms_verbose])
        else:
            starting_mass = 0

        if only_option:
            if ('C' not in only_guide) and ('H' in only_guide):
                h_min = 0
                h_max = mass // 1
                CH_list = []
                for h in range(h_min, h_max + 1):
                    l = []
                    for _ in range(h):
                        l.append('H')
                    CH_list.append(l)
            elif ('H' not in only_guide) and ('C' in only_guide):
                c_min = 0
                c_max = mass // 1
                CH_list = []
                for c in range(c_min, c_max + 1):
                    l = []
                    for _ in range(c):
                        l.append('C')
                    CH_list.append(l)
            elif ('H' not in only_guide) and ('C' not in only_guide):
                return []
        if (not only_option) or ('C' in only_guide and 'H' in only_guide):
            CH_list = []
            for c in range(c_min, c_max + 1): # +1's are to account for range() not comprising last number
                for h in range(h_min, h_max + 1):
                    l = []
                    for _ in range(c):
                        l.append('C')
                    for _ in range(h):
                        l.append('H')
                    CH_list.append(l)
            # for _ in CH_list:
            #     print(_)
            if c_min > c_max:
                print_list.append('\n--> Whoops, those were a few too many Carbon atoms for such a little mass, kid. Try again.\n')
                everything_is_ok = False
            else:
                if extensive:
                    print_list.append(f'--> Extensive sampling : guessing a number of carbon atoms between {c_min} and {c_max}')
                elif not wild:
                    print_list.append(f'--> Guessing a number of carbon atoms between {c_min} and {c_max} (override with "extensive" keyword)')
                if c_max - c_min > 5 and mass - starting_mass > 150:
                    print_list.append('--> Heavy job requested! This is probably going to take a few moments. Have a coffee.')
        return CH_list

    def _mass_of(atom_list):
        mass = 0
        for atom in atom_list:
            mass += pt[atom]
        return mass

    def _pt_trim(pt, overhead):
        pt_local = {}
        pt_local.update(pt)
        for i in pt.keys():
            if pt[i] > overhead:
                pt_local.pop(i)
        return pt_local

    def _flatten(l):
        for i in l: 
            if type(i) == list: 
                _flatten(i)
            else: 
                flattened.append(i) 
        return flattened

    def _cartesian_product(args, repeat=1):
        # cartesian_product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
        # cartesian_product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
        # t1 = time.time()
        pools = [tuple(pool) for pool in args] * repeat
        result = [[]]
        for pool in pools:
            result = [x+[y] for x in result for y in pool]
        for prod in result:
            yield tuple(prod)
        # t2 = time.time()
        # print(f'-> Cartesian Function took {round(t2 - t1, 2)} s')

    
    def _evaluate(CH_guess, candidate, mass):
        global best_gap
        global best_fit
        candidate = list(candidate)
        guess = CH_guess + candidate
        gap = np.abs(np.sum([pt[i] for i in guess]) - mass)
        if gap < best_gap:
            best_gap = gap
            best_fit = guess[:]
            # print('New record, ', best_fit)
        if gap < 1e-1:
            # print('found match ', ''.join(guess))
            matches.append(guess)

    def _saturation_index(atoms_verbose):
        # does not work with P/S(?) - they can be either hypervalent or not
        C = atoms_verbose.count('C')
        H = atoms_verbose.count('H')
        N = atoms_verbose.count('N')
        O = atoms_verbose.count('O')
        F = atoms_verbose.count('F')
        Cl = atoms_verbose.count('Cl')
        Br = atoms_verbose.count('Br')
        I = atoms_verbose.count('I')
        SI = C - (H + F + Cl + Br + I)/2 + N/2 + 1
        proper = True if C + N + O != 0 else False
        if proper:
            return SI
        return 0.5

    def _CH_test(atoms_verb):
        global extensive
        C = atoms_verb.count('C')
        H = atoms_verb.count('H')
        if extensive:
            return True
        if 2 * C < H - 2:             # checks if 'too saturated'
            return False
        elif C > 2 * H:               # checks if 'too unsaturated'
            return False
        else:
            return True

    def _nist_structures(l):
        at_set = set(l)
        formula = ''
        for i in at_set:
            count = l.count(i)
            formula += i
            if count > 1:
                formula += str(count)
        url = f'https://webbook.nist.gov/cgi/cbook.cgi?Formula={formula}&NoIon=on&Units=SI&cMS=on'
        data = requests.get(url).text
        soup = BeautifulSoup(data, 'lxml').body.main.p
        if len(str(soup).split('\n')) == 1:        # if url p1 is a p class there's just one species with that formula
            return '1', url
        num = str(soup).split('\n')[1].split()[0]
        try:
            num = int(num)
        except:
            pass
        if num == 'No':
            return '0', url
        if num == 'NIST':
            return '1', url
        if type(num) is int:
            return str(num), url
        print('------------>', num)

    inp = inp.split()
    try:
        if inp[-1] == 'nist':
            print_list.append('--> NIST keyword - Will look online for matching MS spectras.')
            nist = True
            inp = inp[:-1]

        if 'only_' in inp[-1]:
            only_option = True
            only_guide = inp[-1].split('_')[-1]
            inp = inp[:-1]

        if inp[-1] == 'wild':
            print_list.append('--> Running wild! Requested comprehensive sampling of chemical space.')
            wild = True
            inp = inp[:-1]

        if inp[-1] == 'extensive':
            print_list.append(r'--> Requested extensive sampling of molecules with low %m/m of C.')
            extensive = True
            inp = inp[:-1]

        if inp[-1] == 'mol':
            print_list.append('--> Requested only molecules with even saturation index and common C/H ratio.')
            even_saturation_index = True
            inp = inp[:-1]

        if inp[-1] == 'exact':
            print_list.append('--> Requested exact mass peak - only best match provided.')
            exact_mass = True
            inp = inp[:-1]

        if not exact_mass:
            inp[0] = int(inp[0])
        else:
            inp[0] = float(inp[0])
    except:
        pass
    
    mass = inp[0]
    if len(inp) > 1:  # if len(inp) > 1 it must be guided, so inp[-1] is the guide
        guided = True
        guide = inp[-1]

        parsed = re.findall('[A-Z][^A-Z]*', guide)         #parser for chemical formula interpretation
        parsed = [re.split('(\d+)',i) for i in parsed]
        flattened = []
        parsed = _flatten(parsed)
        for i in parsed:
            if i == '':
                parsed.remove(i)
        # print_list.append('--------> NOT Parsed:', parsed)
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
        # print_list.append('--------> Parsed:', parsed)
        # exit()
        
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
    t_start = time.time()

    if exact_mass:
        pt.update(pt_exact)
    if only_option:
        only_guide = re.findall('[A-Z][^A-Z]*', only_guide)
        # print_list.append(only_guide)
        if all([a in pt for a in only_guide]):
            pt_only = {}
            for element in only_guide:
                pt_only[element] = pt[element]
            pt = pt_only.copy()
            only_elements_string = ''
            for e in range(len(only_guide)):
                only_elements_string += only_guide[e]
                if e == len(only_guide) - 2:
                    only_elements_string += ' and '
                elif e == len(only_guide) - 1:
                    pass
                else:
                    only_elements_string += ', '
            print_list.append(f'--> Using only atoms {only_elements_string} to build structure.')
        else:
            # print_list.append(only_guide)
            unknowns_indexes = np.nonzero([a not in pt for a in only_guide])[0]
            unknowns = [only_guide[int(u)] for u in unknowns_indexes]
            for u in range(len(unknowns) - 1):
                unknowns[u] = str(unknowns[u]) + ' '
            unknowns = ''.join(unknowns)
            print_list.append(f'--> Sorry, what do you mean with {unknowns}? Try again if you misspelled an atom or maybe teach me more chemistry!')
            everything_is_ok = False

    pt_no_CH = {}
    pt_no_CH.update(pt)
    if 'H' in pt_no_CH.keys():
        pt_no_CH.pop('H')
    if 'C' in pt_no_CH.keys():
        pt_no_CH.pop('C')
    if len(pt_no_CH.keys()) > 0:
        lightest_non_CH = min(pt_no_CH.values())
    else:
        lightest_non_CH = 1e6


    if guided:
        if not all([a in pt for a in atoms]):
            guided_unknowns_indexes = np.nonzero([a not in pt for a in atoms_verbose])[0]
            guided_unknowns = [atoms_verbose[int(u)] for u in guided_unknowns_indexes]
            for u in range(len(guided_unknowns) - 1):
                guided_unknowns[u] = str(guided_unknowns[u]) + ' '
            guided_unknowns = ''.join(guided_unknowns)
            print_list.append(f'--> Sorry, what do you mean with {guided_unknowns}? Try again if you misspelled an atom or maybe teach me more chemistry!')
            everything_is_ok = False

        if everything_is_ok and sum([pt[i] for i in atoms_verbose]) > mass:
            print_list.append(f'\n--> Sorry, but the guess you provided ({guide}) weights more than {mass} amu.\n')
            everything_is_ok = False

    if everything_is_ok:
        if guided:
            print_list.append(f'--> Using fragment {guide} as starting guess.')
        # print_list.append(f'Computing ions of {mass} amu...', end='\r')
        CH_list = _CH_list_generator(mass)
        hetero_list = []
        hetero_oh = []
        matches = []
        log = []
        if len(CH_list) > 0:
            for l in CH_list:
                # print_list.append(''.join(l))
                overhead = mass - _mass_of(l)
                if overhead >= lightest_non_CH:
                    hetero_list.append(l)
                    hetero_oh.append(overhead)
                elif 0 <= overhead < 1e-1:
                    matches.append(l)
                    if overhead < best_gap:
                        best_gap = overhead
                        best_fit = l[:]
        else:                                    # if no C nor H, overhead is the whole mass
            hetero_list.append([])
            hetero_oh.append(mass)

        for l in range(len(hetero_list)):
            for number_of_heteroatoms in range(1, 100):
                if lightest_non_CH*number_of_heteroatoms <= hetero_oh[l]:
                    pt_local = _pt_trim(pt_no_CH, hetero_oh[l]) #creates a periodic table without the elements that are too heavy
                    pt_list = [pt_local.keys() for _ in range(number_of_heteroatoms)]
                    # print_list.append(number_of_heteroatoms, ' het between', pt_local.keys())
                    candidates = _cartesian_product(pt_list)
                    for c in candidates:
                        # print_list.append(''.join(hetero_list[l]), c)
                        _evaluate(hetero_list[l], c, mass)
                else:
                    # print_list.append('# Max het were', number_of_heteroatoms - 1)
                    break

        if len(matches) > 0:
            matches_sorted = [sorted(m, key=_weight) for m in matches]
            for m in matches_sorted:
                if even_saturation_index:
                        if not (_saturation_index(m) % 1 < 1e-5 and _CH_test(m)):
                            continue
                if guided:
                    if not all([m.count(i) >= atoms_verbose.count(i) for i in atoms_verbose]):
                        continue
                if m not in log:
                    log.append(m)

            if exact_mass and len(matches) > 0:
                log = []
                log.append(best_fit)

        # print_list.append('\033[K', end='\r')
        # print_list.append('                                                   \n', end='\r')
        if len(log) != 0:
            log.sort(key=_verbose_weight)
            for l in log:
                formula = []
                for letter in list(dict.fromkeys(l)):
                    c = l.count(letter)
                    formula.append(letter)
                    if c > 1:
                        formula.append(str(c))
                    formula.append(' ')
                if nist:
                    NOS = _nist_structures(l)
                    if NOS[0] != '0':
                        string = ' ' * (15 - len(''.join(formula))) + NOS[0] + ' structures found on NIST - ' + NOS[1]
                        formula += string
                print_list.append(''.join(formula))   
            t_end = time.time()
            if exact_mass:
                print_list.append(f'\n--> Best formula found, error {round(best_gap, 6)} amu - Took {round(t_end - t_start, 2)} s')
            else:
                s = 's' if len(log) > 1 else ''
                print_list.append(f'\n--> Found {len(log)} formula{s} - Took {round(t_end - t_start, 2)} s') 
        else:
            print_list.append('Sorry, no combination matches your criteria.')
    print_list.append('**************************************************************************************************************\n')
    return print_list

if __name__ == "__main__":
    print(guide())
    while True:
        input_string = input_alg()
        output = blackjack_alg(input_string)
        pprint(output)