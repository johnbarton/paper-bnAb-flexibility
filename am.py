# Affinity maturation simulation in Python

import sys
import numpy as np                          # numerical tools
from copy import deepcopy                   # deepcopy copies a data structure without any implicit references
from timeit import default_timer as timer   # timer for performance


###### Global parameters ######
    
p_mut        = 0.20                             # probability of mutation per division round
p_CDR        = 0.85                             # probability of mutation in the CDR region
p_CDR_lethal = 0.30                             # probability that a CDR mutation is lethal
p_CDR_silent = 0.50                             # probability that a CDR mutation is silent
p_CDR_affect = 1. - p_CDR_lethal - p_CDR_silent # probability that a CDR mutation affects affinity
p_var        = 0.10                             # probability that a CDR mutation affects the variable region
p_cons       = 1.0 - p_var                      # probability that a CDR mutation affects the conserved (constant) region
p_FR_lethal  = 0.80                             # probability that a framework (FR) mutation is lethal
p_FR_silent  = 0.                               # probability that a FR mutation is silent
p_FR_affect  = 1. - p_FR_lethal - p_FR_silent   # probability that a FR mutation affects flexibility

Ec_start     = 0               # starting value for Ec
nb_Ag        = 10              # number of antigens
conc         = 1.18            # antigen concentration
energy_scale = 0.07            # inverse temperature
E0           = 4.00            # mean binding energy for mixing with flexibility
sigmaQ       = 0.08            # standard deviation for changes in flexibility with FR mutation
help_cutoff  = 0.70            # only B cells in the top (help_cutoff) fraction of binders receive T cell help
p_recycle    = 0.70            # probability that a B cell is recycled
p_exit       = 1. - p_recycle  # probability that a B cell exits the GC

mu     =  1.9   # lognormal mean
sigma  =  0.5   # lognormal standard deviation
corr   =  0.0   # correlation between antigen variable regions
o      =  3.0   # lognormal offset
mumat  = mu * np.ones(nb_Ag)
sigmat = sigma * np.diag(np.ones(nb_Ag))
for i in range(nb_Ag):
    for j in range(i+1,nb_Ag):
        sigmat[i,j] = sigma * corr
        sigmat[j,i] = sigma * corr


###### B cell clone class ######

class BCell:

    def __init__(self, nb = 512, **kwargs):
        """ Initialize clone-specific variables. 
            nb          - population size of the clone
            Ev          - binding energy for each Ag variable region
            Ec          - binding energy for the constant region
            Q           - rigidity parameter (most 0 ---> 1 least flexible)
            nb_FR_mut   - number of accumulated FR mutations
            nb_CDR_mut  - number of accumulated CDR mutations
            last_bound  - number of individual B cells in the clone that last bound each Ag
            generation  - generation in the GC reaction
            history     - history of mutations, generation occurred, and effect on Q/Ec """
        
        self.nb = nb    # default starting population size = 512 (9 rounds of division)
        
        if ('Ev' in kwargs) and ('Ec' in kwargs):
            self.Ev = np.array(kwargs['Ev'])
            self.Ec = kwargs['Ec']
        
        else:
            self.Ev = np.zeros(nb_Ag)
            self.Ec = Ec_start

            self.Ev     = o - np.exp(np.random.multivariate_normal(mumat, sigmat))
            selected_Ag = np.argmax(self.Ev)
            
            if self.Ev[selected_Ag]<0 and np.max(self.Ev)<0:
                self.Ev[selected_Ag] = 0
            
        if 'Q' in kwargs: self.Q = kwargs['Q']
        else:             self.Q = 0.8

        if 'nb_FR_mut' in kwargs: self.nb_FR_mut = kwargs['nb_FR_mut']
        else:                  self.nb_FR_mut = 0
        
        if 'nb_CDR_mut' in kwargs: self.nb_CDR_mut = kwargs['nb_CDR_mut']
        else:                  self.nb_CDR_mut = 0
        
        if 'last_bound' in kwargs: self.last_bound = kwargs['last_bound']
        else:                      self.last_bound = np.random.multinomial(self.nb, pvals = [1/float(nb_Ag)] * nb_Ag)
        
        if 'generation' in kwargs: self.generation = kwargs['generation']
        else:                  self.generation = 0
           
        if 'history' in kwargs: self.history = kwargs['history']
        else:                   self.history = {'generation' : [self.generation], 'nb_FR_mut' : [self.nb_FR_mut], 'nb_CDR_mut' : [self.nb_CDR_mut], 'Q' : [self.Q], 'Ec' : [self.Ec]}

    """ Return a new copy of the input BCell"""
    @classmethod
    def clone(cls, b):
        return cls(1, Ev = deepcopy(b.Ev), Ec = b.Ec, Q = b.Q, generation = b.generation, nb_FR_mut = b.nb_FR_mut, nb_CDR_mut = b.nb_CDR_mut, last_bound = deepcopy(b.last_bound), history = deepcopy(b.history))
    
    def update_history(self):
        """ Add current parameters to the history list. """
        self.history['Q'].append(self.Q)
        self.history['Ec'].append(self.Ec)
        self.history['generation'].append(self.generation)
        self.history['nb_FR_mut'].append(self.nb_FR_mut)
        self.history['nb_CDR_mut'].append(self.nb_CDR_mut)
    
    def bind_to(self, Ag):
        """ Return binding energy with input antigen. """
        return (self.Q * (self.Ec + (p_var * (self.Ev[Ag] - self.Ec)))) + ((1 - self.Q) * E0)

    def divide(self):
        """ Run one round of division. """
        self.nb *= 2
        self.generation += 1
    
    def mutate_CDR(self):
        """ Change in energy due to affinity-affecting CDR mutation. """
        if np.random.rand()<p_var:
            self.Ev += o - np.exp(np.random.multivariate_normal(mumat, sigmat))
        else:
            self.Ec += o - np.exp(np.random.normal(mu, sigma))
        self.nb_CDR_mut += 1
        self.update_history()

    def mutate_FR(self):
        """ Change in flexibility due to affinity-affecting framework (FR) mutation. """
        dQ = np.random.normal(0, sigmaQ)
        if   self.Q + dQ > 1.0:
            self.Q   = 1.0
        elif self.Q + dQ < 0.1:
            self.Q   = 0.1
        else:
            self.Q = self.Q + dQ
        self.nb_FR_mut += 1
        self.update_history()       

    def shm(self):
        """ Run somatic hypermutation and return self + new B cell clones. """
        
        # get number of cells that mutate
        new_clones = []
        n_mut      = np.random.binomial(self.nb, p_mut)
        self.nb   -= n_mut
            
        # get number of CDR vs framework (FR) mutations
        n_CDR = np.random.binomial(n_mut, p_CDR)
        n_FR  = n_mut - n_CDR
            
        # process CDR mutations
        n_die, n_silent, n_affect  = np.random.multinomial(n_CDR, pvals = [p_CDR_lethal, p_CDR_silent, p_CDR_affect])
        self.nb                   += n_silent
        for i in range(n_affect):
            b = BCell.clone(self)
            b.mutate_CDR()
            new_clones.append(b)
        
        # process FR mutations
        n_die, n_silent, n_affect  = np.random.multinomial(n_FR, pvals = [p_FR_lethal, p_FR_silent, p_FR_affect])
        self.nb                   += n_silent
        for i in range(n_affect):
            b = BCell.clone(self)
            b.mutate_FR()
            new_clones.append(b)

        # return the result
        if (self.nb>0): new_clones.append(self)
        return new_clones


###### Main functions ######


def usage():
    print("")


def main(verbose=False):
    """ Simulate the affinity maturation process in a single germinal center (GC) and save the results to a CSV file. """
    
    # Run multiple trials and save all data to file
    
    nb_trial = 100
    start    = timer()
    
    fmem = open('output-memory.csv', 'w')
    ftot = open('output-total.csv',  'w')
    fbig = open('output-largest-clone.csv', 'w')
    
    fmem.write('trial,exit cycle,number,generation,FR_mutations,CDR_mutations,mutations,Q,Ec,'+(','.join(['Ev'+str(i) for i in range(nb_Ag)]))+'\n')
    ftot.write('trial,cycle,number recycled,number exit\n')
    fbig.write('trial,exit cycle,Q,Ec,generation,FR_mutations,CDR_mutations,mutations\n')
    
    for t in range(nb_trial):
    
        print_update(t, nb_trial)   # status check

        # INITIALIZATION - DEFINE DATA STRUCTURES

        recycled_cells = []
        exit_cells     = []
        nb_recycled    = []
        nb_exit        = []

        
        # CYCLES 1 + 2 - CREATE FOUNDERS AND REPLICATE WITHOUT MUTATION
        
        nb_founders = 3                                     # number of founder B cells for a GC
        B_cells     = [BCell() for i in range(nb_founders)]
        
        # Update data
        nb_recycled.append(nb_founders)                     # all founders are recycled
        nb_exit.append(0)                                   # no founders exit the GC
        nb_recycled.append(np.sum([b.nb for b in B_cells])) # all founders replicate and are recycled
        nb_exit.append(0)                                   # no founders exit
        

        # AFFINITY MATURATION
        
        GC_size_max  = np.sum([b.nb for b in B_cells])  # maximum number of cells in the GC (= initial population size)
        nb_cycle_max = 250                              # maximum number of GC cycles
        cycle_number = 2
        
        for cycle_number in range(2, nb_cycle_max):
        
            B_cells, out_cells = run_GC_cycle(B_cells)
            GC_size            = np.sum([b.nb for b in B_cells])       # total number of cells in the GC
            
            if (cycle_number==nb_cycle_max-1) or (GC_size>GC_size_max): # at the end, all B cells exit the GC
                out_cells += B_cells
            else: out_cells = []                                        # avoid recording exit cells except at final step
            
            recycled_cells.append([deepcopy(b) for b in B_cells])
            exit_cells.append(out_cells)
            nb_recycled.append(np.sum([b.nb for b in B_cells]  ))
            nb_exit.append(    np.sum([b.nb for b in out_cells]))

            if (nb_recycled[-1]==0) or (GC_size>GC_size_max): break
        

        # SAVE OUTPUT
        
        for i in range(len(exit_cells)):
            for b in exit_cells[i]:
                if b.nb>100:
                    fmem.write('%d,%d,%d,%d,%d,%d,%d,%lf,%lf' % (t, i+2, b.nb, b.generation, b.nb_FR_mut, b.nb_CDR_mut, b.nb_FR_mut+b.nb_CDR_mut, b.Q, b.Ec))
                    for j in range(nb_Ag): fmem.write(',%lf' % b.Ev[j])
                    fmem.write('\n')
        fmem.flush()

        for i in range(len(nb_recycled)): ftot.write('%d,%d,%d,%d\n' % (t, i+1, nb_recycled[i],nb_exit[i]))
        ftot.flush()
        
        if len(exit_cells[-1])>0:
            idx = np.argmax([b.nb for b in exit_cells[-1]])
            b   = exit_cells[-1][idx]
            for i in range(len(b.history['Q'])):
                fbig.write('%d,%d,%lf,%lf,%d,%d,%d,%d\n' % (t, len(exit_cells)-1, b.history['Q'][i], b.history['Ec'][i], b.history['generation'][i], b.history['nb_FR_mut'][i], b.history['nb_CDR_mut'][i], b.history['nb_FR_mut'][i]+b.history['nb_CDR_mut'][i]))
        fbig.flush()

    # End and output total time
    
    fmem.close()
    ftot.close()
    fbig.close()
    
    end = timer()
    print('\nTotal time: %lfs, average per cycle %lfs' % ((end - start),(end - start)/float(nb_trial)))


def print_update(current, end, bar_length=20):
    """ Print an update of the simulation status. h/t Aravind Voggu on StackOverflow. """
    
    percent = float(current) / end
    dash    = ''.join(['-' for k in range(int(round(percent * bar_length)-1))]) + '>'
    space   = ''.join([' ' for k in range(bar_length - len(dash))])

    sys.stdout.write("\rSimulating: [{0}] {1}%".format(dash + space, int(round(percent * 100))))
    sys.stdout.flush()


def run_dark_zone(B_cells, nb_rounds = 2):
    """ B cells proliferate and undergo SHM in the dark zone. """
    
    for i in range(nb_rounds):
        new_cells = []
        for b in B_cells:
            b.divide()
            new_cells += b.shm()
        B_cells = new_cells
    return B_cells


def run_binding_selection(B_cells):
    """ Select B cells for binding to antigen. """
    
    for b in B_cells:
        
        b.last_bound = np.random.multinomial(b.nb, pvals = [1./float(nb_Ag)] * nb_Ag)
        
        for i in range(nb_Ag):
            
            # compute binding energy and chance of death ( = 1 - chance of survival )
            Ag_bound      = np.exp(energy_scale * b.bind_to(i))
            factor        = conc * Ag_bound
            langmuir_conj = 1. / (1. + factor)
            
            # remove dead cells and update binding details
            n_die            = np.random.binomial(b.last_bound[i], langmuir_conj)
            b.nb            -= n_die
            b.last_bound[i] -= n_die


def run_help_selection(B_cells):
    """ Select B cells to receive T cell help. """
    
    # get binding energies
    binding_energy     = [[b.bind_to(i) for i in range(nb_Ag)] for b in B_cells]
    binding_energy_tot = []
    for i in range(len(B_cells)):
        for j in range(nb_Ag): binding_energy_tot += [binding_energy[i][j]] * B_cells[i].last_bound[j]
    
    # cells in the top (help_cutoff) fraction of binders survive
    if len(binding_energy_tot)>0:
        cut_idx       = np.max([0, int(np.floor(help_cutoff * len(binding_energy_tot)))-1])
        energy_cutoff = np.array(binding_energy_tot)[np.argsort(binding_energy_tot)][::-1][cut_idx]
        n_die_tie     = len(binding_energy_tot) - cut_idx - np.sum(binding_energy_tot < energy_cutoff)

        # kill all B cells below threshold
        for i in np.random.permutation(len(B_cells)):
            for j in np.random.permutation(nb_Ag):
                energy = binding_energy[i][j]
                if energy < energy_cutoff:
                    B_cells[i].nb            -= B_cells[i].last_bound[j]
                    B_cells[i].last_bound[j]  = 0
                elif (energy == energy_cutoff) and (n_die_tie > 0):
                    if B_cells[i].last_bound[j] < n_die_tie:
                        B_cells[i].nb            -= B_cells[i].last_bound[j]
                        n_die_tie                -= B_cells[i].last_bound[j]
                        B_cells[i].last_bound[j]  = 0
                    else:
                        B_cells[i].nb            -= n_die_tie
                        B_cells[i].last_bound[j] -= n_die_tie
                        n_die_tie                 = 0


def run_recycle(B_cells):
    """ Randomly select B cells to be recycled back into the GC or to exit. """

    new_cells  = []                                 # cells that will remain in the GC
    exit_cells = []                                 # cells that will exit the GC
    n_tot      = np.sum([b.nb for b in B_cells])    # total number of cells currently in the GC
    n_exit     = int(np.floor(p_exit * n_tot))      # number of cells that will exit the GC
    b_exit     = []                                 # index of cells that will exit the GC

    if (n_tot > 0) and (n_exit > 0):
        b_exit = np.random.choice(n_tot, n_exit)
    b_exit = np.array(b_exit)

    idx = 0
    for b in B_cells:
    
        # find which cells exit the GC
        n_exit  = np.sum((idx <= b_exit) * (b_exit < idx + b.nb))
        idx    += b.nb
        b.nb   -= n_exit
        
        # add remainder to recycled cells
        if (b.nb>0):
            new_cells.append(b)
    
        # record exit cells
        if (n_exit>0):
            exit_cells.append(deepcopy(b))
            exit_cells[-1].nb = n_exit

    return new_cells, exit_cells

    # STOCHASTIC RECYCLING
#    for b in B_cells:
#        
#        n_exit  = np.random.binomial(b.nb, p_exit)
#        b.nb   -= n_exit
#        
#        if (b.nb>0):
#            new_cells.append(b)
#        
#        if (n_exit>0):
#            exit_cells.append(deepcopy(b))
#            exit_cells[-1].nb = n_exit
#
#    return new_cells, exit_cells


def run_GC_cycle(B_cells):
    """ Run one cycle of the GC reaction. """

    B_cells = run_dark_zone(B_cells)    # DARK  ZONE - two rounds of division + SHM
    run_binding_selection(B_cells)      # LIGHT ZONE - selection for binding to Ag
    run_help_selection(B_cells)         # LIGHT ZONE - selection to receive T cell help
    return run_recycle(B_cells)         # RECYCLE    - randomly pick exiting cells from the surviving B cells


################################################

def test_dark_zone():
    """ Test the run_dark_zone function. """
    
    print('test_dark_zone...')
    if p_CDR<1:
        print('Rewrite this test without assuming p_CDR=1!')
        return 0
    
    n_clones    = []
    n_cells_max = []
    n_tests     = 10000
    n_cells     = 1000
    
    for i in range(n_tests):
        test_cells = [BCell(nb = n_cells)]
        test_cells = run_dark_zone(test_cells, nb_rounds=1)
        
        temp_clones = 0
        temp_max    = 0
        for b in test_cells:
            if b.nb>0:        temp_clones += 1
            if b.nb>temp_max: temp_max     = b.nb
        n_clones.append(temp_clones)
        n_cells_max.append(temp_max)

    E_cells_max     = 2. * n_cells * (1. - p_mut + (p_mut * p_CDR_silent))           # = np
    E_std_cells_max = np.sqrt( E_cells_max * (1. - (E_cells_max / (2. * n_cells))) ) # = sqrt(np(1-p))

    print('E(max cells):     %.3e\tgot: %.3e' % (E_cells_max,    np.mean(n_cells_max)))
    print('E(std max cells): %.3e\tgot: %.3e' % (E_std_cells_max, np.std(n_cells_max)))

    E_clones     = 2. * n_cells * p_mut * p_CDR_affect
    E_std_clones = np.sqrt( E_clones * (1. - (E_clones / (2. * n_cells))))

    print('E(clones):        %.3e\tgot: %.3e' % (E_clones,    np.mean(n_clones)))
    print('E(std clones):    %.3e\tgot: %.3e' % (E_std_clones, np.std(n_clones)))
    print('')


def test_binding_selection():
    """ Test the run_binding_selection function. """
    
    print('test_binding_selection...')
    if nb_Ag!=2:
        print('Rewrite this function for nb_Ag!=2!')
        return 0
    
    n_bound_d = []
    n_survive = []
    n_tests   = 100000
    test_E    = np.array([ 0.,  1.,  2.,  4.])
    n_cells   = np.array([250, 250, 250, 250])
    
    for t in range(n_tests):
        test_cells = [BCell(nb = n_cells[i], Ec = test_E[i], Ev = [0 for k in range(nb_Ag)]) for i in range(len(test_E))]
        run_binding_selection(test_cells)

        n_bound_d.append([b.last_bound[0]-b.last_bound[1] for b in test_cells])
        n_survive.append([b.nb                            for b in test_cells])
        for b in test_cells:
            if b.nb!=np.sum(b.last_bound):
                print('b.nb!=sum(b.last_bound)!')
                return 1

    n_bound_d = np.array(n_bound_d)
    n_survive = np.array(n_survive)

    E_bound_d = [0., 0., 0., 0.]

    print('E(bound difference): '+str(E_bound_d)+'\tgot: '+str(np.mean(n_bound_d, axis=0)))

    E_survive     = n_cells * (conc * np.exp(energy_scale * test_E) / (1. + (conc * np.exp(energy_scale * test_E))))
    E_std_survive = np.sqrt( E_survive * (1. - (E_survive / np.array(n_cells, float))))

    print('E(survive):     '+str(E_survive)    +'\tgot: '+str(np.mean(n_survive, axis=0)))
    print('E(std survive): '+str(E_std_survive)+'\tgot: '+str( np.std(n_survive, axis=0)))
    print('')


def test_help_selection():
    """ Test the run_help_selection function. """
    
    print('test_help_selection...')
    if help_cutoff!=0.7:
        print('Rewrite this for a different value of help_cutoff!')
    
    test_E     = np.array([ 0.,  1.,  2.,  4.])
    n_cells    = np.array([250, 250, 250, 250])
    test_cells = [BCell(nb = n_cells[i], Ec = test_E[i], Ev = [0 for k in range(nb_Ag)]) for i in range(len(test_E))]
    run_help_selection(test_cells)

    n_survive = [b.nb for b in test_cells]
    E_survive = [0, 200, 250, 250]

    print('E(survive): '+str(E_survive)+'  \tgot: '+str(n_survive))

    test_E     = np.array([ 4.,  1.,  2.,  0.])
    n_cells    = np.array([250, 250, 250, 250])
    test_cells = [BCell(nb = n_cells[i], Ec = test_E[i], Ev = [0 for k in range(nb_Ag)]) for i in range(len(test_E))]
    run_help_selection(test_cells)

    n_survive = [b.nb for b in test_cells]
    E_survive = [250, 200, 250, 0]

    print('E(survive): '+str(E_survive)+'  \tgot: '+str(n_survive))

    test_E    = np.array([ 0.,  0.,  0.,  0.])
    n_cells   = np.array([250, 250, 250, 250])
    E_survive = [175, 175, 175, 175]
    n_survive = []
    n_tests   = 10000
    for t in range(n_tests):
        test_cells = [BCell(nb = n_cells[i], Ec = test_E[i], Ev = [0 for k in range(nb_Ag)]) for i in range(len(test_E))]
        run_help_selection(test_cells)
        n_survive.append([b.nb for b in test_cells])

    print('E(survive): '+str(E_survive)+'\tgot: '+str(np.mean(n_survive, axis=0)))
    print('')


def test_recycle():
    """ Test the run_recycle function. """
    
    print('test_recycle...')
    if p_recycle!=0.7:
        print('Rewrite this for a different value of p_recycle!')

    n_cells   = np.array([250, 250, 250, 250])
    E_recycle = [175, 175, 175, 175]
    n_recycle = []
    E_exit    = [75, 75, 75, 75]
    n_exit    = []
    n_tests   = 10000
    for t in range(n_tests):
        test_cells = [BCell(nb = n) for n in n_cells]
        new, exit  = run_recycle(test_cells)
        
        n_recycle.append([b.nb for b in new])
        n_exit.append(   [b.nb for b in exit])

    print('E(recycle): '+str(E_recycle)+'\tgot: '+str(np.mean(n_recycle, axis=0)))
    print('E(exit):    '+str(E_exit)   +'\tgot: '+str(np.mean(n_exit,    axis=0)))
    print('')


def run_tests():
    """ Run diagnostic tests to make sure that the code is functioning as expected. """

    print('running tests\n')
    test_dark_zone()
    test_binding_selection()
    test_help_selection()
    test_recycle()
    print('done')


if __name__ == '__main__': main()

