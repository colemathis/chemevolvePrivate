import numpy as np
from InitializeFunctions import *
import random
import CoreClasses

def check_mass(original_mass, CRS, concentrations):
	''' Checkes conservation of mass
	Arguements
		- original_mass: integer mass of the original system
		- CRS: CRS object
		- concentrations: array of molecule abundances

		 '''
	mass_conserved = False
	test_mass = 0.0
	molecules = CRS.molecule_list
	for m in molecules:
		index =	molecules.index(m)

		molecule_size = len(m)
		molecule_count = np.sum(concentrations[:,:,index])
		mass = molecule_count*molecule_size

		test_mass += mass 
	if test_mass == original_mass:
		mass_conserved = True
	return mass_conserved, test_mass

def get_composition(molecule_list, molecule_ID):
	'''Returns the 'A' and 'B' composition of a molecule string '''
	molecule = molecule_list[molecule_ID]
	nA = molecule.count('A')
	nB = molecule.count('B')
	return nA, nB

def generate_all_replication2_reactions(length, rconstant = 0.001, sconstant = 1E-5, dconstant = 1.0, **kwargs):
	''' Generates all possible replication reactions for binary molecules of size length assigns all reactions the same constant
	forward reactions are generated twice, once for replication, given RCM propensity function, and Spontaneous given standard propensity
	Arguements:
		- length: length of binary replicators 
		- rconstant: replication rate constant
		- sconstant: Spontaneous formation rate of replicators 
		- dconstant: degradation rate constant for replicators
		Keyword Arguements:
		- rep_landscape: If this is passed the replication rate is based on the assigned landscape, either 'gaussian' or 'flat'
		- rmaster: Master sequence for the replication landscape
		- d_landscape: If this is passed the degradation rate is based on the assigned landscape, either 'gaussian' or 'flat'
		- dmaster: Master sequence for the degradation rate

	Fitness landscapes impose a max difference of 1 order of magnitude (ie 10 times faster replication for the master sequence, 
		or 10 times slower degradation rate for master sequence)
		'''
	import CoreClasses
	import itertools
	rxn_IDs = []
	reaction_list = []
	molecule_list = []
	molecule_dict = {}
	rxn_ID = 0
	#molecule_ID = 0
	monomers = 'AB'
	molecule_list.append('A')
	molecule_dict['A'] = molecule_list.index('A')
	molecule_list.append('B')
	molecule_dict['B'] = molecule_list.index('B')
	# Generate all sequences of length l
	original_rconstant = rconstant
	original_dconstant = dconstant
	sequences = itertools.product(monomers, repeat = length)
	for seq in sequences:
		s = ''.join(seq)
		molecule_list.append(s)
		molecule_dict[s] = molecule_list.index(s)
		if 'rep_landscape' in kwargs:
			master_seq = kwargs['rmaster']
			d = hamming_distance(master_seq, s)
			if kwargs['rep_landscape'] == 'gaussian':
				rconstant = original_rconstant*(1.0 + 10.0*gaussian(d, 0, 2) )
			elif kwargs['rep_landscape'] == 'step':
				if d != 0:
					rconstant = original_rconstant
				else:
					rconstant = 10*original_rconstant
		if 'd_landscape' in kwargs:
			master_seq = kwargs['dmaster']
			d = hamming_distance(master_seq, s)
			if kwargs['d_landscape'] == 'gaussian':
				dconstant = original_dconstant*(1.0 - 0.9*gaussian(d, 0, 2) )
			elif kwargs['d_landscape'] == 'step':
				if d != 0:
					dconstant = original_dconstant
				else:
					dconstant = 0.1*original_dconstant

		# Forward Replication Reaction
		rxn_ID = len(reaction_list)
		reactants =  [molecule_dict['A'], molecule_dict['B'] ] 
		reactant_coeff = [ s.count('A'), s.count('B') ]
		products = [ molecule_dict[s] ]
		product_coeff = [1]
		reaction_list.append( CoreClasses.Reaction(rxn_ID, reactants = reactants, reactant_coeff = reactant_coeff, products = products, product_coeff = product_coeff, constant = rconstant, prop = 'RM2') )
		#print 'Reaction List index: ', rxn_ID, 'Reaction ID: ', reaction_list[rxn_ID].ID 
		rxn_IDs.append(rxn_ID)

		# Spontaneous Reaction
		rxn_ID = len(reaction_list)
		reaction_list.append( CoreClasses.Reaction(rxn_ID, reactants = reactants, reactant_coeff = reactant_coeff, products = products, product_coeff = product_coeff, constant = sconstant, prop = 'STD') )
		rxn_IDs.append(rxn_ID)
		# Degradation Reaction
		rxn_ID = len(reaction_list)
		reaction_list.append( CoreClasses.Reaction(rxn_ID, products = reactants, product_coeff = reactant_coeff, reactants = products, reactant_coeff = product_coeff, constant = dconstant, prop = 'STD') )
		#print 'Reaction List index: ', rxn_ID, 'Reaction ID: ', reaction_list[rxn_ID].ID 
		rxn_IDs.append(rxn_ID)
				
	#print molecule_dict
	newCRS = CoreClasses.CRS(molecule_list = molecule_list, molecule_dict = molecule_dict, reactions = reaction_list)
	#print newCRS.molecule_dict
	fname = ''
	if 'rep_landscape' in kwargs:

		master_seq = kwargs['rmaster']
		fname += 'rmaster_' + str(master_seq)
		if kwargs['rep_landscape'] == 'gaussian':
			fname+= '_gaussian_'	
		elif kwargs['rep_landscape'] == 'step':
			fname+= '_step_'
	if 'd_landscape' in kwargs:
		master_seq = kwargs['dmaster']
		fname += 'dmaster_' + str(master_seq)
		if kwargs['d_landscape'] == 'gaussian':
			fname+= '_gaussian_'
		elif kwargs['d_landscape'] == 'step':
			fname+= '_step_'
	fname += 'length_%i.txt' % length
	newCRS.savetxt(fname)
	print "New CRS saved"

def generate_all_replication1_reactions(length, rconstant = 0.001, sconstant = 1E-5, dconstant = 1.0, **kwargs):
	''' Generates all possible replication reactions for binary molecules of size length assigns all reactions the same constant
	forward reactions are generated twice, once for replication, given RCM propensity function, and Spontaneous given standard propensity
	Arguements:
		- length: length of binary replicators 
		- rconstant: replication rate constant
		- sconstant: Spontaneous formation rate of replicators 
		- dconstant: degradation rate constant for replicators
		Keyword Arguements:
		- rep_landscape: If this is passed the replication rate is based on the assigned landscape, either 'gaussian' or 'flat'
		- rmaster: Master sequence for the replication landscape
		- d_landscape: If this is passed the degradation rate is based on the assigned landscape, either 'gaussian' or 'flat'
		- dmaster: Master sequence for the degradation rate

	Fitness landscapes impose a max difference of 1 order of magnitude (ie 10 times faster replication for the master sequence, 
		or 10 times slower degradation rate for master sequence)
		'''
	import CoreClasses
	import itertools
	rxn_IDs = []
	reaction_list = []
	molecule_list = []
	molecule_dict = {}
	rxn_ID = 0
	#molecule_ID = 0
	monomers = 'AB'
	molecule_list.append('A')
	molecule_dict['A'] = molecule_list.index('A')
	molecule_list.append('B')
	molecule_dict['B'] = molecule_list.index('B')
	# Generate all sequences of length l
	original_rconstant = rconstant
	original_dconstant = dconstant
	sequences = itertools.product(monomers, repeat = length)
	for seq in sequences:
		s = ''.join(seq)
		molecule_list.append(s)
		molecule_dict[s] = molecule_list.index(s)
		if 'rep_landscape' in kwargs:
			master_seq = kwargs['rmaster']
			d = hamming_distance(master_seq, s)
			if kwargs['rep_landscape'] == 'gaussian':
				rconstant = original_rconstant*(1.0 + 10.0*gaussian(d, 0, 2) )
			elif kwargs['rep_landscape'] == 'step':
				if d != 0:
					rconstant = original_rconstant
				else:
					rconstant = 10*original_rconstant
		if 'd_landscape' in kwargs:
			master_seq = kwargs['dmaster']
			d = hamming_distance(master_seq, s)
			if kwargs['d_landscape'] == 'gaussian':
				dconstant = original_dconstant*(1.0 - 0.9*gaussian(d, 0, 2) )
			elif kwargs['d_landscape'] == 'step':
				if d != 0:
					dconstant = original_dconstant
				else:
					dconstant = 0.1*original_dconstant

		# Forward Replication Reaction
		rxn_ID = len(reaction_list)
		reactants =  [molecule_dict['A'], molecule_dict['B'] ] 
		reactant_coeff = [ s.count('A'), s.count('B') ]
		products = [ molecule_dict[s] ]
		product_coeff = [1]
		reaction_list.append( CoreClasses.Reaction(rxn_ID, reactants = reactants, reactant_coeff = reactant_coeff, products = products, product_coeff = product_coeff, constant = rconstant, prop = 'RM2') )
		#print 'Reaction List index: ', rxn_ID, 'Reaction ID: ', reaction_list[rxn_ID].ID 
		rxn_IDs.append(rxn_ID)

		# Spontaneous Reaction
		rxn_ID = len(reaction_list)
		reaction_list.append( CoreClasses.Reaction(rxn_ID, reactants = reactants, reactant_coeff = reactant_coeff, products = products, product_coeff = product_coeff, constant = sconstant, prop = 'STD') )
		rxn_IDs.append(rxn_ID)
		# Degradation Reaction
		rxn_ID = len(reaction_list)
		reaction_list.append( CoreClasses.Reaction(rxn_ID, products = reactants, product_coeff = reactant_coeff, reactants = products, reactant_coeff = product_coeff, constant = dconstant, prop = 'STD') )
		#print 'Reaction List index: ', rxn_ID, 'Reaction ID: ', reaction_list[rxn_ID].ID 
		rxn_IDs.append(rxn_ID)
				
	#print molecule_dict
	newCRS = CoreClasses.CRS(molecule_list = molecule_list, molecule_dict = molecule_dict, reactions = reaction_list)
	#print newCRS.molecule_dict
	fname = ''
	if 'rep_landscape' in kwargs:

		master_seq = kwargs['rmaster']
		fname += 'rmaster_' + str(master_seq)
		if kwargs['rep_landscape'] == 'gaussian':
			fname+= '_gaussian_'	
		elif kwargs['rep_landscape'] == 'step':
			fname+= '_step_'
	if 'd_landscape' in kwargs:
		master_seq = kwargs['dmaster']
		fname += 'dmaster_' + str(master_seq)
		if kwargs['d_landscape'] == 'gaussian':
			fname+= '_gaussian_'
		elif kwargs['d_landscape'] == 'step':
			fname+= '_step_'
	fname += 'length_%i.txt' % length
	newCRS.savetxt(fname)
	print "New CRS saved"

def generate_all_replication8_reactions(length, rconstant = 0.001, sconstant = 1E-5, dconstant = 1.0, **kwargs):
	''' Generates all possible replication reactions for binary molecules of size length assigns all reactions the same constant
	forward reactions are generated twice, once for replication, given RCM propensity function, and Spontaneous given standard propensity
	Arguements:
		- length: length of binary replicators 
		- rconstant: replication rate constant
		- sconstant: Spontaneous formation rate of replicators 
		- dconstant: degradation rate constant for replicators
		Keyword Arguements:
		- rep_landscape: If this is passed the replication rate is based on the assigned landscape, either 'gaussian' or 'flat'
		- rmaster: Master sequence for the replication landscape
		- d_landscape: If this is passed the degradation rate is based on the assigned landscape, either 'gaussian' or 'flat'
		- dmaster: Master sequence for the degradation rate

	Fitness landscapes impose a max difference of 1 order of magnitude (ie 10 times faster replication for the master sequence, 
		or 10 times slower degradation rate for master sequence)
		'''
	import CoreClasses
	import itertools
	rxn_IDs = []
	reaction_list = []
	molecule_list = []
	molecule_dict = {}
	rxn_ID = 0
	#molecule_ID = 0
	monomers = 'AB'
	molecule_list.append('A')
	molecule_dict['A'] = molecule_list.index('A')
	molecule_list.append('B')
	molecule_dict['B'] = molecule_list.index('B')
	# Generate all sequences of length l
	original_rconstant = rconstant
	original_dconstant = dconstant
	sequences = itertools.product(monomers, repeat = length)
	for seq in sequences:
		s = ''.join(seq)
		molecule_list.append(s)
		molecule_dict[s] = molecule_list.index(s)
		if 'rep_landscape' in kwargs:
			master_seq = kwargs['rmaster']
			d = hamming_distance(master_seq, s)
			if kwargs['rep_landscape'] == 'gaussian':
				rconstant = original_rconstant*(1.0 + 10.0*gaussian(d, 0, 2) )
			elif kwargs['rep_landscape'] == 'step':
				if d != 0:
					rconstant = original_rconstant
				else:
					rconstant = 10*original_rconstant
		if 'd_landscape' in kwargs:
			master_seq = kwargs['dmaster']
			d = hamming_distance(master_seq, s)
			if kwargs['d_landscape'] == 'gaussian':
				dconstant = original_dconstant*(1.0 - 0.9*gaussian(d, 0, 2) )
			elif kwargs['d_landscape'] == 'step':
				if d != 0:
					dconstant = original_dconstant
				else:
					dconstant = 0.1*original_dconstant

		# Forward Replication Reaction
		rxn_ID = len(reaction_list)
		reactants =  [molecule_dict['A'], molecule_dict['B'] ] 
		reactant_coeff = [ s.count('A'), s.count('B') ]
		products = [ molecule_dict[s] ]
		product_coeff = [1]
		reaction_list.append( CoreClasses.Reaction(rxn_ID, reactants = reactants, reactant_coeff = reactant_coeff, products = products, product_coeff = product_coeff, constant = rconstant, prop = 'RM8') )
		#print 'Reaction List index: ', rxn_ID, 'Reaction ID: ', reaction_list[rxn_ID].ID 
		rxn_IDs.append(rxn_ID)

		# Spontaneous Reaction
		rxn_ID = len(reaction_list)
		reaction_list.append( CoreClasses.Reaction(rxn_ID, reactants = reactants, reactant_coeff = reactant_coeff, products = products, product_coeff = product_coeff, constant = sconstant, prop = 'STD') )
		rxn_IDs.append(rxn_ID)
		# Degradation Reaction
		rxn_ID = len(reaction_list)
		reaction_list.append( CoreClasses.Reaction(rxn_ID, products = reactants, product_coeff = reactant_coeff, reactants = products, reactant_coeff = product_coeff, constant = dconstant, prop = 'STD') )
		#print 'Reaction List index: ', rxn_ID, 'Reaction ID: ', reaction_list[rxn_ID].ID 
		rxn_IDs.append(rxn_ID)
				
	#print molecule_dict
	newCRS = CoreClasses.CRS(molecule_list = molecule_list, molecule_dict = molecule_dict, reactions = reaction_list)
	#print newCRS.molecule_dict
	fname = ''
	if 'rep_landscape' in kwargs:

		master_seq = kwargs['rmaster']
		fname += 'rmaster_' + str(master_seq)
		if kwargs['rep_landscape'] == 'gaussian':
			fname+= '_gaussian_'	
		elif kwargs['rep_landscape'] == 'step':
			fname+= '_step_'
	if 'd_landscape' in kwargs:
		master_seq = kwargs['dmaster']
		fname += 'dmaster_' + str(master_seq)
		if kwargs['d_landscape'] == 'gaussian':
			fname+= '_gaussian_'
		elif kwargs['d_landscape'] == 'step':
			fname+= '_step_'
	fname += 'length_%i.txt' % length
	newCRS.savetxt(fname)
	print "New CRS saved"


def pick_replicator(dice_roll, rxn, CRS, concentrations, mu = 0.001):
	'''Given a dice_roll and a replication reaction, determine the mutation outcome, return rxn object 
	Arguements:
		- dice_roll: random number between 0 and total mutation propensity
		- rxn: original replication reaction
		- CRS: CRS object
		- concentrations: concentration array containing all replicators and monomer concentrations 
		- mu: per-base mutation rate 
	Return:
		- picked_rxn: a Reaction object containing the new sequence to be produced and monomers to be consumed
		 If not enough resources present to replicate, a null Reaction object is returned

	'''

	checkpoint = 0.0
	seq_found = False
	seq = CRS.molecule_list[rxn.products[0]]
	R_L = len(seq)
	#print "Trying to replicate: ", seq
	reactant_concentrations = concentrations[rxn.reactants]
	replicator_concentration = concentrations[rxn.products]
	reactant_coeff = rxn.reactant_coeff
	#catalyzed_constants = rxn.catalyzed_constants

	#Calculate Propensity
	Ap = rxn.constant 

	nA = reactant_coeff[0] # If you're reading this you should confirm that 'A' is stored at index 0
	nB = reactant_coeff[1] # If you're reading this you should confirm that 'B' is stored at index 1


	binomialA = 0    #Used for calculating the contribution from copying A-residues
	binomialB = 0   #Used for calculating the intermediate of contribution from copying A-residues and B-residues
	q_error = 0.0

	for eA in range(0, nA + 1):
	    #Here eA is the number of errors in copying A-residues
	    if seq_found == True:
	        break

	    binomialA = (math.factorial(nA)/(math.factorial(nA - eA)*math.factorial(eA)))*pow(rxn.constant*reactant_concentrations[0], nA - eA)*pow(rxn.constant*reactant_concentrations[1], eA)  #calculates number of sequences with eA errors in copying A and the resource contribution to these sequences

	    for eB in range(0, nB + 1):
	        # Here eB is the number of errors in copying B-residues

	        if eA == 0 and eB == 0:
	            # Keeps perfect copying probability seperate from copies made with errors
	            
	            q_p = pow(1 - mu, R_L)*pow(rxn.constant*reactant_concentrations[0], nA)*pow(rxn.constant*reactant_concentrations[1], nB)
	            checkpoint += q_p*replicator_concentration
	        else:
	            binomialB = (math.factorial(nB)/(math.factorial(nB - eB)*math.factorial(eB)))*pow(rxn.constant*reactant_concentrations[1], nB - eB)*pow(rxn.constant*reactant_concentrations[0], eB) #adds number of mutants with eB B-errors
	            
	            q_error += pow(mu, eA + eB)*pow(1 - mu, R_L - eA - eB)*binomialA*binomialB
	            checkpoint += q_error*replicator_concentration
	        if checkpoint >= dice_roll:
	        	A_errors = eA
	        	B_errors = eB 
	        	seq_found = True
	        	break


	Astring = 'B'*A_errors + 'A'*(nA - A_errors)
	Bstring = 'A'*B_errors + 'B'*(nB - B_errors)

	Alist = list(Astring)
	Blist = list(Bstring)

	random.shuffle(Alist)
	random.shuffle(Blist)

	new_seq = ''
	for i in range(R_L):
	    if seq[i] == 'A':
	        new_seq += Alist.pop()

	    elif seq[i] == 'B':
	        new_seq += Blist.pop()
	Acount = new_seq.count('A')
	Bcount = new_seq.count('B')
	#print "Mutated Seq: ", new_seq
	if (Acount != 0 and Acount > reactant_concentrations[0]) or (Bcount != 0 and Bcount > reactant_concentrations[1]):

		print 'New check: Not enough food to replicate'
		picked_rxn = CoreClasses.Reaction(-1,products = [0,1], product_coeff = [0,0], reactants =[0, 1], reactant_coeff = [0, 0], prop = 'RCM')

	else:

		new_seq_ID  = CRS.molecule_dict[new_seq]
		picked_rxn  = CoreClasses.Reaction(-1,products = [new_seq_ID], product_coeff = [1], reactants =[0, 1], reactant_coeff = [Acount, Bcount], prop = 'RCM')
	#raw_input("Enter")
	return picked_rxn


def hamming_distance(seq_a, seq_b):
	h = 0
	if len(seq_a) != len(seq_b):

		h += abs(len(seq_a) - len(seq_b) )
		if len(seq_a) > len(seq_b):
			l = len(seq_b)
			for i in range(l):
				if seq_a[i] != seq_b[i]:
					h += 1

		else:
			l = len(seq_a)
			for i in range(l):
				if seq_a[i] != seq_b[i]:
					h += 1

	else:
		l = len(seq_a)
		for i in range(l):
			if seq_a[i] != seq_b[i]:
				h += 1

	return h

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
