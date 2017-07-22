import numpy as np
import PropensityFunctions as Propensity
import CoreClasses as Core 
import random
import OutputFunctions as Out
import InitializeFunctions as Init

####################################################
### Load C library
####################################################
from ctypes import cdll
from ctypes import byref, c_int, c_ulong, c_double, POINTER
def get_libpath():
    """
    Get the library path of the the distributed SSA library.
    """
    import os
    import re
    from os.path import dirname, abspath, realpath, join
    from platform import system

    root = dirname(abspath(realpath(__file__)))
   
    if system() == 'Linux':
        library = 'Linux-SSA.so'
    elif system() == 'Darwin':
        library = 'OSX-SSA.so'
    elif system() == 'Windows':
        library = "Win-SSA.so"
    else:
        raise RuntimeError("unsupported platform - \"{}\"".format(system()))

    return os.path.join(root, 'clibs', library)

_SSA_LIB = cdll.LoadLibrary(get_libpath())
#                               current_t, next_t,  r_seed, max_x, max_y, num_m, num_r, concentrations,     constants         propensity_ints, reaction_arr,  catalyst_arr
_SSA_LIB.SSA_update.argtypes = (c_double, c_double, c_int, c_int, c_int, c_int, c_int, POINTER(c_double), POINTER(c_double), POINTER(c_int), POINTER(c_int), POINTER(c_double))
_SSA_LIB.SSA_update.restype = c_double
SSA_update = _SSA_LIB.SSA_update ### Renaming function for convinence
####################################################
####################################################
def pick_reaction(dice_roll, CRS, concentrations, **kwargs):
	''' Picks a reaction to occur stochastically 

	Arguements:
		- dice_roll: float which should be a number between zero and the total propensity of reactions
		- CRS: the CRS object which contains all possible reactions and molecules
		- concentrations: the list of concentrations indexed by molecule ID
		- propensity_function: which propensity function to use, default: standard

	Return:
		- rxn: a Reaction object'''

	checkpoint = 0.0
	for rxn in CRS.reactions:
		reactant_concentrations = [concentrations[i] for i in rxn.reactants]
		catalyst_concentrations = [concentrations[i] for i in rxn.catalysts]
		reactant_coeff = rxn.reactant_coeff
		catalyzed_constants = rxn.catalyzed_constants
		#print rxn.catalysts
		if rxn.prop == 'STD':
			# print "Reactant concentrations: ", reactant_concentrations
			# print 'Product ID numbers: ',rxn.products
			checkpoint += Propensity.standard_propensity(rxn, CRS, concentrations)
			#print "dice_roll: ", dice_roll, ' checkpoint: ', checkpoint
			if checkpoint >= dice_roll:
				break
			
		elif rxn.prop == 'RCM':
			mu = kwargs['mu']
			checkpoint += Propensity.replicator_composition_propensity_envMutation(rxn, CRS, concentrations, mu = mu)
			if checkpoint >= dice_roll:
				mutation_dice = checkpoint - dice_roll
				rxn = pick_replicator(mutation_dice, rxn,CRS, concentrations, mu)
				break
		elif rxn.prop[:2] == 'MM':
			expon = int(rxn.prop[2])
			kcat = 10**expon
			checkpoint += Propensity.MM_kinetics(rxn, CRS, concentrations, kcat)
			if checkpoint >= dice_roll:
				break
				
		
		
	#raw_input("Enter")
	return rxn
####################################################
def execute_rxn(rxn, CRS, concentrations):
	''' Executes a single reaction instance

	Arguements:
		- rxn: Reaction object to execute_rxn
		- CRS: CRS object containing the entire system
		- concentrations: list of molecule concentrations indexed by ID

	Return:
		- concentrations: updated list of molecule concentrations indexed by ID '''
	num_reactants = len(rxn.reactants)
	
	num_products = len(rxn.products)
	
	# Reduce Reactants
	for i in range(num_reactants):
		reactant_index = rxn.reactants[i]
		concentrations[reactant_index] -= rxn.reactant_coeff[i]
		
	# Increase Products	
	for i in range(num_products):
		product_index =rxn.products[i]
		concentrations[product_index] += rxn.product_coeff[i]
		
	return concentrations
####################################################
def SSA_evolve(tau, tau_max, concentrations, CRS, random_seed, output_prefix= None,  t_out= None):

	if (output_prefix != None and t_out == None):
		raise ValueError('Output file prefix specified but no output frequency given, please provide an output time frequency')
		
	elif (output_prefix == None and type(t_out) == float):
		raise ValueError('Output frequency provided but output file prefix was not provided, please provide a file prefix name')
		
	import sys
	import random
	from ctypes import c_int,  c_double, POINTER
	constants, propensity_ints, reaction_arr, catalyst_arr = Init.convert_CRS_to_npArrays(CRS)
	concentrations_ptr, constants_ptr, propensity_ints_ptr, reaction_arr_ptr, catalyst_arr_ptr= Init.get_c_pointers(concentrations, constants, propensity_ints, reaction_arr, catalyst_arr)
	freq_counter = 0.0
	random.seed(random_seed)
	while tau < tau_max:
		# Get seed
		r_seed = random.randint(0, sys.maxint)
		# Update concentrations in place using C function
		c_tau = SSA_update(c_double(tau), c_double(freq_counter),r_seed, c_int(1),c_int(1), c_int(len(CRS.molecule_list)), c_int(len(constants)), concentrations_ptr, constants_ptr, propensity_ints_ptr, reaction_arr_ptr, catalyst_arr_ptr )
		# Update Time
		tau = c_tau
		# Update random seed
		random.jumpahead(tau-freq_counter)
		print tau
		# Output data
		Out.output_concentrations(concentrations, 'tutorial_data',time = freq_counter)
		freq_counter += t_out
	Out.tidy_timeseries(CRS.molecule_list, 'tutorial_data', delete_dat = True)

	return concentrations

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
