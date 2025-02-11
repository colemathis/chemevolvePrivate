import numpy as np
import math

def standard_propensity(rxn, CRS, concentrations):
	''' Standard Propensity function calculates propensity as the concentrations of the reactants raised to their coefficients 

	Arguements:
		- rxn: Reaction object
		- CRS: CRS object for system
		- concentrations: list of molecule concentrations indexed by ID

	Return:
		- Ap: float, propensity of rxn given the current concentrations

	'''
	# Get data out of rxn and concentration objects
	reactant_concentrations = concentrations[rxn.reactants]
	catalyst_concentrations = concentrations[rxn.catalysts]
	reactant_coeff = rxn.reactant_coeff
	catalyzed_constants = rxn.catalyzed_constants

	#Calculate Propensity
	Ap = rxn.constant #*np.prod( np.power(reactant_concentrations, reactant_coeff) )
	num_reactants = len(reactant_concentrations)
	for i in range(num_reactants):
		Ap = Ap*np.power(reactant_concentrations[i],reactant_coeff[i])
	
	# Ap = rxn.constant
	#Ap = Ap*np.prod(std_propensity_helper1(reactant_concentrations, reactant_coeff))
	
	enhancement = 0.0
	#if catalyst_concentrations != [] and sum(catalyst_concentrations) != 0.0:
	
	num_cats = len(catalyst_concentrations)
	for i in range(num_cats):
		#print 'catalyst: ', rxn.catalysts[i]

		enhancement += catalyst_concentrations[i]*catalyzed_constants[i]
		#print 'enhancement: ', enhancement
		#enhancement += np.sum( std_propensity_helper2(catalyst_concentrations,catalyzed_constants) )
	Ap = Ap*(1+enhancement)
	return Ap


def MM_kinetics(rxn, CRS, concentrations, kcat):
	''' Michaelis Menten kinetic Propensity function calculates propensity according the Michaelis-Menten assumptions
	Only supports unimolecular reactions with a single catalyst
	IMPORTANTLY: value of catalyzed constants now corresponds to K_m from MM kinetic equation
	Arguements:
		- rxn: Reaction object
		- CRS: CRS object for system
		- concentrations: list of molecule concentrations indexed by ID

	Return:
		- Ap: float, propensity of rxn given the current concentrations

	'''
	reactant_concentrations = concentrations[rxn.reactants[0]]
	catalyst_concentrations = concentrations[rxn.catalysts]
	reactant_coeff = rxn.reactant_coeff[0]
	catalyzed_constants = rxn.catalyzed_constants
	if len(rxn.reactants) > 1:
		raise ValueError('Michaelis Menten kinetics only supports unimolecular reactions, # Rectants input %i' % len(rxn.reactants))
	elif reactant_coeff > 1: 
		raise ValueError('Michaelis Menten kinetics only supports unimolecular reactions, Stoichimetric coefficient %i' %len(reactant_coeff[0]))
	# Get data out of rxn and concentration objects
	
	#Calculate Propensity
	Ap = rxn.constant*reactant_concentrations #*np.prod( np.power(reactant_concentrations, reactant_coeff) )
	enhancement = 0.0
	
	num_cats = len(catalyst_concentrations)
	if reactant_concentrations > 0:
		for i in range(num_cats):
			#print 'catalyst: ', rxn.catalysts[i]

			enhancement += kcat*catalyst_concentrations[i]/(catalyzed_constants[i] + reactant_concentrations)
			#print 'enhancement: ', enhancement
			#enhancement += np.sum( std_propensity_helper2(catalyst_concentrations,catalyzed_constants) )
		Ap = Ap*(1+enhancement)
	
	return Ap

def sort_propensities(CRS, concentrations, **kwargs):
	lattice_size = np.product(concentrations.shape[:-1])
	lattice_shape = concentrations.shape[:-1]
	propensity_arr = np.zeros(lattice_size).reshape(lattice_shape)
	
	# Iterate over each lattice site. Site is the propensity and index gives its location
	for site_index, site_Ap in np.ndenumerate(propensity_arr):
		(x,y) = site_index
		# Iterate over each reactation
		for rxn in CRS.reactions:

			if rxn.prop == 'STD':
				Ap = standard_propensity(rxn, CRS, concentrations[site_index])
				
			elif rxn.prop == 'RM8':
				mu = kwargs['mu']
				Ap = replicator_composition_propensity_envMutation8(rxn, CRS, concentrations[site_index], mu = mu)
			elif rxn.prop == 'RM2':
				mu = kwargs['mu']
				Ap = replicator_composition_propensity_envMutation2(rxn, CRS, concentrations[site_index], mu = mu)
			elif rxn.prop == 'RM1':
				mu = kwargs['mu']
				Ap = replicator_composition_propensity_envMutation1(rxn, CRS, concentrations[site_index], mu = mu)
			elif rxn.prop[:2] == 'MM':
				expon = int(rxn.prop[2])
				kcat = 10**expon
				Ap = MM_kinetics(rxn, CRS, concentrations[site_index], kcat)				
			site_Ap += Ap

		propensity_arr[site_index] = site_Ap
		#print site_Ap
	max_Ap = max(propensity_arr)
	for site_index, site_Ap in np.ndenumerate(propensity_arr):
		(x,y) = site_index
		if site_Ap == max_Ap:
			break
	Ap_rxns = []
	for rxn in CRS.reactions:

		if rxn.prop == 'STD':
			Ap = standard_propensity(rxn, CRS, concentrations[site_index])
			Ap_rxns.append( (Ap, rxn))
			
		elif rxn.prop == 'RM2':
			mu = kwargs['mu']
			Ap = replicator_composition_propensity_envMutation2(rxn, CRS, concentrations[site_index], mu = mu)
			Ap_rxns.append( (Ap, rxn))
		elif rxn.prop == 'RM8':
			mu = kwargs['mu']
			Ap = replicator_composition_propensity_envMutation8(rxn, CRS, concentrations[site_index], mu = mu)
			Ap_rxns.append( (Ap, rxn))
		elif rxn.prop == 'RM1':
			mu = kwargs['mu']
			Ap = replicator_composition_propensity_envMutation1(rxn, CRS, concentrations[site_index], mu = mu)
			Ap_rxns.append( (Ap, rxn))

		elif rxn.prop[:2] == 'MM':
			expon = int(rxn.prop[2])
			kcat = 10**expon
			Ap = MM_kinetics(rxn, CRS, concentrations[site_index], kcat)
			Ap_rxns.append( (Ap, rxn))
	Ap, rxn_list = zip(*sorted(Ap_rxns, key=lambda x: x[0], reverse = True)	)
	CRS.reactions = rxn_list		

	return CRS
	
def calculate_propensities(CRS, concentrations, **kwargs):
	''' Calculate the propensity of a reaction according to the concentrations and propensity function

	Arguements: 
		- CRS: CRS object
		- concentrations: array of molecule concentrations indexed by (position,ID) 
		- propensity_function: function option used to calculate reaction propensities

	Return:
		propensity_arr: an array of floats giving the total reaction propensity at each point in the system '''
	# Initialize an Empty Array to store propensities for each lattice site
	lattice_size = np.product(concentrations.shape[:-1])
	lattice_shape = concentrations.shape[:-1]
	propensity_arr = np.zeros(lattice_size).reshape(lattice_shape)
	
	# Iterate over each lattice site. Site is the propensity and index gives its location
	for site_index, site_Ap in np.ndenumerate(propensity_arr):
		(x,y) = site_index
		site_Ap = 0.0
		# Iterate over each reactation
		for rxn in CRS.reactions:
			Ap = 0.0
			if rxn.prop == 'STD':
				Ap = standard_propensity(rxn, CRS, concentrations[site_index])
				
			elif rxn.prop == 'RM2':
				mu = kwargs['mu']
				Ap = replicator_composition_propensity_envMutation2(rxn, CRS, concentrations[site_index], mu = mu)

			elif rxn.prop == 'RM8':
				mu = kwargs['mu']
				Ap = replicator_composition_propensity_envMutation8(rxn, CRS, concentrations[site_index], mu = mu)
			elif rxn.prop == 'RM1':
				mu = kwargs['mu']
				Ap = replicator_composition_propensity_envMutation1(rxn, CRS, concentrations[site_index], mu = mu)

			elif rxn.prop[:2] == 'MM':
				expon = int(rxn.prop[2])
				kcat = 10**expon
				Ap = MM_kinetics(rxn, CRS, concentrations[site_index], kcat)
			site_Ap += Ap

		propensity_arr[site_index] = site_Ap
		#print site_Ap
	
	return propensity_arr


def replicator_composition_propensity_envMutation8(rxn, CRS, concentrations, mu = 0.001):
	''' Replication Propensity function calculates propensity as the concentrations of the replicator and the composition of the enivornment 
		This propensity function calcuates the mutation propensity as a function of the resources availible
	Arguements:
		- rxn: Reaction object
		- CRS: CRS object for system
		- concentrations: list of molecule concentrations indexed by ID

	Return:
		- Ap: float, propensity of rxn given the current concentrations

	'''
	#from ReplicatorFunctions import get_composition

	# Get data out of rxn and concentration objects
	reactant_concentrations = concentrations[rxn.reactants]
	replicator_concentration = concentrations[rxn.products]
	reactant_coeff = rxn.reactant_coeff

	#catalyzed_constants = rxn.catalyzed_constants

	#Calculate Propensity
	Ap = rxn.constant 

	nA = reactant_coeff[0] # If you're reading this you should confirm that 'A' is stored at index 0
	nB = reactant_coeff[1] # If you're reading this you should confirm that 'B' is stored at index 1
	R_L = nA + nB
	if mu != 0:

	    binomialA = 0    #Used for calculating the contribution from copying A-residues
	    binomialB = 0   #Used for calculating the intermediate of contribution from copying A-residues and B-residues
	    q_error = 0.0
	    
	    for eA in range(0, nA + 1):
	        #Here eA is the number of errors in copying A-residues
	        binomialA = (math.factorial(nA)/(math.factorial(nA - eA)*math.factorial(eA)))*pow(reactant_concentrations[0], nA - eA)*pow(reactant_concentrations[1], eA)  #calculates number of sequences with eA errors in copying A and the resource contribution to these sequences

	        for eB in range(0, nB + 1):
	            # Here eB is the number of errors in copying B-residues
	            if eA == 0 and eB == 0:
	                # Keeps perfect copying probability seperate from copies made with errors
	                q_p = pow(1 - mu, R_L)*pow(reactant_concentrations[0], nA)*pow(reactant_concentrations[1], nB)
	                
	            else:
	                binomialB = (math.factorial(nB)/(math.factorial(nB - eB)*math.factorial(eB)))*pow(reactant_concentrations[1], nB - eB)*pow(reactant_concentrations[0], eB) #adds number of mutants with eB B-errors
	                
	                q_error += pow(mu, eA + eB)*pow(1 - mu, R_L - eA - eB)*binomialA*binomialB

	elif mu == 0:
		q_p = pow(1 - mu, R_L)*pow(reactant_concentrations[0], nA)*pow(reactant_concentrations[1], nB)
		q_error = 0

	Ap = Ap*(q_p + q_error)*replicator_concentration 

	
	return Ap

def replicator_composition_propensity_envMutation2(rxn, CRS, concentrations, mu = 0.001):
	''' Replication Propensity function calculates propensity as the concentrations of the replicator and the composition of the enivornment 
		This propensity function calcuates the mutation propensity as a function of the resources availible
	Arguements:
		- rxn: Reaction object
		- CRS: CRS object for system
		- concentrations: list of molecule concentrations indexed by ID

	Return:
		- Ap: float, propensity of rxn given the current concentrations

	'''
	#from ReplicatorFunctions import get_composition

	# Get data out of rxn and concentration objects
	reactant_concentrations = concentrations[rxn.reactants]
	replicator_concentration = concentrations[rxn.products]
	reactant_coeff = rxn.reactant_coeff

	#catalyzed_constants = rxn.catalyzed_constants

	#Calculate Propensity
	Ap = rxn.constant 

	nA = reactant_coeff[0] # If you're reading this you should confirm that 'A' is stored at index 0
	nB = reactant_coeff[1] # If you're reading this you should confirm that 'B' is stored at index 1
	R_L = nA + nB
	if mu != 0:

	    binomialA = 0    #Used for calculating the contribution from copying A-residues
	    binomialB = 0   #Used for calculating the intermediate of contribution from copying A-residues and B-residues
	    q_error = 0.0
	    
	    for eA in range(0, nA + 1):
	        #Here eA is the number of errors in copying A-residues
	        binomialA = (math.factorial(nA)/(math.factorial(nA - eA)*math.factorial(eA)))
	        for eB in range(0, nB + 1):
	            # Here eB is the number of errors in copying B-residues
	            if eA == 0 and eB == 0:
	                # Keeps perfect copying probability seperate from copies made with errors
	                q_p = pow(1 - mu, R_L)*(reactant_concentrations[0]*nA)*(reactant_concentrations[1]*nB)
	                
	            else:
	                binomialB = (math.factorial(nB)/(math.factorial(nB - eB)*math.factorial(eB))) #adds number of mutants with eB B-errors
	                q_error += pow(mu, eA + eB)*pow(1 - mu, R_L - eA - eB)*binomialA*binomialB*( reactant_concentrations[0]*(nA - eA +eB)*reactant_concentrations[1]*(nB - eB + eA ) )

	elif mu == 0:
		q_p = pow(1 - mu, R_L)*(reactant_concentrations[0]*nA)*(reactant_concentrations[1]*nB)
		q_error = 0

	Ap = Ap*(q_p + q_error)*replicator_concentration 

	
	return Ap

def replicator_composition_propensity_envMutation1(rxn, CRS, concentrations, mu = 0.001):
	''' Replication Propensity function calculates propensity as the concentrations of the replicator and the composition of the enivornment 
		This propensity function calcuates the mutation propensity as a function of the resources availible
	Arguements:
		- rxn: Reaction object
		- CRS: CRS object for system
		- concentrations: list of molecule concentrations indexed by ID

	Return:
		- Ap: float, propensity of rxn given the current concentrations

	'''
	#from ReplicatorFunctions import get_composition

	# Get data out of rxn and concentration objects
	reactant_concentrations = concentrations[rxn.reactants]
	replicator_concentration = concentrations[rxn.products]
	reactant_coeff = rxn.reactant_coeff

	#catalyzed_constants = rxn.catalyzed_constants

	#Calculate Propensity
	Ap = rxn.constant 

	nA = reactant_coeff[0] # If you're reading this you should confirm that 'A' is stored at index 0
	nB = reactant_coeff[1] # If you're reading this you should confirm that 'B' is stored at index 1
	R_L = float(nA + nB)
	if mu != 0:

	    binomialA = 0    #Used for calculating the contribution from copying A-residues
	    binomialB = 0   #Used for calculating the intermediate of contribution from copying A-residues and B-residues
	    q_error = 0.0
	    
	    for eA in range(0, nA + 1):
	        #Here eA is the number of errors in copying A-residues
	        binomialA = (math.factorial(nA)/(math.factorial(nA - eA)*math.factorial(eA)))*pow(reactant_concentrations[0], (nA - eA)/R_L)*pow(reactant_concentrations[1], eA/R_L)  #calculates number of sequences with eA errors in copying A and the resource contribution to these sequences

	        for eB in range(0, nB + 1):
	            # Here eB is the number of errors in copying B-residues
	            if eA == 0 and eB == 0:
	                # Keeps perfect copying probability seperate from copies made with errors
	                q_p = pow(1 - mu, R_L)*pow(reactant_concentrations[0], nA/R_L)*pow(reactant_concentrations[1], nB/R_L)
	                
	            else:
	                binomialB = (math.factorial(nB)/(math.factorial(nB - eB)*math.factorial(eB)))*pow(reactant_concentrations[1], (nB - eB)/R_L )*pow(reactant_concentrations[0], eB/R_L) #adds number of mutants with eB B-errors
	                
	                q_error += pow(mu, eA + eB)*pow(1 - mu, R_L - eA - eB)*binomialA*binomialB

	elif mu == 0:
		q_p = pow(1 - mu, R_L)*pow(reactant_concentrations[0], nA/R_L)*pow(reactant_concentrations[1], nB/R_L)
		q_error = 0

	Ap = Ap*(q_p + q_error)*replicator_concentration 

	
	return Ap