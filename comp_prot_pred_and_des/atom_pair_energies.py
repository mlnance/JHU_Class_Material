# @file: atom_pair_energies.py
# @author: Mike Pacella and Rebecca Alford
# @brief: Print atom pair energies for physical terms
# @note: Tested 2/22/13

# By Michael Pacella
from pyrosetta import rosetta
import csv

def etable_atom_pair_energies(res1, atom_index_1, res2, atom_index_2, sfxn):
    """
    Usage: lj_atr, lj_rep, solv=etable_atom_pair_energies(atom1, atom2, sfxn)
        Description: given a pair of atoms and scorefunction, use the precomputed
        'etable' to return LJ attractive, LJ repulsive, and LK solvation energies
    """
    score_manager = rosetta.core.scoring.ScoringManager.get_instance()
    etable_ptr = score_manager.etable( sfxn.energy_method_options().etable_type() )
    etable = etable_ptr.lock()
    etable_energy = rosetta.core.scoring.etable.AnalyticEtableEnergy(etable,
                                                  sfxn.energy_method_options())
        
        # Construct coulomb class for calculating fa_elec energies                                           
    coulomb = rosetta.core.scoring.etable.coulomb.Coulomb(sfxn.energy_method_options())

        # Construct AtomPairEnergy container to hold computed energies.
    ape = rosetta.core.scoring.etable.AtomPairEnergy()

        # Set all energies in the AtomPairEnergy to zero prior to calculation.
    ape.attractive, ape.bead_bead_interaction, ape.repulsive, ape.solvation = \
                                                             0.0, 0.0, 0.0, 0.0

        # Calculate the distance squared and set it in the AtomPairEnergy.
    ape.distance_squared = res1.xyz(atom_index_1).distance_squared(res2.xyz(atom_index_2))

        # Evaluate energies from pre-calculated etable, using a weight of 1.0
        # in order to match the raw energies from eval_ci_2b.
    atom1 = res1.atom(atom_index_1)
    atom2 = res2.atom(atom_index_2) 
    etable_energy.atom_pair_energy(atom1, atom2, 1.0, ape)

        # Calculate atom-atom scores.
    lj_atr = ape.attractive
    lj_rep = ape.repulsive
    solv = ape.solvation
    fa_elec = coulomb.eval_atom_atom_fa_elecE(res1.xyz(atom_index_1),res1.atomic_charge(atom_index_1), \
    			res2.xyz(atom_index_2), res2.atomic_charge(atom_index_2))
    
    return lj_atr, lj_rep, solv, fa_elec
    
def print_atom_pair_energy_table(sfxn, score_type, residue_1, residue_2, output_filename):
	 
	atr_total, rep_total, solv_total, fa_elec_total = 0.0, 0.0, 0.0, 0.0

	list_of_res1_atoms = []
	header_list = []
	header_list.append(' ')
	
	for atom in range(residue_2.natoms()):
		header_list.append(residue_2.atom_name(atom+1))
	list_of_res1_atoms.append(header_list)
	
	for i in range(residue_1.natoms()):
		list_of_interactions = []
		list_of_interactions.append(residue_1.atom_name(i+1))
	
		for j in range(residue_2.natoms()):
		
			atr, rep ,solv, fa_elec = etable_atom_pair_energies(residue_1, i+1, \
							residue_2, j+1, sfxn)
			
			if score_type == 'fa_elec':
				list_of_interactions.append(fa_elec)
			elif score_type == 'fa_atr':
				list_of_interactions.append(atr)
			elif score_type == 'fa_rep':
				list_of_interactions.append(rep)
			elif score_type == 'fa_sol':
				list_of_interactions.append(solv)
			else:
				print("please enter a valid score_type: fa_elec, fa_atr, fa_rep, or fa_sol")
				return 
		list_of_res1_atoms.append(list_of_interactions)
	
	with open(output_filename + "_" + str(score_type) + ".csv", "wb") as f:
		writer = csv.writer(f)
		writer.writerows(list_of_res1_atoms)
		
def print_residue_pair_energies(res, pose, sfxn, score_type, output_filename):
	list_of_scores = []
	sfxn.score(pose)

	for i in range(1, pose.total_residue()+1):
		emap = rosetta.core.scoring.EMapVector()
		sfxn.eval_ci_2b(pose.residue(res),pose.residue(i),pose,emap)
		if abs(emap[score_type]) > 0.0:
			list_of_scores.append([i,emap[score_type]])
		
	
		with open(output_filename + "_" + str(score_type) + ".csv", "wb") as f:
			writer = csv.writer(f)
			writer.writerows(list_of_scores)	
	
