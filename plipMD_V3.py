import MDAnalysis as md

import pandas as pd

from plip.structure.preparation import PDBComplex

import progressbar

import os
import sys
import warnings
warnings.filterwarnings("ignore")


def plipmd(topol=None,traj=None):

	traj=list(traj.strip('[]').split(','))
	
	u = md.Universe(topol,traj)

	ligand_name=input ('\n\n1) Type the ResName of your Ligand (must be 3 letter code -example: LIG -):\n>')

	sol_name=input ('\n2) Type the ResName of your Water (must be 3-4 letter code -example: WAT or SOL or TIP3P -):\n>')

	for res in u.residues:
		if res.resname==ligand_name:
			res.resname='LIG'
		if res.resname==sol_name:
			res.resname='HOH'
		if 'HI' in res.resname:
			res.resname='HIS'
		if 'CY' in res.resname:
			res.resname='CYS'
	for atom in u.atoms:
		if atom.name=='OH2':
			atom.name='OW'

	System=u.select_atoms('protein or (resname LIG or resname HOH)')
	System=System.select_atoms('protein or resname LIG or (around 7 resname LIG)',updating=True)

	
	for ts in u.trajectory[0:1]:
		name='frame_tmp.pdb'
		PDB= md.Writer(name, multiframe=False)
		PDB.write(System)
		plip_job = PDBComplex()
		plip_job.load_pdb(name) 
		plip_job.analyze()
		print ('\nINFO:',plip_job,'\n')
		ligand=input('3) Type the name of the ligand in trajectory to analyze (- example: LIG:S:152 -):\n>')
	os.remove(name)

	table=pd.DataFrame()
	index=0
	print ('\nINFO: Your trajectory lenght is:{} steps\n'.format(range(len(u.trajectory))))
	start=int(input('4) Type the starting STEP to analyze:\n>'))
	finish=int(input('\n5) Type the ending STEP to analyze:\n>'))
	bar=progressbar.ProgressBar(max_value=finish)
	print ('\n\n-----  -----  -----  RUNNING THE ANALYSIS  -----  -----  -----\n\n')
	for i in range(start,finish):
		name='frame_tmp.pdb'
		PDB= md.Writer(name, multiframe=False)
		for ts in u.trajectory[i:i+1]:
			PDB.write(System)
			plip_job = PDBComplex()
			plip_job.load_pdb(name) 
			plip_job.analyze()
			interactions = plip_job.interaction_sets[ligand]
			for interaction in interactions.all_itypes:
				interaction_type=str(type(interaction)).split('.')[-1].replace("'>","")
				table.loc[index,'Frame']=ts.frame
				table.loc[index,'Time']=ts.time
				table.loc[index,'Type']=interaction_type
				table.loc[index,'Residue']=interaction.restype+str(interaction.resnr)
				table.loc[index,'Ligand']=interaction.restype_l+str(interaction.resnr_l)
				
				if interaction_type == 'hbond':
					table.loc[index,'Acceptor']=interaction.atype
					table.loc[index,'AcceptorIdx']=interaction.a.idx
					table.loc[index,'Donor']=interaction.dtype
					table.loc[index,'DonorIdx']=interaction.d.idx
					table.loc[index,'DistanceAD']=interaction.distance_ad
					table.loc[index,'DistanceAH']=interaction.distance_ah
					table.loc[index,'Angle']=interaction.angle
					table.loc[index,'Force']=interaction.type
					table.loc[index,'ProtIsDon']=interaction.protisdon
				
				elif interaction_type == 'pication':
					table.loc[index,'Charge']=interaction.charge.type
					table.loc[index,'ChargedAtoms']=",".join([i.type for i in interaction.charge.atoms])
					table.loc[index,'Force']=interaction.type
					table.loc[index,'RingType']=interaction.ring.type
					table.loc[index,'RingAtoms']=",".join([i.type for i in interaction.ring.atoms])
					table.loc[index,'RingAtomsIdx']=",".join([str(i.idx) for i in interaction.ring.atoms])        
				
				elif interaction_type=='saltbridge':
					table.loc[index,'NegAtoms']=",".join([i.type for i in interaction.negative.atoms])
					table.loc[index,'NegAtomsIdx']=",".join([str(i.idx) for i in interaction.negative.atoms])
					table.loc[index,'PosAtoms']=",".join([i.type for i in interaction.positive.atoms])
					table.loc[index,'PosAtomsIdx']=",".join([str(i.idx) for i in interaction.positive.atoms])
					table.loc[index,'Distance']=interaction.distance
					table.loc[index,'ProtIsPos']=interaction.protispos
					
				elif interaction_type == 'hydroph_interaction':
					table.loc[index,'RecAtom']=interaction.bsatom.type
					table.loc[index,'RecAtomIdx']=interaction.bsatom.idx
					table.loc[index,'LigAtom']=interaction.ligatom.type
					table.loc[index,'LigAtomIdx']=interaction.ligatom.idx
					table.loc[index,'Distance']=interaction.distance
					
				elif interaction_type == 'pistack':
					table.loc[index,'RecRingType']=interaction.proteinring.type
					table.loc[index,'LigRingType']=interaction.ligandring.type
					table.loc[index,'RecRingAtoms']=",".join([i.type for i in interaction.proteinring.atoms])
					table.loc[index,'RecAtomsIdx']=",".join([str(i.idx) for i in interaction.proteinring.atoms])
					table.loc[index,'LigRingAtoms']=",".join([i.type for i in interaction.ligandring.atoms])
					table.loc[index,'LigRingAtomsIdx']=",".join([str(i.idx) for i in interaction.ligandring.atoms])
					table.loc[index,'Distance']=interaction.distance
					table.loc[index,'Angle']=interaction.angle
					table.loc[index,'Offset']=interaction.offset
					
				elif interaction_type == 'waterbridge':
					table.loc[index,'AccType']=interaction.atype
					table.loc[index,'DonType']=interaction.dtype
					table.loc[index,'WaterIdx']=interaction.water_orig_idx
					table.loc[index,'DistanceAWat']=interaction.distance_aw
					table.loc[index,'DistanceDWat']=interaction.distance_dw
					table.loc[index,'AngleDon']=interaction.d_angle
					table.loc[index,'AngleWat']=interaction.w_angle
					table.loc[index,'ProtIsDon']=interaction.protisdon

				elif interaction_type == 'halogenbond':
					table.loc[index,'Acceptor']=interaction.acctype
					table.loc[index,'Donor']=interaction.acctype
					table.loc[index,'Distance']=interaction.distance
					table.loc[index,'DonAngle']=interaction.don_angle
					table.loc[index,'AccAngle']=interaction.acc_angle
	  
				elif interaction_type=='metal_complex':
					table.loc[index,'MetalType']=interaction.metal.type
					table.loc[index,'Idx']=interaction.metal.idx
					table.loc[index,'TargetType']=interaction.target_type
					table.loc[index,'FunctGroup']=interaction.target.fgroup
					table.loc[index,'Geometry']=interaction.geometry
					table.loc[index,'Distance']=interaction.distance
					table.loc[index,'Location']=interaction.location
				
				index=index+1    
		bar.update(i+1)
		os.remove(name)
		
	print ('\n\n-----  -----  -----  SAVING THE RESULTS, PLEASE WAIT  -----  -----  -----\n\n')	
	table.set_index(['Frame','Time'], inplace=True)
	table.sort_index(inplace=True)
	table.to_excel('Interactions_Table.xlsx')
	print ('\n\n***** ***** ***** ALL DONE, DATA SAVED ON: Interactions_Table.xlsx ***** ***** *****\n\n')	 
if __name__ == "__main__":
	plipmd(sys.argv[1],sys.argv[2])
