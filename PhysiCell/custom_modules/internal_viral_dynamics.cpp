#include "./internal_viral_dynamics.h" 

using namespace PhysiCell; 

Submodel_Information internal_viral_dynamics_info; 

void internal_virus_model_setup( void )
{
	internal_viral_dynamics_info.name = "internal viral dynamics"; 
	internal_viral_dynamics_info.version = "0.2.0";
	internal_viral_dynamics_info.main_function= internal_virus_model; 

	// what custom data do I need? 
	
	internal_viral_dynamics_info.cell_variables.push_back( "virion" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "uncoated virion" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "viral RNA" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "viral protein" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "assembled virion" ); 

	// submodel_registry.register_model( internal_viral_dynamics_info ); 
	internal_viral_dynamics_info.register_model();
	
	return; 
}
void internal_virus_model( Cell* pCell, Phenotype& phenotype, double dt )
{
	// bookkeeping -- find microenvironment variables we need

	static int nE = microenvironment.find_density_index( "virion" ); 
	static int nUV = microenvironment.find_density_index( "uncoated virion" ); 
	static int nR = microenvironment.find_density_index( "viral RNA" ); 
	static int nP = microenvironment.find_density_index( "viral protein" ); 
	static int nA = microenvironment.find_density_index( "assembled virion" ); 

	// bookkeeping -- find custom data we need 
	
	
	// actual model goes here 
	
	// uncoat endocytosed virus
	double dE = dt * pCell->custom_data["virion_uncoating_rate"] *  ( phenotype.molecular.internalized_total_substrates[nE] ); 
	if( dE > phenotype.molecular.internalized_total_substrates[nE] )
	{ dE = phenotype.molecular.internalized_total_substrates[nE]; } 
	phenotype.molecular.internalized_total_substrates[nE] -= dE; 
	phenotype.molecular.internalized_total_substrates[nUV] += dE; 
	
	// convert uncoated virus to usable mRNA 
	double dR = dt * pCell->custom_data["uncoated_to_RNA_rate"] *  ( phenotype.molecular.internalized_total_substrates[nUV] ); 
	if( dR > phenotype.molecular.internalized_total_substrates[nUV] )
	{ dR = phenotype.molecular.internalized_total_substrates[nUV]; } 
	phenotype.molecular.internalized_total_substrates[nUV] -= dR; 
	phenotype.molecular.internalized_total_substrates[nR] += dR; 
	
	// use mRNA to create viral protein 
	double dP = dt * pCell->custom_data["protein_synthesis_rate"] *  ( phenotype.molecular.internalized_total_substrates[nR] ); 
	phenotype.molecular.internalized_total_substrates[nP] += dP; 
	
	// degrade mRNA 


	// degrade protein 
	
	
	// assemble virus 
	double dA = dt * pCell->custom_data["virion_assembly_rate"] *  ( phenotype.molecular.internalized_total_substrates[nP] ); 
	if( dA > phenotype.molecular.internalized_total_substrates[nP] )
	{ dA = phenotype.molecular.internalized_total_substrates[nP]; } 
	phenotype.molecular.internalized_total_substrates[nP] -= dA; 
	phenotype.molecular.internalized_total_substrates[nA] += dA; 
	
	
	// set export rate 
	
	phenotype.secretion.net_export_rates[nA] = pCell->custom_data["virion_export_rate" ] *  ( phenotype.molecular.internalized_total_substrates[nA] ); 
	
	
	
	// record data: 
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	return; 
}
	
