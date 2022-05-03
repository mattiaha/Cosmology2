#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"
#include "Perturbations.h"

int main(int argc, char** argv) {
	Utils::StartTiming("Everything");

	//=========================================================================
	// Parameters
	//=========================================================================
	
	// Background parameters
	double h = 0.67;
	double OmegaB = 0.05;
	double OmegaCDM = 0.267;
	double OmegaK = 0.0;
	double Neff = 0.0;
	double TCMB = 2.7255;

	// Recombination parameters
	double Yp = 0.0;

	// Power-spectrum parameters
	double A_s = 2.1e-9;
	double n_s = 0.965;
	double kpivot_mpc = 0.05;

	//=========================================================================
	// Module I
	//=========================================================================

	// Set up and solve the background
	BackgroundCosmology cosmo(h, OmegaB, OmegaCDM, OmegaK, Neff, TCMB);
	cosmo.solve();
	cosmo.info();

	// Output background evolution quantities
	cosmo.output("cosmology.txt");

	//=========================================================================
	// Module II
	//=========================================================================
	RecombinationHistory rec(&cosmo, Yp);
	rec.solve();
	rec.info();

	// Output recombination quantities
	rec.output("recombination.txt");

	// Solve the perturbations
	Perturbations pert(&cosmo, &rec);
	pert.solve();
	pert.info();

	// Output perturbation quantities
	double kvalue1 = 0.001 / Constants.Mpc;
	pert.output(kvalue1, "perturbations_k0.001.txt");
	double kvalue2 = 0.01 / Constants.Mpc;
	pert.output(kvalue2, "perturbations_k0.01.txt");
	double kvalue3 = 0.1 / Constants.Mpc;
	pert.output(kvalue3, "perturbations_k0.1.txt");
	// Remove when module is completed
	return 0;
	

	Utils::EndTiming("Everything");
}