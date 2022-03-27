#include "Utils.h"
#include "BackgroundCosmology.h"


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
	

	
	

	Utils::EndTiming("Everything");
}
