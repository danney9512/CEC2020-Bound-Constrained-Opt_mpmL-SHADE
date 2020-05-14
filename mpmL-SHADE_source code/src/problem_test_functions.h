#ifndef TEST_FUNCTIONS__
#define TEST_FUNCTIONS__

#include <string>
#include "alg_individual.h"
#include "problem.h"

const double INF = 1.0e99;
const double EPS = 1.0e-14;
const double E = 2.7182818284590452353602874713526625;
const double PI = 3.1415926535897932384626433832795029;

const std::string FunctionName[10] =  { "Storn_Chebyshev", "Inverse_Hilbert", "Lennard_Jones", "Rastrihin", "Griewangk", "Weierstrass", "Modified_Schwefel", "Expanded_Schaffer_F6", "Happy_Cat", "Ackley" };

// (F1-F3 in CEC 2019 100-Digit Competition).
double Storn_Chebyshev(const Individual::GeneVec &gene_vec);
double Inverse_Hilbert(const Individual::GeneVec &gene_vec);
double Lennard_Jones(const Individual::GeneVec &gene_vec);

// Rotatable and shiftable functions (F4-F10 in CEC 2019 100-Digit Competition).
double Rastrihin(const Individual::GeneVec &gene_vec, const CProblem &prob);
double Griewangk(const Individual::GeneVec &gene_vec, const CProblem &prob);
double Weierstrass(const Individual::GeneVec &gene_vec, const CProblem &prob);
double Modified_Schwefel(const Individual::GeneVec &gene_vec, const CProblem &prob);
double Expanded_Schaffer_F6(const Individual::GeneVec &gene_vec, const CProblem &prob);
double Happy_Cat(const Individual::GeneVec &gene_vec, const CProblem &prob);
double Ackley(const Individual::GeneVec &gene_vec, const CProblem &prob);

void CEC2019_100Digit_Evaluate(Individual *indiv, const CProblem &prob);
const int CEC2019_100Digit_Score(Individual *indiv, const double &optimum);
void CEC2014_Evaluate(Individual *indiv, const CProblem &prob);

// CEC 2014
double Ellips(const Individual::GeneVec &gene_vec, const CProblem &prob);
double Bent_Cigar(const Individual::GeneVec &gene_vec, const CProblem &prob);
double Discus(const Individual::GeneVec &gene_vec, const CProblem &prob);
double Rosenbrock(const Individual::GeneVec &gene_vec, const CProblem &prob);
double Katsuura(const Individual::GeneVec &gene_vec, const CProblem &prob);
double Hgbat(const Individual::GeneVec &gene_vec, const CProblem &prob);
double Grie_Rosen(const Individual::GeneVec &gene_vec, const CProblem &prob);
double hf01(const Individual::GeneVec &gene_vec, const CProblem &prob);
double hf02(const Individual::GeneVec &gene_vec, const CProblem &prob);
double hf03(const Individual::GeneVec &gene_vec, const CProblem &prob);
double hf04(const Individual::GeneVec &gene_vec, const CProblem &prob);
double hf05(const Individual::GeneVec &gene_vec, const CProblem &prob);
double hf06(const Individual::GeneVec &gene_vec, const CProblem &prob);
double cf01(const Individual::GeneVec &gene_vec, const CProblem &prob);
double cf02(const Individual::GeneVec &gene_vec, const CProblem &prob);
double cf03(const Individual::GeneVec &gene_vec, const CProblem &prob);
double cf04(const Individual::GeneVec &gene_vec, const CProblem &prob);
double cf05(const Individual::GeneVec &gene_vec, const CProblem &prob);
double cf06(const Individual::GeneVec &gene_vec, const CProblem &prob);
double cf07(const Individual::GeneVec &gene_vec, const CProblem &prob);
double cf08(const Individual::GeneVec &gene_vec, const CProblem &prob);


// Foundation functions (shift and rotate).
Individual::GeneVec shiftfunc(const Individual::GeneVec &gene_vec, const CProblem::RNV &shift_vec);
Individual::GeneVec rotatefunc(const Individual::GeneVec &gene_vec, const std::vector<CProblem::RNV> &rotate_matrix);
Individual::GeneVec sr_func(const Individual::GeneVec &gene_vec, const CProblem::RNV &shift_vec, const std::vector<CProblem::RNV> &rotate_matrix, const double& shrink_rate, const bool& shift_flag, const bool& rotate_flag);
Individual::GeneVec shuffle_func(const Individual::GeneVec &gene_vec, const std::vector<int> &shuffle_vec);
double cf_cal(const Individual::GeneVec &gene_vec, const CProblem::RNV &shift_vec, const std::vector<double> &fit_vec, const std::vector<double> &delta, const std::vector<double> &bias);



// CEC2020
// add at 2019/12/24
namespace CEC2020 
{
	void CEC2020_BoundConstrained_Evaluate(Individual* indiv, const CProblem& prob);

	double Lunacek_bi_Rastrigin(const Individual::GeneVec& gene_vec, const CProblem& prob);
	double hf01(const Individual::GeneVec& gene_vec, const CProblem& prob);
	double hf05(const Individual::GeneVec& gene_vec, const CProblem& prob);
	double hf06(const Individual::GeneVec& gene_vec, const CProblem& prob);
	double cf02(const Individual::GeneVec& gene_vec, const CProblem& prob);
	double cf04(const Individual::GeneVec& gene_vec, const CProblem& prob);
	double cf05(const Individual::GeneVec& gene_vec, const CProblem& prob);
}



#endif // !TEST_FUNCTIONS__

