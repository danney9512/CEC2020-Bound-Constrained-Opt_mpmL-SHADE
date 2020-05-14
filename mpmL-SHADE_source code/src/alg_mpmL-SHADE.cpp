#include "alg_mpmL-SHADE.h"
#include "alg_population.h"
#include "alg_initialization.h"
#include "alg_mutation.h"
#include "alg_crossover.h"
#include "alg_math.h"
#include "problem.h"
#include "problem_test_functions.h"
#include "alg_log.h"

#include <string>
#include <cstdlib>
#include <cmath>
#include <random>
#include <iostream>
#include <climits>
using namespace std;

//-------------------------------
// mpmL-SHADE
//-------------------------------
void mpmL_SHADE::Setup(std::ifstream& ifile)
{
	if (!ifile) return;

	// Loading and Setting
	string dummy;
	ifile >> dummy >> dummy >> name_;
	ifile >> dummy >> dummy >> rNinit_;
	ifile >> dummy >> dummy >> nmin_;
	ifile >> dummy >> dummy >> max_nfe_;
	ifile >> dummy >> dummy >> finit_;
	ifile >> dummy >> dummy >> crinit_;
	ifile >> dummy >> dummy >> rarc_;
	ifile >> dummy >> dummy >> p_;
	ifile >> dummy >> dummy >> h_;
	ifile >> dummy >> dummy >> scalemax_;
}

void mpmL_SHADE::MemorySystem::update_memory(const std::vector<Memory>& success_parameter, const std::vector<double>& success_fit_dif)
{
	if (success_parameter.size() != 0)
	{
		// Calculate weights for each success parameter
		double total_improve = 0.0;
		for (size_t i = 0; i < success_fit_dif.size(); ++i) total_improve += success_fit_dif[i];
		std::vector<double> weight(success_fit_dif.size());
		for (size_t i = 0; i < success_fit_dif.size(); ++i) weight[i] = success_fit_dif[i] / total_improve;

		// weighted Lehmer mean
		double pow_sum_f = 0.0, sum_f = 0.0, pow_sum_cr = 0.0, sum_cr = 0.0;
		for (size_t i = 0; i < success_parameter.size(); ++i)
		{
			pow_sum_f += weight[i] * success_parameter[i].F() * success_parameter[i].F();
			sum_f += weight[i] * success_parameter[i].F();

			pow_sum_cr += weight[i] * success_parameter[i].CR() * success_parameter[i].CR();
			sum_cr += weight[i] * success_parameter[i].CR();
		}

		// Update memory
		if (sum_f == 0) memories_[k_].F() = 0.0;
		else memories_[k_].F() = pow_sum_f / sum_f;
		if (sum_cr == 0) memories_[k_].CR() = 0.0;
		else memories_[k_].CR() = pow_sum_cr / sum_cr;
		++k_;
		if (k_ >= h_) k_ = 0;
	}
}

void mpmL_SHADE::MemorySystem::pertub_memory()
{
	memories_[k_].F() = alg_math::randDouble(0.0, 1.0);
	memories_[k_].CR() = alg_math::randDouble(0.0, 1.0);
	++k_;
	if (k_ >= h_) k_ = 0;
}

void mpmL_SHADE::SubPopGeneration(SubPopulation& sub_pop, const CProblem& prob, const std::size_t& NP, const std::size_t& A, default_random_engine& generator, unsigned long long int& nfe)
{

	const unsigned long long int Max_NFE = max_nfe_;
	const double pbest = p_, ub = prob.upper_bound(), lb = prob.lower_bound(), scalemax = scalemax_;
	const size_t H = h_, Gene_Len = prob.dim();

	sub_pop.sample_parameter().clear();

	if (!sub_pop.success_parameter().empty())
	{
		sub_pop.success_parameter().clear();
		sub_pop.success_fit_dif().clear();
	}

	double scale = 0.1 + (scalemax - 0.1) * (nfe * 1.0 / Max_NFE);
	// For each individual in current popultaion
	for (size_t i = 0; i < NP; ++i)
	{
		// Generate CRi and Fi
		int r = alg_math::randInt(0, (int)(H - 1));

		std::normal_distribution<double> randn(sub_pop.memory_sys()[r].CR(), 0.1);
		double CRi = std::min(1.0, std::fabs(randn(generator)));

		std::cauchy_distribution<double> randc(sub_pop.memory_sys()[r].F(), scale);
		double Fi = randc(generator);
		while (Fi <= 0.0) Fi = randc(generator);
		if (Fi > 1.0) Fi = 1.0;

		sub_pop.sample_parameter().push_back(Memory(Fi, CRi));

		// cur-to-pbest mutation
		Individual::GeneVec donor_vector;
		donor_vector = CurtopBest_DonorVec((int)i, pbest, Fi, sub_pop.pop(), sub_pop.archive_pop());

		// DE crossover
		int R = alg_math::randInt(0, (int)(Gene_Len - 1));
		for (size_t k = 0; k < Gene_Len; ++k)
		{
			sub_pop.children()[i].gene()[k] = (alg_math::randDouble(0.0, 1.0) <= CRi || k == R) ? donor_vector[k] : sub_pop.pop()[i].gene()[k];

			// Repair scheme
			if (sub_pop.children()[i].gene()[k] > ub) sub_pop.children()[i].gene()[k] = (sub_pop.pop()[i].gene()[k] + ub) / 2;
			else if (sub_pop.children()[i].gene()[k] < lb) sub_pop.children()[i].gene()[k] = (sub_pop.pop()[i].gene()[k] + lb) / 2;
		}
		// Evaluate the individual's fitness
		CEC2020::CEC2020_BoundConstrained_Evaluate(&sub_pop.children()[i], prob);
		++nfe;

		if (nfe >= Max_NFE)
		{
			sub_pop.stop_flag() = true;
			sub_pop.stop_idx() = i;
			break;
		}
	}

	// DE Selection
	for (size_t i = 0; i < NP; ++i)
	{
		if (sub_pop.stop_flag() && i > sub_pop.stop_idx()) break;

		if (sub_pop.children()[i].fitness() <= sub_pop.pop()[i].fitness())
		{
			if (sub_pop.children()[i].fitness() < sub_pop.pop()[i].fitness())
			{
				sub_pop.archive_pop().push_back(sub_pop.pop()[i]);
				sub_pop.success_parameter().push_back(Memory(sub_pop.sample_parameter()[i].F(), sub_pop.sample_parameter()[i].CR()));
				sub_pop.success_fit_dif().push_back(std::fabs(sub_pop.children()[i].fitness() - sub_pop.pop()[i].fitness()));
			}
			sub_pop.pop()[i] = sub_pop.children()[i];
		}
	}


	//PolyMutation
	PolynomialMutation PolyMut20(1.0 / prob.dim(), 20); //dig
	PolynomialMutation PolyMut5(1.0 / prob.dim(), 5);  //exploration

	bool Improvement = true;
	while (Improvement == true && sub_pop.stop_flag() == false)
	{
		Improvement = false;
		sub_pop.pop().sort();
		for (size_t i = 1; i < sub_pop.pop().size(); i += 1)
		{
			if (sub_pop.pop()[i].fitness() == sub_pop.pop()[i - 1].fitness())
			{
				if (alg_math::randDouble(0.0, 1.0) <= nfe * 1.0 / Max_NFE)
				{
					while (PolyMut20(&sub_pop.pop()[i]) == false);
				}
				else
				{
					while (PolyMut5(&sub_pop.pop()[i]) == false);
				}
				CEC2020::CEC2020_BoundConstrained_Evaluate(&sub_pop.pop()[i], prob);
				nfe += 1;
				if (nfe >= Max_NFE)
				{
					sub_pop.stop_flag() = true;
					break;
				}
				Improvement = true;
			}
		}
	}
	sub_pop.pop().sort();
	
	// Resize archive 
	if (sub_pop.archive_pop().size() > A)
	{
		sub_pop.archive_pop().shuffle();
		sub_pop.archive_pop().resize(A);
	}
	// Using the success individuals' F and CR to calculate the new lehmer mean
	if (sub_pop.success_parameter().size() > 0)
	{
		sub_pop.memory_sys().update_memory(sub_pop.success_parameter(), sub_pop.success_fit_dif());
		sub_pop.total_no_success() += sub_pop.no_success_cnt();
		sub_pop.no_success_cnt() = 0;
	}
	else
	{
		sub_pop.no_success_cnt() += NP;
		if (sub_pop.max_no_success_cnt() < sub_pop.no_success_cnt()) sub_pop.max_no_success_cnt() = sub_pop.no_success_cnt();

		if (alg_math::randDouble(0.0, 1.0) <= sub_pop.no_success_cnt() * 1.0 / nfe)
		{
			sub_pop.memory_sys().pertub_memory();
			++sub_pop.perturb_cnt();
			// Reset the counter
			sub_pop.total_no_success() += sub_pop.no_success_cnt();
			sub_pop.no_success_cnt() = 0;
		}
	}
}

double EuclideanDistance(const Individual& ind1, const Individual& ind2)
{
	if (ind1.gene().size() == ind2.gene().size())
	{
		double EucDis = 0;
		for (size_t i = 0; i < ind1.gene().size(); i += 1)
		{
			EucDis += (ind1.gene()[i] - ind2.gene()[i]) * (ind1.gene()[i] - ind2.gene()[i]);
		}
		EucDis = sqrt(EucDis);
		return EucDis;
	}
	else
		return 0;
}

void mpmL_SHADE::Solve(Population* solutions, const CProblem& prob)
{
	// Set log system and target problem 
	Log log(name() + '_' + to_string(prob.id()) + "_D" + to_string(prob.dim()) + "_" + prob.name(), max_nfe_, prob.dim(), prob.global_optimum());
	Individual::SetTargetProblem(prob);

	// Set parameter and random generator
	size_t sub_num = ceil(prob.dim() * 0.2);
	//const size_t NPinit = (size_t)round(prob.dim() * rNinit_), Gene_Len = prob.dim(), H = h_, Nmin = nmin_;
	const size_t NPinit = (size_t)round(prob.dim() * 18), Gene_Len = prob.dim(), H = h_, Nmin = nmin_;
	const unsigned long long int Max_NFE = max_nfe_;
	unsigned long long int nfe = 0;
	size_t NP = NPinit, A = (size_t)(NPinit * rarc_);
	const double pbest = p_, ub = prob.upper_bound(), lb = prob.lower_bound(), scalemax = scalemax_;
	std::default_random_engine generator;

	// Initialize and evaluate population  
	Population big_pop(NP * sub_num);
    
	vector<SubPopulation> sub_pops(sub_num, SubPopulation(NP, H));

	RandomInitialization(&big_pop, prob);
	for (size_t i = 0; i < sub_num * NP; ++i)
	{
		CEC2020::CEC2020_BoundConstrained_Evaluate(&big_pop[i], prob);
		++nfe;
	}
	big_pop.sort();
	
	//log.store_pop(&big_pop, (int)nfe);

	//Randomly generate reference point
	big_pop.shuffle();

	//Euclidean distance version
	vector<bool> isassign_bigpop(NP * sub_num, false);
	vector<vector<double>> EucDis(NP * sub_num, vector<double>(NP * sub_num, 0.0));
	for (size_t i = 0; i < NP * sub_num; i += 1)
	{
		for (size_t j = i + 1; j < NP * sub_num; j += 1)
		{
			EucDis[i][j] = EuclideanDistance(big_pop[i], big_pop[j]);
			EucDis[j][i] = EucDis[i][j];
		}
	}

	for (size_t i = 0; i < sub_num; i += 1)
	{
		sub_pops[i].pop()[0] = big_pop[i];
		isassign_bigpop[i] = true;
	}

	for (size_t i = 0; i < sub_num; i += 1)
	{
		size_t BasePoint = i;
		for (size_t j = 1; j < NP; j += 1)
		{
			size_t nearest = BasePoint;
			for (size_t k = 0; k < NP * sub_num; k += 1)
			{
				if (isassign_bigpop[k] == false)
				{
					if (nearest == BasePoint)
					{
						nearest = k;
					}
					else if (EucDis[BasePoint][k] < EucDis[BasePoint][nearest])
					{
						nearest = k;
					}
				}
			}
			sub_pops[i].pop()[j] = big_pop[nearest];
			isassign_bigpop[nearest] = true;
		}
		sub_pops[i].pop().sort();
	}

	// Initialize memory system
	for (size_t i = 0; i < sub_num; ++i)
	{
		for (size_t j = 0; j < H; ++j)
		{
			sub_pops[i].memory_sys()[j] = Memory(finit_, crinit_);
		}
	}

	//log.store(&pop, memory_sys, sample_parameter, success_parameter, (int)nfe);

	size_t cycle_cnt = 0;

	while (nfe < Max_NFE)
	{
		cycle_cnt += 1;
		for (size_t i = 0; i < sub_num; ++i)
		{
			SubPopGeneration(sub_pops[i], prob, NP, A, generator, nfe);
			
			if (sub_pops[i].stop_flag() == true)
			{
				for (size_t j = 0; j < sub_num; ++j)
				{
					sub_pops[j].stop_flag() = true;
				}
				break;
			}
		}

		size_t best_sub_num = 0;
		for (size_t i = 1; i < sub_num; i += 1)
		{
			if (sub_pops[i].pop()[0].fitness() < sub_pops[best_sub_num].pop()[0].fitness())
				best_sub_num = i;
		}

		log.store_bestfitness(&sub_pops[best_sub_num].pop(), (int)nfe);

		if (log.record_point((int)nfe))
		{
			log.store_errorvalue(&sub_pops[best_sub_num].pop());

			/* // store every sub-populations
			for (int i = 0; i < sub_num; i += 1)
			{
				log.store_subpop(&sub_pops[i].pop(), i);
			}
			*/
		}

		//log.store(&pop, memory_sys, sample_parameter, success_parameter, (int)nfe);
		//log.store_pop(&pop, &archive_pop, (int)nfe);


		// LPSR
		NP = (size_t)round(((1.0 * (int)(Nmin - NPinit) / (double)Max_NFE) * nfe) + NPinit);
		A = (size_t)NP * rarc_;

		for (size_t i = 0; i < sub_num; i++)
		{
			sub_pops[i].pop().sort();
			sub_pops[i].pop().resize(NP);
			if (sub_pops[i].archive_pop().size() > A)
			{
				sub_pops[i].archive_pop().shuffle();
				sub_pops[i].archive_pop().resize(A);
			}
		}
	}

	for (size_t i = 0; i < sub_num; i += 1)
	{
		sub_pops[i].pop().sort();
	}

	size_t best_sub_num = 0;
	for (size_t i = 1; i < sub_num; i += 1)
	{
		if (sub_pops[i].pop()[0].fitness() < sub_pops[best_sub_num].pop()[0].fitness())
			best_sub_num = i;
	}


	std::cout << "Cycle_cnt: " << cycle_cnt << " ";
	
	*solutions = sub_pops[best_sub_num].pop();
	log.close();
}

Individual::GeneVec mpmL_SHADE::CurtopBest_DonorVec(int target_idx, double p, double f, const Population& pop, const Population& archive)
{
	size_t Pop_Size = pop.size();
	int best_idx_max = round(Pop_Size * p);
	if (best_idx_max < 2) best_idx_max = 2;

	int x_best = alg_math::randInt(0, best_idx_max - 1),
		x_r1 = alg_math::randInt(0, (int)(Pop_Size - 1)),
		x_r2 = alg_math::randInt(0, (int)(Pop_Size + archive.size() - 1));
	while (x_r1 == target_idx || x_r2 == target_idx || x_r1 == x_r2)
	{
		x_r1 = rand() % Pop_Size;
		x_r2 = rand() % (Pop_Size + archive.size());
	}
	size_t gene_len = pop[target_idx].gene().size();
	Individual::GeneVec g_cur = pop[target_idx].gene(),
		g_best = pop[x_best].gene(),
		g_r1 = pop[x_r1].gene(),
		g_r2 = (x_r2 >= Pop_Size) ? archive[x_r2 - Pop_Size].gene() : pop[x_r2].gene();

	Individual::GeneVec donor_vector(gene_len);
	for (size_t j = 0; j < gene_len; ++j)
	{
		donor_vector[j] = g_cur[j] + f * (g_best[j] - g_cur[j] + g_r1[j] - g_r2[j]);
	}
	return donor_vector;
}