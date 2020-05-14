//=========================================================================================
// Multi-population Modified L-SHADE for Single Objective Bound Constrained Optimization
// mpmL-SHADE
//=========================================================================================
// Authors     : Yann-Chern Jou
//               Shuo-Ying Wang
//               Jia-Fong Yeh
//				 Tsung-Che Chiang
// Version     : v1.2
// Created on  : Mar 13, 2020
// Language	   : C++ (std C++17)
//
// More details on the following paper:
//
// CEC2020
// "Multi-population Modified L-SHADE for Single Objective Bound Constrained Optimization"
//  
//=========================================================================================

#include "problem.h"
#include "problem_set.h"
#include "problem_test_functions.h"
#include "alg_population.h"
#include "alg_base.h"
#include "alg_L-SHADE.h"
#include "alg_mL-SHADE.h"
#include "alg_mpmL-SHADE.h"
#include "alg_log.h"
#include "alg_result.h"
#include "experiment.h"
#include "self_implement_functions.h"


#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <string>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <vector>

using namespace std;

const int problem_D (const int dim)
{
	return (dim / 5) - 1;
}


int main()
{
	// Set random seed
	srand(time(0));

	// Compute time initialize
	double START, END;

	//Load problem sets with 4 different dimemsions (D5, D10, D15, D20)
	constexpr int all_dim = 4;
	vector<CPromblemSet> prob_sets;
	for (int i = 0; i < all_dim; i += 1)
	{
		string problem_filename = "problem_list_D" + IntToStr((i + 1) * 5, 2) + ".ini";

		ifstream prob_set_ini(problem_filename);
		if (!prob_set_ini) cout << "Error msg: Can't load the problems' definition." << endl;

		CPromblemSet prob_set(prob_set_ini);
		prob_sets.push_back(prob_set);
	}

	// set experiment file names
	string experiment_filename = "experiment_list.ini";

	// Load experiment's setting
	ifstream exp_list_ini(experiment_filename);
	if (!exp_list_ini) cout << "Error msg: Can't load the experiments." << endl;
	ExperimentSet exp_set(exp_list_ini);


	BaseEA* ea = nullptr;
	string algo_name, algo_ini_str;
	size_t problem_id;
	int problem_dim;
	constexpr size_t NumOfRuns = 30;

	cout << fixed << setprecision(15) << "Experiment Started." << endl;
	cout << "If you get any error message above this line, terminate the program." << endl << endl;

	// Exp start
	START = clock();

	for (size_t i = 0; i < exp_set.size(); ++i)
	{
		algo_name = exp_set[i].algo_name();
		problem_id = exp_set[i].problem_id();
		problem_dim = exp_set[i].dim();

		cout << "Experiment in F" << IntToStr(problem_id,2)  << " D" << IntToStr(problem_dim, 2) << ", " << algo_name << " solves " << prob_sets[problem_D(problem_dim)][problem_id - 1].name()  << endl;

		// ********************************************* //
		// You can add your own algorithm name here ! //
		// ********************************************* //
		if (algo_name == "L-SHADE")
		{
			ea = new L_SHADE();
		}
		else if (algo_name == "mL-SHADE")
		{
			ea = new mL_SHADE();
		}
		else if (algo_name == "mpmL-SHADE")
		{
			ea = new mpmL_SHADE();
		}
		else // Skip this experiment if we can't recognize the alorithm
		{
			cout << "Error msg: Can't recognize the algorithm." << endl;
			continue;
		}

		// Loading algorithm parameter's setting
		algo_ini_str = "experiments\\" + algo_name + "\\exp_CEC2020_" + IntToStr(problem_id, 2) + "_D" + IntToStr(problem_dim, 2) + ".ini";
		ifstream exp_algo_para_ini(algo_ini_str);

		if (!exp_algo_para_ini)
		{
			cout << "Error msg: Can't load the algorithm's parameters, check the .ini file's path." << endl;
		}
		else
		{
			// setting algorithm parameter
			ea->Setup(exp_algo_para_ini);
			
			string result_filename = algo_name + "_F" + IntToStr(problem_id, 2) + "_D" + IntToStr(problem_dim, 2) + ".csv";
			Result result(problem_id, NumOfRuns, result_filename);

			// Running EA
			for (size_t k = 0; k < NumOfRuns; ++k)
			{
				Log::set_k((int)k);
				Population solutions;
				ea->Solve(&solutions, prob_sets[problem_D(problem_dim)][problem_id - 1]);

				result[k] = solutions[0].fitness();
				cout << "Run " << k << ", best fitness: " << solutions[0].fitness() <<  endl;
			}

			result.compute(prob_sets[problem_D(problem_dim)][problem_id - 1].global_optimum());
			result.OutputToFile(prob_sets[problem_D(problem_dim)][problem_id - 1].global_optimum());
			result.close();
			cout << "Experiment " << i + 1 << " ended." << endl << endl;
			

			if (ea != nullptr) delete ea;
		}

	}

	END = clock();
	cout << "All Experiment for is ended." << endl;
	cout << "Execution Time: " << (END - START) / CLOCKS_PER_SEC << " sec" << endl << endl;


	system("pause");
	return 0;
}