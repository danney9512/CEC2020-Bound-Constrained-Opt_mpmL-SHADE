
#include "experiment.h"
#include <string>

Experiment::Experiment(std::string algo_str, size_t pid, int dim)
{
	algo_name_ = algo_str;
	problem_id_ = pid;
	dim_ = dim;
}

ExperimentSet::ExperimentSet(std::ifstream &ifile)
{
	std::string algo_str;
	size_t pid;

	int dim;

	std::string dummy;
	while (ifile >> dummy)
	{
		ifile >> dummy >> dummy >> pid;
		ifile >> dummy >> dummy >> algo_str;
		ifile >> dummy >> dummy >> dim;
		exps.push_back(Experiment(algo_str, pid, dim));
	}
}