#ifndef EXPERIMENT__
#define EXPERIMENT__

#include <fstream>
#include "alg_base.h"
#include <vector>

class Experiment
{
public:
	explicit Experiment(std::string algo_str = "", size_t pid = 0, int dim = 0);
	const std::string & algo_name() { return algo_name_; }
	const size_t & problem_id() { return problem_id_; }
	const int & dim() { return dim_; }

private:
	std::string algo_name_;
	size_t problem_id_;
	int dim_;
};

class ExperimentSet
{
public:
	explicit ExperimentSet(std::ifstream &ifile);

	Experiment & operator[](std::size_t i) { return exps[i]; }
	const Experiment & operator[](std::size_t i) const { return exps[i]; }

	std::size_t size() const { return exps.size(); }
	bool empty() const { return size() == 0; }
	void resize(std::size_t t) { exps.resize(t); }
	void push_back(const Experiment &exp) { exps.push_back(exp); }
	void clear() { exps.clear(); }
private:
	std::vector<Experiment> exps;
};


#endif // !EXPERIMENT__
