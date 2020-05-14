#include "problem_set.h"
#include "problem.h"

#include <string>

CPromblemSet::CPromblemSet(std::ifstream &ifile)
{
	std::string name;
	int ID, dim;
	double lower_bound, upper_bound, global_optimum;
	bool s_flag, r_flag, sh_flag, cp_flag;

	std::string dummy; 
	while (ifile >> dummy)
	{
		ifile >> dummy >> dummy >> ID;
		ifile >> dummy >> dummy >> name;
		ifile >> dummy >> dummy >> dim;
		ifile >> dummy >> dummy >> lower_bound >> upper_bound;
		ifile >> dummy >> dummy >> global_optimum;
		ifile >> dummy >> dummy >> s_flag;
		ifile >> dummy >> dummy >> r_flag;
		ifile >> dummy >> dummy >> sh_flag;
		ifile >> dummy >> dummy >> cp_flag;
		problems_.push_back(CProblem(ID, name, dim, lower_bound, upper_bound, global_optimum, s_flag, r_flag, sh_flag, cp_flag));
	}
}