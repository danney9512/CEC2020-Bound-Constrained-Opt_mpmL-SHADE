#include "problem.h"
#include "self_implement_functions.h"
#include <sstream>
#include <fstream>
#include <cstdio>

CProblem::CProblem(int id, std::string name, int dim, double l_bd, double u_bd, double g_opm, bool s_flag, bool r_flag, bool sh_flag, bool cp_flag)
{
	id_ = id;
	func_name_ = name;
	dimension_ = dim;
	lower_bound_ = l_bd;
	upper_bound_ = u_bd;
	global_optimum_ = g_opm;
	shift_flag_ = s_flag;
	rotate_flag_ = r_flag;
	shuffle_flag_ = sh_flag;
	comp_flag_ = cp_flag;

	size_t times = 1;
	if (comp_flag_) times = 10;

	if (shift_flag_)
	{
		std::string file_name = "input_data\\shift_data_", dump;
		file_name += IntToStr(id_, 2) + ".txt";

		double val;
		std::ifstream ip(file_name);
		if (ip)
		{
			for (int k = 0; k < times; ++k)
			{
				for (int i = 0; i < dimension_; ++i)
				{
					ip >> val;
					shift_val_vec_.push_back(val);
				}
				std::getline(ip, dump);
			}
		}
		else
		{
			printf("F%d Can't load shift data!\n", id_);
		}
		ip.close();
	}
	if (rotate_flag_)
	{
		std::string file_name = "input_data\\M_" + IntToStr(id_, 2) + "_D" + IntToStr(dimension_, 2) + ".txt";
		double val;

		std::ifstream ip(file_name);
		if (ip)
		{
			for (int i = 0; i < times * dimension_; ++i)
			{
				RNV tmp_r_vec;
				for (int j = 0; j < dimension_; ++j)
				{
					ip >> val;
					tmp_r_vec.push_back(val);
				}
				rotate_matrix_.push_back(tmp_r_vec);
			}
		}
		else
		{
			printf("F%d Can't load rotate matrix!\n", id_);
		}
		ip.close();
	}

	if (shuffle_flag_)
	{	
		std::string file_name = "input_data\\shuffle_data_" + IntToStr(id_, 2) + "_D" + IntToStr(dimension_, 2) + ".txt";
		int val;

		std::ifstream ip(file_name);
		if (ip)
		{
			for (int i = 0; i < times * dimension_; ++i)
			{
				ip >> val;
				shuffle_val_vec_.push_back(val);
			}
		}
		else
		{
			printf("F%d Can't load shuffle vector!\n", id_);
		}
		ip.close();
	}

}

std::ostream & operator << (std::ostream &os, const CProblem &problem)
{
	os << "Problem's Info :" << std::endl;
	os << "Id : " << problem.id() << " / Name : " << problem.name() << std::endl;
	os << "Dim : " << problem.dim() << std::endl;
	os << "Gloabal Optimum : " << problem.global_optimum() << std::endl;
	os << "Search Range : [" << problem.lower_bound() << "," << problem.upper_bound() << "]" << std::endl;
	os << "Shiftable : " << problem.shiftable() << " / Rotatable : " << problem.rotatable() << std::endl;
	return os;
}