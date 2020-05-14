#ifndef PROBLEM__
#define PROBLEM__

#include <iostream>
#include <ostream>
#include <string>
#include <vector>

class CProblem
{
public:
	typedef std::vector<double> RNV;
	explicit CProblem(int id = -1, std::string name = "", int dim = 0, double l_bd = 0.0, double u_bd = 0.0, double g_opm = 0.0, bool s_flag = false, bool r_flag = false, bool sh_flag = false, bool cp_flag = false);

	const int & id() const { return id_; }
	const std::string & name() const { return func_name_; }
	const int & dim() const { return dimension_; }

	const double & lower_bound() const { return lower_bound_; }
	const double & upper_bound() const { return upper_bound_; }
	const double & global_optimum() const { return global_optimum_; }

	const bool & shiftable() const { return shift_flag_; }
	const bool & rotatable() const { return rotate_flag_; }
	const bool & shuffleable() const { return shuffle_flag_; }
	const bool & composition() const { return comp_flag_; }

	const std::vector<RNV> & rotate_matrix() const { return rotate_matrix_; }
	const RNV & shift_val_vec() const { return shift_val_vec_; }
	const std::vector<int> & shuffle_val_vec() const { return shuffle_val_vec_; }
	
	void set_id(const size_t new_id) { id_ = new_id; }
	void set_name(const std::string new_name) { func_name_ = new_name; }
	void set_dim(const size_t new_dim) { dimension_ = new_dim; }

	void set_lower_bound(const double new_lb) { lower_bound_ = new_lb; }
	void set_upper_bound(const double new_ub) { upper_bound_ = new_ub; }
	void set_global_optimum(const double new_go) { global_optimum_ = new_go; }

	void set_shift_flag(const bool new_sf) { shift_flag_ = new_sf; }
	void set_rotate_flag(const bool new_rf) { rotate_flag_ = new_rf; }
	void set_shuffle_flag(const bool new_shf) { shuffle_flag_ = new_shf; }
	void set_comp_flag(const bool new_cf) { comp_flag_ = new_cf; }

	void set_rotate_matrix(const std::vector<RNV> & new_rm) { rotate_matrix_ = new_rm; }
	void set_shift_vector(const RNV & new_sv) { shift_val_vec_ = new_sv; }
	void set_shuffle_vector(const std::vector<int> & new_shv) { shuffle_val_vec_ = new_shv; }

private:
	int id_, dimension_;
	std::string func_name_;
	double global_optimum_, lower_bound_, upper_bound_;
	bool shift_flag_, rotate_flag_, shuffle_flag_;
	bool comp_flag_;

	std::vector<RNV> rotate_matrix_;
	RNV shift_val_vec_;
	std::vector<int> shuffle_val_vec_;
	
};

std::ostream & operator << (std::ostream &os, const CProblem &problem);

#endif // !PROBLEM__
