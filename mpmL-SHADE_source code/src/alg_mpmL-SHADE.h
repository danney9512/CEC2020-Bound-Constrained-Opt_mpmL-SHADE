#ifndef ALG_mpmLSHADE__
#define ALG_mpmLSHADE__

#include "alg_individual.h"
#include "alg_population.h"
#include "alg_base.h"
#include <cstddef>
#include <string>
#include <fstream>
#include <vector>

class mpmL_SHADE : public BaseEA
{
public:
	mpmL_SHADE() : BaseEA("mpmL-SHADE") {}

	virtual void Setup(std::ifstream& ifile);
	virtual void Solve(Population* solutions, const CProblem& prob);


	class Memory
	{
	public:
		explicit Memory(double f = 0.0, double cr = 0.0) : f_(f), cr_(cr) {}

		double& F() { return f_; }
		const double& F() const { return f_; }
		double& CR() { return cr_; }
		const double& CR() const { return cr_; }

	private:
		double f_,
			cr_;
	};
	class MemorySystem
	{
	public:
		explicit MemorySystem(std::size_t H = 0) : memories_(H), k_(0), h_(H) {}

		size_t& K() { return k_; }
		const size_t& K() const { return k_; }
		size_t& H() { return h_; }
		const size_t& H() const { return h_; }
		Memory& operator[](std::size_t i) { return memories_[i]; }
		const Memory& operator[](std::size_t i) const { return memories_[i]; }

		std::size_t size() const { return memories_.size(); }
		bool empty() const { return size() == 0; }
		void push_back(const Memory& mem)
		{
			memories_.push_back(mem);
			h_ += 1;
		}
		void resize(std::size_t t) { memories_.resize(t); }
		void update_memory(const std::vector<Memory>& success_parameter, const std::vector<double>& success_fit_dif);
		void pertub_memory();
	private:
		size_t k_, h_;
		std::vector<Memory> memories_;
	};

	class SubPopulation
	{
	public:
		typedef std::vector<double> DoubleVec;
		typedef std::vector<Memory> MemoryVec;

		explicit SubPopulation(std::size_t NP = 0, std::size_t H = 0) : pop_(NP), children_(NP), memory_sys_(H), stop_flag_(false), stop_idx_(0), no_success_cnt_(0), max_no_success_cnt_(0), perturb_cnt_(0), total_no_success_(0), fitness_stop_nfe_(0), past_best_fitness_(0), past_best_stop_nfe_(0){}

		Population& pop() { return pop_; }
		const Population& pop() const { return pop_; }
		Population& archive_pop() { return archive_pop_; }
		const Population& archive_pop() const { return archive_pop_; }
		Population& children() { return children_; }
		const Population& children() const { return children_; }

		MemorySystem& memory_sys() { return memory_sys_; }
		const MemorySystem& memory_sys() const { return memory_sys_; }

		MemoryVec& sample_parameter() { return sample_parameter_; }
		const MemoryVec& sample_parameter() const { return sample_parameter_; }
		MemoryVec& success_parameter() { return success_parameter_; }
		const MemoryVec& success_parameter() const { return success_parameter_; }
		DoubleVec& success_fit_dif() { return success_fit_dif_; }
		const DoubleVec success_fit_dif()const { return success_fit_dif_; }

		bool& stop_flag() { return stop_flag_; }
		const bool& stop_flag() const { return stop_flag_; }

		size_t& stop_idx() { return stop_idx_; }
		const size_t& stop_idx() const { return stop_idx_; }
		size_t& no_success_cnt() { return no_success_cnt_; }
		const size_t& no_success_cnt() const { return no_success_cnt_; }
		size_t& max_no_success_cnt() { return max_no_success_cnt_; }
		const size_t& max_no_success_cnt() const { return max_no_success_cnt_; }
		size_t& perturb_cnt() { return perturb_cnt_; }
		const size_t& perturb_cnt() const { return perturb_cnt_; }
		size_t& total_no_success() { return total_no_success_; }
		const size_t& total_no_success() const { return total_no_success_; }

		unsigned long long int& fitness_stop_nfe() { return fitness_stop_nfe_; }
		const unsigned long long int& fitness_stop_nfe() const { return fitness_stop_nfe_; }
		double& past_best_fitness() { return past_best_fitness_; }
		const double& past_best_fitness() const { return past_best_fitness_; }
		unsigned long long int& past_best_stop_nfe() { return past_best_stop_nfe_; }
		const unsigned long long int& past_best_stop_nfe() const { return past_best_stop_nfe_; }

	private:
		Population pop_;
		Population archive_pop_;
		Population children_;
		
		MemorySystem memory_sys_;
		MemoryVec sample_parameter_, success_parameter_;
		DoubleVec success_fit_dif_;

		// Set variables for terminal criteria
		bool stop_flag_;
		size_t stop_idx_, no_success_cnt_, max_no_success_cnt_, perturb_cnt_, total_no_success_;

		//fitness stop parameter
		unsigned long long int fitness_stop_nfe_;
		unsigned long long int past_best_stop_nfe_;
		double past_best_fitness_;
	};

	void SubPopGeneration(SubPopulation& sub_pop, const CProblem& prob, const std::size_t& NP, const std::size_t& A, std::default_random_engine& generator, unsigned long long int& nfe);

private:
	double p_,
		rarc_,
		finit_,
		crinit_,
		scalemax_;
	size_t h_,
		rNinit_,
		nmin_;
	unsigned long long int max_nfe_;
	Individual::GeneVec CurtopBest_DonorVec(int target_idx, double p, double f, const Population& pop, const Population& archive);
};
#endif // !ALG_mpmLSHADE__


