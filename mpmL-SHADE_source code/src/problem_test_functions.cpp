#include "problem_test_functions.h"
#include <cmath>
#include <cstdio>

double Storn_Chebyshev(const Individual::GeneVec &gene_vec)
{
	/* Valid for any D>2
	constraints: unconstrained
	type: multi-modal with one global minimum; non-separable
	initial upper bound = 2^D, initial lower bound = -D^n
	value-to-reach = f(x*)+1.0e-8
	f(x*)=0.0; x*=(128,0,-256,0,160,0,-32,0,1) (n=9)
	x*=(32768,0,-131072,0,212992,0,-180224,0,84480,0,-21504,0,2688,0,-128,0,1) (n=17)
	*/
	
	static int sample;
	static double dx, dy;
	double a = 1.0, b = 1.2, px, y = -1, sum = 0.0;

	for (int i = 0; i < gene_vec.size()-2; ++i)
	{
		dx = 2.4 * b - a;
		a = b; 
		b = dx;
	}

	sample =  32 * gene_vec.size();
	dy = 2.0 / (double)sample;

	for (int i = 0; i <= sample; ++i)
	{
		px = gene_vec[0];
		for (int j = 1; j < gene_vec.size(); ++j)
		{
			px = y*px + gene_vec[j];
		}
		if (px < -1 || px > 1) sum += (1. - fabs(px))*(1. - fabs(px));
		y += dy;
	}
	for (int i = -1; i <= 1; i += 2)
	{
		px = gene_vec[0];
		for (int j = 1; j <  gene_vec.size(); ++j)
		{
			px = 1.2*px + gene_vec[j];
		}
		if (px < dx) sum += px * px;
	}
	return sum;
}

double Inverse_Hilbert(const Individual::GeneVec &gene_vec)
{
	/* valid for any dimension, n=k*k, k=2,3,4,...
	constraints: unconstrained
	type: multi-modal with one global minimum; non-separable
	initial upper bound = 2^n, initial lower bound = -2^n
	value-to-reach = f(x*)+1.0e-8
	f(x*) = 0.0; x*={{9,-36,30},{-36,192,-180},{30,-180,180}} (n=9)
	x*={{16,-120,240,-140},{-120,1200,-2700,1680},{240,-2700,6480,4200},{-140,1680,-4200,2800}} (n=16)
	*/
	int b = (int)( sqrt( (double)(gene_vec.size()) ) );
	double sum = 0.0;
	static double hilbert[10][10], y[10][10];

	for (int i = 0; i < b; ++i)
	{
		for (int j = 0; j < b; ++j)
		{
			hilbert[i][j] = 1. / (double)(i + j + 1);
		}
	}
	for (int j = 0; j < b; ++j)
	{
		for (int k = 0; k < b; ++k)
		{
			y[j][k] = 0;
			for (int i = 0; i < b; ++i)
			{
				y[j][k] += hilbert[j][i] * gene_vec[k + b * i];
			}
		}
	}
	for (int i = 0; i < b; ++i)
	{
		for (int j = 0; j < b; ++j)
		{
			if (i == j) sum += fabs(y[i][j] - 1);
			else sum += fabs(y[i][j]);
		}
	}
	return sum;
}

double Lennard_Jones(const Individual::GeneVec &gene_vec)
{
	/* valid for any dimension, D=3*k, k=2,3,4,...,25.   k is the number of atoms in 3-D space
	constraints: unconstrained
	type: multi-modal with one global minimum; non-separable
	initial upper bound = 4, initial lower bound = -4
	value-to-reach = minima[k-2]+.0001
	f(x*) = minima[k-2]; see array of minima below; additional minima available at the
	Cambridge cluster database: http://www-wales.ch.cam.ac.uk/~jon/structures/LJ/tables.150.html
	*/
	int D = (int) gene_vec.size(), a, b;
	double xd, yd, zd, ed, ud, sum = 0.0;

	static double minima[] = { -1.,-3.,-6.,-9.103852,-12.712062,-16.505384,-19.821489,-24.113360,
		-28.422532,-32.765970,-37.967600,-44.326801,-47.845157,-52.322627,-56.815742,-61.317995,
		-66.530949,-72.659782,-77.1777043,-81.684571,-86.809782,-02.844472,-97.348815,-102.372663 };

	int k = D / 3;
	if (k < 2)
	{
		k = 2;
		D = 6;
	}

	for (int i = 0; i < k - 1; ++i)
	{
		for (int j = i + 1; j < k; ++j)
		{
			a = 3 * i;
			b = 3 * j;
			xd = gene_vec[a] - gene_vec[b];
			yd = gene_vec[a + 1] - gene_vec[b + 1];
			zd = gene_vec[a + 2] - gene_vec[b + 2];
			ed = xd*xd + yd*yd + zd*zd;
			ud = ed*ed*ed;
			if (ed>0) sum += (1.0 / ud - 2.0) / ud;
			else sum += 1.0e20;
		}
	}
	return sum + 12.7120622568;
}

double Rastrihin(const Individual::GeneVec &gene_vec, const CProblem &prob)
{
	Individual::GeneVec sr_gene_vec(gene_vec);
	sr_gene_vec = sr_func(sr_gene_vec, prob.shift_val_vec(), prob.rotate_matrix(), 5.12 / 100.0, prob.shiftable(), prob.rotatable());

	double sum = 0.0;
	for (int i = 0; i<sr_gene_vec.size(); ++i)
	{
		sum += (sr_gene_vec[i] * sr_gene_vec[i] - 10.0*cos(2.0*PI*sr_gene_vec[i]) + 10.0);
	}
	return sum;
}

double Griewangk(const Individual::GeneVec &gene_vec, const CProblem &prob)
{
	Individual::GeneVec sr_gene_vec(gene_vec);
	sr_gene_vec = sr_func(sr_gene_vec, prob.shift_val_vec(), prob.rotate_matrix(), 600.0 / 100.0, prob.shiftable(), prob.rotatable());

	double s = 0.0, p = 1.0, sum = 0.0;

	for (int i = 0; i < sr_gene_vec.size(); ++i)
	{
		s += (sr_gene_vec[i] * sr_gene_vec[i]);
		p *= cos(sr_gene_vec[i] / sqrt(1.0 + i));
	}
	sum = 1.0 + s / 4000.0 - p;
	return sum;
}

double Weierstrass(const Individual::GeneVec &gene_vec, const CProblem &prob)
{
	Individual::GeneVec sr_gene_vec(gene_vec);
	sr_gene_vec = sr_func(sr_gene_vec, prob.shift_val_vec(), prob.rotate_matrix(), 0.5 / 100.0, prob.shiftable(), prob.rotatable());

	int k_max = 20;
	double sum = 0.0, sum1, sum2, a = 0.5, b = 3.0, two_times_pi = 2.0*PI;

	std::vector<double> a_pow_vec, b_pow_vec, sum2_vec, two_times_pi_times_b_vec;
	for (int i = 0; i <= k_max; ++i)
	{
		a_pow_vec.push_back(pow(a, i));
		b_pow_vec.push_back(pow(b, i));
		sum2_vec.push_back(a_pow_vec[i] * cos(2.0*PI*b_pow_vec[i] * 0.5));
		two_times_pi_times_b_vec.push_back(two_times_pi*b_pow_vec[i]);
	}


	for (int i = 0; i < sr_gene_vec.size(); ++i)
	{
		sum1 = 0.0;
		sum2 = 0.0;
		for (int j = 0; j <= k_max; ++j)
		{
			sum1 += a_pow_vec[j]*cos(two_times_pi_times_b_vec[j]*(sr_gene_vec[i] + 0.5));
			sum2 += sum2_vec[j];//pow(a, j)*cos(2.0*PI*pow(b, j)*0.5);
		}
		sum += sum1;
	}
	return sum - sr_gene_vec.size() * sum2;
}

double Modified_Schwefel(const Individual::GeneVec &gene_vec, const CProblem &prob)
{
	Individual::GeneVec sr_gene_vec(gene_vec);
	sr_gene_vec = sr_func(sr_gene_vec, prob.shift_val_vec(), prob.rotate_matrix(), 1000.0 / 100.0, prob.shiftable(), prob.rotatable());

	double sum = 0.0, tmp;
	for (int i = 0; i < sr_gene_vec.size(); ++i)
	{
		sr_gene_vec[i] += 4.209687462275036e+002;
		if (sr_gene_vec[i] > 500)
		{
			sum -= (500.0 - fmod(sr_gene_vec[i], 500))*sin(pow(500.0 - fmod(sr_gene_vec[i], 500), 0.5));
			tmp = (sr_gene_vec[i] - 500.0) / 100;
			sum += tmp*tmp / sr_gene_vec.size();
		}
		else if (sr_gene_vec[i] < -500)
		{
			sum -= (-500.0 + fmod(fabs(sr_gene_vec[i]), 500))*sin(pow(500.0 - fmod(fabs(sr_gene_vec[i]), 500), 0.5));
			tmp = (sr_gene_vec[i] + 500.0) / 100;
			sum += tmp*tmp / sr_gene_vec.size();
		}
		else
		{
			sum -= sr_gene_vec[i] * sin(pow(fabs(sr_gene_vec[i]), 0.5));
		}
	}
	sum += 4.189828872724338e+002 * sr_gene_vec.size();
	return sum;
}

double Expanded_Schaffer_F6(const Individual::GeneVec &gene_vec, const CProblem &prob)
{
	Individual::GeneVec sr_gene_vec(gene_vec);
	sr_gene_vec = sr_func(sr_gene_vec, prob.shift_val_vec(), prob.rotate_matrix(), 1.0, prob.shiftable(), prob.rotatable());
	double temp1, temp2, sum = 0.0;
	for (int i = 0; i < (int)sr_gene_vec.size() - 1; ++i)
	{	
		temp1 = sin(sqrt(sr_gene_vec[i] * sr_gene_vec[i] + sr_gene_vec[i + 1] * sr_gene_vec[i + 1]));
		temp1 = temp1*temp1;
		temp2 = 1.0 + 0.001*(sr_gene_vec[i] * sr_gene_vec[i] + sr_gene_vec[i + 1] * sr_gene_vec[i + 1]);
		sum += 0.5 + (temp1 - 0.5) / (temp2*temp2);
	}
	
	temp1 = sin(sqrt(sr_gene_vec[sr_gene_vec.size() - 1] * sr_gene_vec[sr_gene_vec.size() - 1] + sr_gene_vec[0] * sr_gene_vec[0]));
	temp1 = temp1*temp1;
	temp2 = 1.0 + 0.001*(sr_gene_vec[sr_gene_vec.size() - 1] * sr_gene_vec[sr_gene_vec.size() - 1] + sr_gene_vec[0] * sr_gene_vec[0]);
	sum += 0.5 + (temp1 - 0.5) / (temp2*temp2);
	
	return sum;
}

double Happy_Cat(const Individual::GeneVec &gene_vec, const CProblem &prob)
{
	Individual::GeneVec sr_gene_vec(gene_vec);
	sr_gene_vec = sr_func(sr_gene_vec, prob.shift_val_vec(), prob.rotate_matrix(), 5.0 / 100.0, prob.shiftable(), prob.rotatable());

	double alpha = 1.0 / 8.0, r2 = 0.0, sum_z = 0.0, sum = 0.0;

	for (int i = 0; i < sr_gene_vec.size(); ++i)
	{
		sr_gene_vec[i] = sr_gene_vec[i] - 1.0;
		r2 += sr_gene_vec[i] * sr_gene_vec[i];
		sum_z += sr_gene_vec[i];
	}
	sum = pow(fabs(r2 - sr_gene_vec.size()), 2 * alpha) + (0.5*r2 + sum_z) / sr_gene_vec.size() + 0.5;
	return sum;
}

double Ackley(const Individual::GeneVec &gene_vec, const CProblem &prob)
{
	Individual::GeneVec sr_gene_vec(gene_vec);
	sr_gene_vec = sr_func(sr_gene_vec, prob.shift_val_vec(), prob.rotate_matrix(), 1.0, prob.shiftable(), prob.rotatable());
	
	double sum1 = 0.0, sum2 = 0.0, sum = 0.0;
	for (int i = 0; i < sr_gene_vec.size(); ++i)
	{
		sum1 += sr_gene_vec[i] * sr_gene_vec[i];
		sum2 += cos(2.0*PI*sr_gene_vec[i]);
	}
	sum1 = -0.2*sqrt(sum1 / sr_gene_vec.size());
	sum2 /= sr_gene_vec.size();
	sum = E - 20.0*exp(sum1) - exp(sum2) + 20.0;
	return sum;
}

Individual::GeneVec shiftfunc(const Individual::GeneVec &gene_vec, const CProblem::RNV &shift_vec)
{
	if (gene_vec.size() != shift_vec.size()) return Individual::GeneVec(gene_vec);

	Individual::GeneVec shifted_gene_vec(gene_vec.size());
	for (int i = 0; i < gene_vec.size(); ++i)
	{
		shifted_gene_vec[i] = (gene_vec[i] - shift_vec[i]);
	}
	return shifted_gene_vec;
}

Individual::GeneVec rotatefunc(const Individual::GeneVec &gene_vec, const std::vector<CProblem::RNV> &rotate_matrix)
{
	if( !(gene_vec.size() == rotate_matrix.size() && gene_vec.size())) return Individual::GeneVec(gene_vec);

	Individual::GeneVec rotated_gene_vec(gene_vec.size());

	for (int i = 0; i < gene_vec.size(); ++i)
	{
		rotated_gene_vec[i] = 0.0;
		for (int j = 0; j < gene_vec.size(); ++j)
		{
			rotated_gene_vec[i] += gene_vec[j] * rotate_matrix[i][j];
		}
	}
	return rotated_gene_vec;
}

Individual::GeneVec sr_func(const Individual::GeneVec &gene_vec, const CProblem::RNV &shift_vec, const std::vector<CProblem::RNV> &rotate_matrix, const double& shrink_rate, const bool& shift_flag, const bool& rotate_flag)
{
	Individual::GeneVec sr_gene_vec(gene_vec);
	if (shift_flag)
	{
		sr_gene_vec = shiftfunc(sr_gene_vec, shift_vec);
	}

	for (int i = 0; i < sr_gene_vec.size(); ++i)
	{
		sr_gene_vec[i] *= shrink_rate;
	}

	if (rotate_flag)
	{
		sr_gene_vec = rotatefunc(sr_gene_vec, rotate_matrix);
	}
	return sr_gene_vec;
}

Individual::GeneVec shuffle_func(const Individual::GeneVec &gene_vec, const std::vector<int> &shuffle_vec)
{
	Individual::GeneVec shuffle_gene_vec(gene_vec);
	for (int i = 0; i < gene_vec.size(); ++i)
	{
		shuffle_gene_vec[i] = gene_vec[shuffle_vec[i] - 1];
	}
	return shuffle_gene_vec;
}


void CEC2019_100Digit_Evaluate(Individual *indiv, const CProblem &prob)
{
	if (prob.name() == "Storn_Chebyshev")              indiv->fitness() = Storn_Chebyshev(indiv->gene()) + 1.0;
	else if (prob.name() == "Inverse_Hilbert")         indiv->fitness() = Inverse_Hilbert(indiv->gene()) + 1.0;
	else if (prob.name() == "Lennard_Jones")           indiv->fitness() = Lennard_Jones(indiv->gene()) + 1.0;
	else if (prob.name() == "Rastrihin")               indiv->fitness() = Rastrihin(indiv->gene(), prob) + 1.0;
	else if (prob.name() == "Griewangk")               indiv->fitness() = Griewangk(indiv->gene(), prob) + 1.0;
	else if (prob.name() == "Weierstrass")             indiv->fitness() = Weierstrass(indiv->gene(), prob) + 1.0;
	else if (prob.name() == "Modified_Schwefel")       indiv->fitness() = Modified_Schwefel(indiv->gene(), prob) + 1.0;
	else if (prob.name() == "Expanded_Schaffer_F6")    indiv->fitness() = Expanded_Schaffer_F6(indiv->gene(), prob) + 1.0;
	else if (prob.name() == "Happy_Cat")               indiv->fitness() = Happy_Cat(indiv->gene(), prob) + 1.0;
	else if (prob.name() == "Ackley")                  indiv->fitness() = Ackley(indiv->gene(), prob) + 1.0;
	else                                               printf("Can't find the corresponding function!\n");
}

void CEC2014_Evaluate(Individual *indiv, const CProblem &prob)
{
	switch (prob.id())
	{
	case 1:
		indiv->fitness() = Ellips(indiv->gene(), prob) + 100.0;
		break;
	case 2:
		indiv->fitness() = Bent_Cigar(indiv->gene(), prob) + 200.0;
		break;
	case 3:
		indiv->fitness() = Discus(indiv->gene(), prob) + 300.0;
		break;
	case 4:
		indiv->fitness() = Rosenbrock(indiv->gene(), prob) + 400.0;
		break;
	case 5:
		indiv->fitness() = Ackley(indiv->gene(), prob) + 500.0;
		break;
	case 6:
		indiv->fitness() = Weierstrass(indiv->gene(), prob) + 600.0;
		break;
	case 7:
		indiv->fitness() = Griewangk(indiv->gene(), prob) + 700.0;
		break;
	case 8:
		indiv->fitness() = Rastrihin(indiv->gene(), prob) + 800.0;
		break;
	case 9:
		indiv->fitness() = Rastrihin(indiv->gene(), prob) + 900.0;
		break;
	case 10:
		indiv->fitness() = Modified_Schwefel(indiv->gene(), prob) + 1000.0;
		break;
	case 11:
		indiv->fitness() = Modified_Schwefel(indiv->gene(), prob) + 1100.0;
		break;
	case 12:
		indiv->fitness() = Katsuura(indiv->gene(), prob) + 1200.0;
		break;
	case 13:
		indiv->fitness() = Happy_Cat(indiv->gene(), prob) + 1300.0;
		break;
	case 14:
		indiv->fitness() = Hgbat(indiv->gene(), prob) + 1400.0;
		break;
	case 15:
		indiv->fitness() = Grie_Rosen(indiv->gene(), prob) + 1500.0;
		break;
	case 16:
		indiv->fitness() = Expanded_Schaffer_F6(indiv->gene(), prob) + 1600.0;
		break;
	case 17:
		indiv->fitness() = hf01(indiv->gene(), prob) + 1700.0;
		break;
	case 18:
		indiv->fitness() = hf02(indiv->gene(), prob) + 1800.0;
		break;
	case 19:
		indiv->fitness() = hf03(indiv->gene(), prob) + 1900.0;
		break;
	case 20:
		indiv->fitness() = hf04(indiv->gene(), prob) + 2000.0;
		break;
	case 21:
		indiv->fitness() = hf05(indiv->gene(), prob) + 2100.0;
		break;
	case 22:
		indiv->fitness() = hf06(indiv->gene(), prob) + 2200.0;
		break;
	case 23:
		indiv->fitness() = cf01(indiv->gene(), prob) + 2300.0;
		break;
	case 24:
		indiv->fitness() = cf02(indiv->gene(), prob) + 2400.0;
		break;
	case 25:
		indiv->fitness() = cf03(indiv->gene(), prob) + 2500.0;
		break;
	case 26:
		indiv->fitness() = cf04(indiv->gene(), prob) + 2600.0;
		break;
	case 27:
		indiv->fitness() = cf05(indiv->gene(), prob) + 2700.0;
		break;
	case 28:
		indiv->fitness() = cf06(indiv->gene(), prob) + 2800.0;
		break;
	case 29:
		indiv->fitness() = cf07(indiv->gene(), prob) + 2900.0;
		break;
	case 30:
		indiv->fitness() = cf08(indiv->gene(), prob) + 3000.0;
		break;
	default:
		printf("Can't find the corresponding function!\n");
		break;
	}                                             
}


const int CEC2019_100Digit_Score(Individual *indiv, const double &optimum)
{
	double fitness = indiv->fitness(); 
	for (int i = 0; i <= 9; ++i)
	{
		if (fitness - optimum >= pow(10.0, -i))
		{
			return i;
		}
	}
	return 10;
}


double Ellips(const Individual::GeneVec &gene_vec, const CProblem &prob) /* Ellipsoidal */
{
	Individual::GeneVec sr_gene_vec(gene_vec);
	sr_gene_vec = sr_func(sr_gene_vec, prob.shift_val_vec(), prob.rotate_matrix(), 1.0, prob.shiftable(), prob.rotatable());

	double sum = 0.0;
	size_t size_sub_one = sr_gene_vec.size() - 1;
	for (int i = 0; i < sr_gene_vec.size(); ++i)
	{
		sum += pow(10.0, 6.0 * i / size_sub_one) * sr_gene_vec[i] * sr_gene_vec[i];
	}
	return sum;
}

double Bent_Cigar(const Individual::GeneVec &gene_vec, const CProblem &prob) /* Bent_Cigar */
{
	Individual::GeneVec sr_gene_vec(gene_vec);
	sr_gene_vec = sr_func(sr_gene_vec, prob.shift_val_vec(), prob.rotate_matrix(), 1.0, prob.shiftable(), prob.rotatable());

	double sum = sr_gene_vec[0] * sr_gene_vec[0];
	double ten_pow_six = pow(10.0, 6.0);
	for (int i = 1; i < sr_gene_vec.size(); ++i)
	{
		sum += ten_pow_six * sr_gene_vec[i] * sr_gene_vec[i];
	}
	return sum;
}

double Discus(const Individual::GeneVec &gene_vec, const CProblem &prob) /* Discus */
{
	Individual::GeneVec sr_gene_vec(gene_vec);
	sr_gene_vec = sr_func(sr_gene_vec, prob.shift_val_vec(), prob.rotate_matrix(), 1.0, prob.shiftable(), prob.rotatable());

	double sum = pow(10.0, 6.0) * sr_gene_vec[0] * sr_gene_vec[0];
	for (int i = 1; i < sr_gene_vec.size(); ++i)
	{
		sum += sr_gene_vec[i] * sr_gene_vec[i];
	}
	return sum;
}

double Rosenbrock(const Individual::GeneVec &gene_vec, const CProblem &prob) /* Rosenbrock's */
{
	Individual::GeneVec sr_gene_vec(gene_vec);
	sr_gene_vec = sr_func(sr_gene_vec, prob.shift_val_vec(), prob.rotate_matrix(), 2.048 / 100.0, prob.shiftable(), prob.rotatable());

	double tmp1, tmp2, sum = 0.0;

	sr_gene_vec[0] += 1.0;
	for (int i = 0; i < sr_gene_vec.size() - 1; ++i)
	{
		sr_gene_vec[i + 1] += 1.0;

		tmp1 = sr_gene_vec[i] * sr_gene_vec[i] - sr_gene_vec[i + 1];
		tmp2 = sr_gene_vec[i] - 1.0;
		sum += 100.0 * tmp1 * tmp1 + tmp2 * tmp2;
	}
	return sum;
}

/*«ÝÀu¤Æ*/
double Katsuura(const Individual::GeneVec &gene_vec, const CProblem &prob) /* Katsuura  */
{
	Individual::GeneVec sr_gene_vec(gene_vec);
	sr_gene_vec = sr_func(sr_gene_vec, prob.shift_val_vec(), prob.rotate_matrix(), 5.0 / 100.0, prob.shiftable(), prob.rotatable());

	double temp, tmp1, tmp2, tmp3 = pow(1.0 * sr_gene_vec.size(), 1.2), sum = 1.0;
	std::vector<double> pow_j(33);
	for (int i = 0; i < pow_j.size(); ++i)
	{
		pow_j[i] = pow(2.0, i);
	}

	for (int i = 0; i < sr_gene_vec.size(); ++i)
	{
		temp = 0.0;
		for (int j = 1; j <= 32; ++j)
		{
			tmp2 = pow_j[j] * sr_gene_vec[i];
			temp += fabs(tmp2 - floor(tmp2 + 0.5)) / pow_j[j];
		}
		sum *= pow(1.0 + (i + 1)*temp, 10.0 / tmp3);
	}
	tmp1 = 10.0 / sr_gene_vec.size() / sr_gene_vec.size();
	sum = sum * tmp1 - tmp1;
	return sum;
}

double Hgbat(const Individual::GeneVec &gene_vec, const CProblem &prob) /* HGBat, provdided by Hans-Georg Beyer (HGB)*/
/* original global optimum: [-1,-1,...,-1] */
{
	Individual::GeneVec sr_gene_vec(gene_vec);
	sr_gene_vec = sr_func(sr_gene_vec, prob.shift_val_vec(), prob.rotate_matrix(), 5.0 / 100.0, prob.shiftable(), prob.rotatable());

	double alpha = 1.0 / 4.0, r2 = 0.0, sum_z = 0.0, sum = 0.0;
	double two_times_alpha = 2 * alpha;
	for (int i = 0; i < sr_gene_vec.size(); ++i)
	{
		sr_gene_vec[i] = sr_gene_vec[i] - 1.0;//shift to orgin
		r2 += sr_gene_vec[i] * sr_gene_vec[i];
		sum_z += sr_gene_vec[i];
	}
	sum = pow(fabs(pow(r2, 2.0) - pow(sum_z, 2.0)), two_times_alpha) + (0.5*r2 + sum_z) / sr_gene_vec.size() + 0.5;
	return sum;
}


double Grie_Rosen(const Individual::GeneVec &gene_vec, const CProblem &prob) /* Griewank-Rosenbrock  */
{
	Individual::GeneVec sr_gene_vec(gene_vec);
	sr_gene_vec = sr_func(sr_gene_vec, prob.shift_val_vec(), prob.rotate_matrix(), 5.0 / 100.0, prob.shiftable(), prob.rotatable());

	double temp, tmp1, tmp2, sum = 0.0;

	sr_gene_vec[0] += 1.0;
	for (int i = 0; i < sr_gene_vec.size() - 1; ++i)
	{
		sr_gene_vec[i + 1] += 1.0;
		tmp1 = sr_gene_vec[i] * sr_gene_vec[i] - sr_gene_vec[i + 1];
		tmp2 = sr_gene_vec[i] - 1.0;
		temp = 100.0*tmp1*tmp1 + tmp2 * tmp2;
		sum += (temp*temp) / 4000.0 - cos(temp) + 1.0;
	}
	tmp1 = sr_gene_vec[sr_gene_vec.size() - 1] * sr_gene_vec[sr_gene_vec.size() - 1] - sr_gene_vec[0];
	tmp2 = sr_gene_vec[sr_gene_vec.size() - 1] - 1.0;
	temp = 100.0 * tmp1 * tmp1 + tmp2 * tmp2;
	sum += (temp*temp) / 4000.0 - cos(temp) + 1.0;
	return sum;
}

double hf01(const Individual::GeneVec &gene_vec, const CProblem &prob) /* hybrid function 01 */
{
	int tmp, cf_num = 3, G[3], G_nx[3];
	double fit[3], Gp[3] = { 0.3, 0.3, 0.4 };

	tmp = 0;
	for (int i = 0; i < cf_num - 1; ++i)
	{
		G_nx[i] = ceil(Gp[i] * prob.dim());
		tmp += G_nx[i];
	}
	G_nx[cf_num - 1] = prob.dim() - tmp;

	G[0] = 0;
	for (int i = 1; i < cf_num; ++i)
	{
		G[i] = G[i - 1] + G_nx[i - 1];
	}
	
	Individual::GeneVec sr_gene_vec(gene_vec);
	sr_gene_vec = sr_func(sr_gene_vec, prob.shift_val_vec(), prob.rotate_matrix(), 1.0, prob.shiftable(), prob.rotatable());
	if(prob.shuffleable()) 	sr_gene_vec = shuffle_func(sr_gene_vec, prob.shuffle_val_vec());

	CProblem hybrid_prob = prob;
	hybrid_prob.set_shift_flag(false);
	hybrid_prob.set_rotate_flag(false);
	
	fit[0] = Modified_Schwefel(Individual::GeneVec(sr_gene_vec.begin() + G[0], sr_gene_vec.begin() + G[0] + G_nx[0]), hybrid_prob);
	fit[1] = Rastrihin(Individual::GeneVec(sr_gene_vec.begin() + G[1], sr_gene_vec.begin() + G[1] + G_nx[1]), hybrid_prob);
	fit[2] = Ellips(Individual::GeneVec(sr_gene_vec.begin() + G[2], sr_gene_vec.begin() + G[2] + G_nx[2]), hybrid_prob);

	double sum = 0.0;
	for (int i = 0; i < cf_num; ++i) sum += fit[i];
	return sum;
}

double hf02(const Individual::GeneVec &gene_vec, const CProblem &prob) /* hybrid function 02 */
{
	int tmp, cf_num = 3, G[3], G_nx[3];
	double fit[3], Gp[3] = { 0.3, 0.3, 0.4 };

	tmp = 0;
	for (int i = 0; i < cf_num - 1; ++i)
	{
		G_nx[i] = ceil(Gp[i] * prob.dim());
		tmp += G_nx[i];
	}
	G_nx[cf_num - 1] = prob.dim() - tmp;

	G[0] = 0;
	for (int i = 1; i < cf_num; ++i)
	{
		G[i] = G[i - 1] + G_nx[i - 1];
	}

	Individual::GeneVec sr_gene_vec(gene_vec);
	sr_gene_vec = sr_func(sr_gene_vec, prob.shift_val_vec(), prob.rotate_matrix(), 1.0, prob.shiftable(), prob.rotatable());
	if (prob.shuffleable()) 	sr_gene_vec = shuffle_func(sr_gene_vec, prob.shuffle_val_vec());

	CProblem hybrid_prob = prob;
	hybrid_prob.set_shift_flag(false);
	hybrid_prob.set_rotate_flag(false);

	fit[0] = Bent_Cigar(Individual::GeneVec(sr_gene_vec.begin() + G[0], sr_gene_vec.begin() + G[0] + G_nx[0]), hybrid_prob);
	fit[1] = Hgbat(Individual::GeneVec(sr_gene_vec.begin() + G[1], sr_gene_vec.begin() + G[1] + G_nx[1]), hybrid_prob);
	fit[2] = Rastrihin(Individual::GeneVec(sr_gene_vec.begin() + G[2], sr_gene_vec.begin() + G[2] + G_nx[2]), hybrid_prob);

	double sum = 0.0;
	for (int i = 0; i < cf_num; ++i) sum += fit[i];
	return sum;
}

double hf03(const Individual::GeneVec &gene_vec, const CProblem &prob) /* hybrid function 03 */
{
	int tmp, cf_num = 4, G[4], G_nx[4];
	double fit[4], Gp[4] = { 0.2, 0.2, 0.3, 0.3 };

	tmp = 0;
	for (int i = 0; i < cf_num - 1; ++i)
	{
		G_nx[i] = ceil(Gp[i] * prob.dim());
		tmp += G_nx[i];
	}
	G_nx[cf_num - 1] = prob.dim() - tmp;

	G[0] = 0;
	for (int i = 1; i < cf_num; ++i)
	{
		G[i] = G[i - 1] + G_nx[i - 1];
	}

	Individual::GeneVec sr_gene_vec(gene_vec);
	sr_gene_vec = sr_func(sr_gene_vec, prob.shift_val_vec(), prob.rotate_matrix(), 1.0, prob.shiftable(), prob.rotatable());
	if (prob.shuffleable()) 	sr_gene_vec = shuffle_func(sr_gene_vec, prob.shuffle_val_vec());

	CProblem hybrid_prob = prob;
	hybrid_prob.set_shift_flag(false);
	hybrid_prob.set_rotate_flag(false);

	fit[0] = Griewangk(Individual::GeneVec(sr_gene_vec.begin() + G[0], sr_gene_vec.begin() + G[0] + G_nx[0]), hybrid_prob);
	fit[1] = Weierstrass(Individual::GeneVec(sr_gene_vec.begin() + G[1], sr_gene_vec.begin() + G[1] + G_nx[1]), hybrid_prob);
	fit[2] = Rosenbrock(Individual::GeneVec(sr_gene_vec.begin() + G[2], sr_gene_vec.begin() + G[2] + G_nx[2]), hybrid_prob);
	fit[3] = Expanded_Schaffer_F6(Individual::GeneVec(sr_gene_vec.begin() + G[3], sr_gene_vec.begin() + G[3] + G_nx[3]), hybrid_prob);

	double sum = 0.0;
	for (int i = 0; i < cf_num; ++i) sum += fit[i];
	return sum;
}

double hf04(const Individual::GeneVec &gene_vec, const CProblem &prob) /* hybrid function 04 */
{
	int tmp, cf_num = 4, G[4], G_nx[4];
	double fit[4], Gp[4] = { 0.2, 0.2, 0.3, 0.3 };

	tmp = 0;
	for (int i = 0; i < cf_num - 1; ++i)
	{
		G_nx[i] = ceil(Gp[i] * prob.dim());
		tmp += G_nx[i];
	}
	G_nx[cf_num - 1] = prob.dim() - tmp;

	G[0] = 0;
	for (int i = 1; i < cf_num; ++i)
	{
		G[i] = G[i - 1] + G_nx[i - 1];
	}

	Individual::GeneVec sr_gene_vec(gene_vec);
	sr_gene_vec = sr_func(sr_gene_vec, prob.shift_val_vec(), prob.rotate_matrix(), 1.0, prob.shiftable(), prob.rotatable());
	if (prob.shuffleable()) 	sr_gene_vec = shuffle_func(sr_gene_vec, prob.shuffle_val_vec());

	CProblem hybrid_prob = prob;
	hybrid_prob.set_shift_flag(false);
	hybrid_prob.set_rotate_flag(false);

	fit[0] = Hgbat(Individual::GeneVec(sr_gene_vec.begin() + G[0], sr_gene_vec.begin() + G[0] + G_nx[0]), hybrid_prob);
	fit[1] = Discus(Individual::GeneVec(sr_gene_vec.begin() + G[1], sr_gene_vec.begin() + G[1] + G_nx[1]), hybrid_prob);
	fit[2] = Grie_Rosen(Individual::GeneVec(sr_gene_vec.begin() + G[2], sr_gene_vec.begin() + G[2] + G_nx[2]), hybrid_prob);
	fit[3] = Rastrihin(Individual::GeneVec(sr_gene_vec.begin() + G[3], sr_gene_vec.begin() + G[3] + G_nx[3]), hybrid_prob);

	double sum = 0.0;
	for (int i = 0; i < cf_num; ++i) sum += fit[i];
	return sum;
}

double hf05(const Individual::GeneVec &gene_vec, const CProblem &prob) /* hybrid function 05 */
{
	int tmp, cf_num = 5, G[5], G_nx[5];
	double fit[5], Gp[5] = { 0.1, 0.2, 0.2, 0.2, 0.3 };

	tmp = 0;
	for (int i = 0; i < cf_num - 1; ++i)
	{
		G_nx[i] = ceil(Gp[i] * prob.dim());
		tmp += G_nx[i];
	}
	G_nx[cf_num - 1] = prob.dim() - tmp;

	G[0] = 0;
	for (int i = 1; i < cf_num; ++i)
	{
		G[i] = G[i - 1] + G_nx[i - 1];
	}

	Individual::GeneVec sr_gene_vec(gene_vec);
	sr_gene_vec = sr_func(sr_gene_vec, prob.shift_val_vec(), prob.rotate_matrix(), 1.0, prob.shiftable(), prob.rotatable());
	if (prob.shuffleable()) 	sr_gene_vec = shuffle_func(sr_gene_vec, prob.shuffle_val_vec());

	CProblem hybrid_prob = prob;
	hybrid_prob.set_shift_flag(false);
	hybrid_prob.set_rotate_flag(false);

	fit[0] = Expanded_Schaffer_F6(Individual::GeneVec(sr_gene_vec.begin() + G[0], sr_gene_vec.begin() + G[0] + G_nx[0]), hybrid_prob);
	fit[1] = Hgbat(Individual::GeneVec(sr_gene_vec.begin() + G[1], sr_gene_vec.begin() + G[1] + G_nx[1]), hybrid_prob);
	fit[2] = Rosenbrock(Individual::GeneVec(sr_gene_vec.begin() + G[2], sr_gene_vec.begin() + G[2] + G_nx[2]), hybrid_prob);
	fit[3] = Modified_Schwefel(Individual::GeneVec(sr_gene_vec.begin() + G[3], sr_gene_vec.begin() + G[3] + G_nx[3]), hybrid_prob);
	fit[4] = Ellips(Individual::GeneVec(sr_gene_vec.begin() + G[4], sr_gene_vec.begin() + G[4] + G_nx[4]), hybrid_prob);

	double sum = 0.0;
	for (int i = 0; i < cf_num; ++i) sum += fit[i];
	return sum;
}

double hf06(const Individual::GeneVec &gene_vec, const CProblem &prob) /* hybrid function 06 */
{
	int tmp, cf_num = 5, G[5], G_nx[5];
	double fit[5], Gp[5] = { 0.1, 0.2, 0.2, 0.2, 0.3 };

	tmp = 0;
	for (int i = 0; i < cf_num - 1; ++i)
	{
		G_nx[i] = ceil(Gp[i] * prob.dim());
		tmp += G_nx[i];
	}
	G_nx[cf_num - 1] = prob.dim() - tmp;

	G[0] = 0;
	for (int i = 1; i < cf_num; ++i)
	{
		G[i] = G[i - 1] + G_nx[i - 1];
	}

	Individual::GeneVec sr_gene_vec(gene_vec);
	sr_gene_vec = sr_func(sr_gene_vec, prob.shift_val_vec(), prob.rotate_matrix(), 1.0, prob.shiftable(), prob.rotatable());
	if (prob.shuffleable()) 	sr_gene_vec = shuffle_func(sr_gene_vec, prob.shuffle_val_vec());

	CProblem hybrid_prob = prob;
	hybrid_prob.set_shift_flag(false);
	hybrid_prob.set_rotate_flag(false);

	fit[0] = Katsuura(Individual::GeneVec(sr_gene_vec.begin() + G[0], sr_gene_vec.begin() + G[0] + G_nx[0]), hybrid_prob);
	fit[1] = Happy_Cat(Individual::GeneVec(sr_gene_vec.begin() + G[1], sr_gene_vec.begin() + G[1] + G_nx[1]), hybrid_prob);
	fit[2] = Grie_Rosen(Individual::GeneVec(sr_gene_vec.begin() + G[2], sr_gene_vec.begin() + G[2] + G_nx[2]), hybrid_prob);
	fit[3] = Modified_Schwefel(Individual::GeneVec(sr_gene_vec.begin() + G[3], sr_gene_vec.begin() + G[3] + G_nx[3]), hybrid_prob);
	fit[4] = Ackley(Individual::GeneVec(sr_gene_vec.begin() + G[4], sr_gene_vec.begin() + G[4] + G_nx[4]), hybrid_prob);

	double sum = 0.0;
	for (int i = 0; i < cf_num; ++i) sum += fit[i];
	return sum;
}


double cf01(const Individual::GeneVec &gene_vec, const CProblem &prob) /* composition function 01 */
{
	int cf_num = 5;
	std::vector<double> fit(5);
	std::vector<double> delta = { 10, 20, 30, 40, 50 };
	std::vector<double> bias = { 0, 100, 200, 300, 400 };

	CProblem comp_prob = prob;
	comp_prob.set_rotate_flag(prob.rotatable());
	comp_prob.set_shift_flag(true);

	std::vector<double>::const_iterator sfit = prob.shift_val_vec().begin(), seit = prob.shift_val_vec().begin() + prob.dim();
	std::vector<CProblem::RNV>::const_iterator rfit = prob.rotate_matrix().begin(), reit = prob.rotate_matrix().begin() + prob.dim();

	comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
	comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
	fit[0] = 10000 * Rosenbrock(gene_vec, comp_prob) / 1e+4;
	
	sfit = seit; seit += prob.dim();
	rfit = reit; reit += prob.dim();

	comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
	comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
	fit[1] = 10000 * Ellips(gene_vec, comp_prob) / 1e+10;
	
	sfit = seit; seit += prob.dim();
	rfit = reit; reit += prob.dim();

	comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
	comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
	fit[2] = 10000 * Bent_Cigar(gene_vec, comp_prob) / 1e+30;
	
	sfit = seit; seit += prob.dim();
	rfit = reit; reit += prob.dim();

	comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
	comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
	fit[3] = 10000 * Discus(gene_vec, comp_prob) / 1e+10;
	
	comp_prob.set_rotate_flag(false);

	sfit = seit; seit += prob.dim();
	rfit = reit; reit += prob.dim();

	comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
	comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
	fit[4] = 10000 * Ellips(gene_vec, comp_prob) / 1e+10;
	return cf_cal(gene_vec, prob.shift_val_vec(), fit, delta, bias);
}

double cf02(const Individual::GeneVec &gene_vec, const CProblem &prob) /* composition function 02 */
{
	int cf_num = 3;
	std::vector<double> fit(3);
	std::vector<double> delta = { 20, 20, 20 };
	std::vector<double> bias = { 0, 100, 200 };

	CProblem comp_prob = prob;
	comp_prob.set_rotate_flag(false);
	comp_prob.set_shift_flag(true);

	std::vector<double>::const_iterator sfit = prob.shift_val_vec().begin(), seit = prob.shift_val_vec().begin() + prob.dim();
	std::vector<CProblem::RNV>::const_iterator rfit = prob.rotate_matrix().begin(), reit = prob.rotate_matrix().begin() + prob.dim();

	comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
	comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
	fit[0] = Modified_Schwefel(gene_vec, comp_prob);

	comp_prob.set_rotate_flag(prob.rotatable());

	sfit = seit; seit += prob.dim();
	rfit = reit; reit += prob.dim();

	comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
	comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
	fit[1] = Rastrihin(gene_vec, comp_prob);

	sfit = seit; seit += prob.dim();
	rfit = reit; reit += prob.dim();

	comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
	comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
	fit[2] = Hgbat(gene_vec, comp_prob);

	return cf_cal(gene_vec, prob.shift_val_vec(), fit, delta, bias);
}

double cf03(const Individual::GeneVec &gene_vec, const CProblem &prob) /* composition function 03 */
{
	int cf_num = 3;
	std::vector<double> fit(3);
	std::vector<double> delta = { 10, 30, 50 };
	std::vector<double> bias = { 0, 100, 200 };

	CProblem comp_prob = prob;
	comp_prob.set_rotate_flag(prob.rotatable());
	comp_prob.set_shift_flag(true);

	std::vector<double>::const_iterator sfit = prob.shift_val_vec().begin(), seit = prob.shift_val_vec().begin() + prob.dim();
	std::vector<CProblem::RNV>::const_iterator rfit = prob.rotate_matrix().begin(), reit = prob.rotate_matrix().begin() + prob.dim();

	comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
	comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
	fit[0] = 1000 * Modified_Schwefel(gene_vec, comp_prob) / 4e+3;

	sfit = seit; seit += prob.dim();
	rfit = reit; reit += prob.dim();

	comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
	comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
	fit[1] = 1000 * Rastrihin(gene_vec, comp_prob) / 1e+3;

	sfit = seit; seit += prob.dim();
	rfit = reit; reit += prob.dim();

	comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
	comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
	fit[2] = 1000 * Ellips(gene_vec, comp_prob) / 1e+10;

	return cf_cal(gene_vec, prob.shift_val_vec(), fit, delta, bias);
}

double cf04(const Individual::GeneVec &gene_vec, const CProblem &prob) /* composition function 04 */
{
	int cf_num = 5;
	std::vector<double> fit(5);
	std::vector<double> delta = { 10, 10, 10, 10, 10 };
	std::vector<double> bias = { 0, 100, 200, 300, 400 };

	CProblem comp_prob = prob;
	comp_prob.set_rotate_flag(prob.rotatable());
	comp_prob.set_shift_flag(true);

	std::vector<double>::const_iterator sfit = prob.shift_val_vec().begin(), seit = prob.shift_val_vec().begin() + prob.dim();
	std::vector<CProblem::RNV>::const_iterator rfit = prob.rotate_matrix().begin(), reit = prob.rotate_matrix().begin() + prob.dim();

	comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
	comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
	fit[0] = 1000 * Modified_Schwefel(gene_vec, comp_prob) / 4e+3;
	sfit = seit; seit += prob.dim();
	rfit = reit; reit += prob.dim();

	comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
	comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
	fit[1] = 1000 * Happy_Cat(gene_vec, comp_prob) / 1e+3;

	sfit = seit; seit += prob.dim();
	rfit = reit; reit += prob.dim();

	comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
	comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
	fit[2] = 1000 * Ellips(gene_vec, comp_prob) / 1e+10;

	sfit = seit; seit += prob.dim();
	rfit = reit; reit += prob.dim();

	comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
	comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
	fit[3] = 1000 * Weierstrass(gene_vec, comp_prob) / 400;

	sfit = seit; seit += prob.dim();
	rfit = reit; reit += prob.dim();

	comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
	comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
	fit[4] = 1000 * Griewangk(gene_vec, comp_prob) / 100;

	return cf_cal(gene_vec, prob.shift_val_vec(), fit, delta, bias);
}

double cf05(const Individual::GeneVec &gene_vec, const CProblem &prob) /* composition function 05 */
{
	int cf_num = 5;
	std::vector<double> fit(5);
	std::vector<double> delta = { 10, 10, 10, 20, 20 };
	std::vector<double> bias = { 0, 100, 200, 300, 400 };

	CProblem comp_prob = prob;
	comp_prob.set_rotate_flag(prob.rotatable());
	comp_prob.set_shift_flag(true);

	std::vector<double>::const_iterator sfit = prob.shift_val_vec().begin(), seit = prob.shift_val_vec().begin() + prob.dim();
	std::vector<CProblem::RNV>::const_iterator rfit = prob.rotate_matrix().begin(), reit = prob.rotate_matrix().begin() + prob.dim();

	comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
	comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
	fit[0] = 10000 * Hgbat(gene_vec, comp_prob) / 1000;

	sfit = seit; seit += prob.dim();
	rfit = reit; reit += prob.dim();

	comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
	comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
	fit[1] = 10000 * Rastrihin(gene_vec, comp_prob) / 1e+3;

	sfit = seit; seit += prob.dim();
	rfit = reit; reit += prob.dim();

	comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
	comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
	fit[2] = 10000 * Modified_Schwefel(gene_vec, comp_prob) / 4e+3;

	sfit = seit; seit += prob.dim();
	rfit = reit; reit += prob.dim();

	comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
	comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
	fit[3] = 10000 * Weierstrass(gene_vec, comp_prob) / 400;

	sfit = seit; seit += prob.dim();
	rfit = reit; reit += prob.dim();

	comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
	comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
	fit[4] = 10000 * Ellips(gene_vec, comp_prob) / 1e+10;

	return cf_cal(gene_vec, prob.shift_val_vec(), fit, delta, bias);
}

double cf06(const Individual::GeneVec &gene_vec, const CProblem &prob) /* composition function 06 */
{
	int cf_num = 5;
	std::vector<double> fit(5);
	std::vector<double> delta = { 10, 20, 30, 40, 50 };
	std::vector<double> bias = { 0, 100, 200, 300, 400 };

	CProblem comp_prob = prob;
	comp_prob.set_rotate_flag(prob.rotatable());
	comp_prob.set_shift_flag(true);

	std::vector<double>::const_iterator sfit = prob.shift_val_vec().begin(), seit = prob.shift_val_vec().begin() + prob.dim();
	std::vector<CProblem::RNV>::const_iterator rfit = prob.rotate_matrix().begin(), reit = prob.rotate_matrix().begin() + prob.dim();

	comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
	comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
	fit[0] = 10000 * Grie_Rosen(gene_vec, comp_prob) / 4e+3;

	sfit = seit; seit += prob.dim();
	rfit = reit; reit += prob.dim();

	comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
	comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
	fit[1] = 10000 * Happy_Cat(gene_vec, comp_prob) / 1e+3;

	sfit = seit; seit += prob.dim();
	rfit = reit; reit += prob.dim();

	comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
	comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
	fit[2] = 10000 * Modified_Schwefel(gene_vec, comp_prob) / 4e+3;

	sfit = seit; seit += prob.dim();
	rfit = reit; reit += prob.dim();

	comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
	comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
	fit[3] = 10000 * Expanded_Schaffer_F6(gene_vec, comp_prob) / 2e+7;

	sfit = seit; seit += prob.dim();
	rfit = reit; reit += prob.dim();

	comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
	comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
	fit[4] = 10000 * Ellips(gene_vec, comp_prob) / 1e+10;

	return cf_cal(gene_vec, prob.shift_val_vec(), fit, delta, bias);
}

double cf07(const Individual::GeneVec &gene_vec, const CProblem &prob) /* composition function 07 */
{
	int cf_num = 3;
	std::vector<double> fit(3);
	std::vector<double> delta = { 10, 30, 50 };
	std::vector<double> bias = { 0, 100, 200 };

	CProblem comp_prob = prob;
	comp_prob.set_rotate_flag(prob.rotatable());
	comp_prob.set_shift_flag(true);

	std::vector<double>::const_iterator sfit = prob.shift_val_vec().begin(), seit = prob.shift_val_vec().begin() + prob.dim();
	std::vector<CProblem::RNV>::const_iterator rfit = prob.rotate_matrix().begin(), reit = prob.rotate_matrix().begin() + prob.dim();
	std::vector<int>::const_iterator shfit = prob.shuffle_val_vec().begin(), sheit = prob.shuffle_val_vec().begin() + prob.dim();

	comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
	comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
	comp_prob.set_shuffle_vector(std::vector<int>(shfit, sheit));
	fit[0] = hf01(gene_vec, comp_prob);

	sfit = seit; seit += prob.dim();
	rfit = reit; reit += prob.dim();
	shfit = sheit; sheit += prob.dim();

	comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
	comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
	comp_prob.set_shuffle_vector(std::vector<int>(shfit, sheit));
	fit[1] = hf02(gene_vec, comp_prob);

	sfit = seit; seit += prob.dim();
	rfit = reit; reit += prob.dim();
	shfit = sheit; sheit += prob.dim();

	comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
	comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
	comp_prob.set_shuffle_vector(std::vector<int>(shfit, sheit));
	fit[2] = hf03(gene_vec, comp_prob);

	return cf_cal(gene_vec, prob.shift_val_vec(), fit, delta, bias);
}

double cf08(const Individual::GeneVec &gene_vec, const CProblem &prob) /* composition function 08 */
{
	int cf_num = 3;
	std::vector<double> fit(3);
	std::vector<double> delta = { 10, 30, 50 };
	std::vector<double> bias = { 0, 100, 200 };

	CProblem comp_prob = prob;
	comp_prob.set_rotate_flag(prob.rotatable());
	comp_prob.set_shift_flag(true);

	std::vector<double>::const_iterator sfit = prob.shift_val_vec().begin(), seit = prob.shift_val_vec().begin() + prob.dim();
	std::vector<CProblem::RNV>::const_iterator rfit = prob.rotate_matrix().begin(), reit = prob.rotate_matrix().begin() + prob.dim();
	std::vector<int>::const_iterator shfit = prob.shuffle_val_vec().begin(), sheit = prob.shuffle_val_vec().begin() + prob.dim();

	comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
	comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
	comp_prob.set_shuffle_vector(std::vector<int>(shfit, sheit));
	fit[0] = hf04(gene_vec, comp_prob);

	sfit = seit; seit += prob.dim();
	rfit = reit; reit += prob.dim();
	shfit = sheit; sheit += prob.dim();

	comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
	comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
	comp_prob.set_shuffle_vector(std::vector<int>(shfit, sheit));
	fit[1] = hf05(gene_vec, comp_prob);

	sfit = seit; seit += prob.dim();
	rfit = reit; reit += prob.dim();
	shfit = sheit; sheit += prob.dim();

	comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
	comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
	comp_prob.set_shuffle_vector(std::vector<int>(shfit, sheit));
	fit[2] = hf06(gene_vec, comp_prob);

	return cf_cal(gene_vec, prob.shift_val_vec(), fit, delta, bias);
}

double cf_cal(const Individual::GeneVec &gene_vec, const CProblem::RNV &shift_vec, const std::vector<double> &fit_vec, const std::vector<double> &delta, const std::vector<double> &bias)
{
	std::vector<double> w(fit_vec.size());
	double w_max = 0.0, w_sum = 0.0;
	std::vector<double> new_fit_vec = fit_vec;

	for (int i = 0; i < fit_vec.size(); ++i)
	{
		new_fit_vec[i] += bias[i];
		w[i] = 0;
		for (int j = 0; j < gene_vec.size(); ++j)
		{
			w[i] += pow(gene_vec[j] - shift_vec[i*gene_vec.size() + j], 2.0);
		}
		if (w[i] != 0) w[i] = pow(1.0 / w[i], 0.5) * exp(-w[i] / 2.0 / gene_vec.size() / pow(delta[i], 2.0));
		else           w[i] = 1.0e99;
		if (w[i] > w_max) w_max = w[i];
	}

	for (int i = 0; i < fit_vec.size(); ++i)
	{
		w_sum += w[i];
	}

	if (w_max == 0)
	{
		for (int i = 0; i < fit_vec.size(); ++i) w[i] = 1;
		w_sum = fit_vec.size();
	}
	double sum = 0.0;
	for (int i = 0; i < fit_vec.size(); ++i)
	{
		sum += w[i] / w_sum * new_fit_vec[i];
	}
	return sum;
}




/******************************************************************************************/
//	2019/12/24 ·s¼Wfunction

namespace CEC2020 {

	void CEC2020_BoundConstrained_Evaluate(Individual* indiv, const CProblem& prob)
	{
		if (prob.id() == 1)				indiv->fitness() = Bent_Cigar(indiv->gene(), prob) + prob.global_optimum();
		else if (prob.id() == 2)		indiv->fitness() = Modified_Schwefel(indiv->gene(), prob) + prob.global_optimum();
		else if (prob.id() == 3)		indiv->fitness() = CEC2020::Lunacek_bi_Rastrigin(indiv->gene(), prob) + prob.global_optimum();
		else if (prob.id() == 4)		indiv->fitness() = Grie_Rosen(indiv->gene(), prob) + prob.global_optimum();
		else if (prob.id() == 5)		indiv->fitness() = CEC2020::hf01(indiv->gene(), prob) + prob.global_optimum();
		else if (prob.id() == 6)		indiv->fitness() = CEC2020::hf06(indiv->gene(), prob) + prob.global_optimum();
		else if (prob.id() == 7)		indiv->fitness() = CEC2020::hf05(indiv->gene(), prob) + prob.global_optimum();
		else if (prob.id() == 8)		indiv->fitness() = CEC2020::cf02(indiv->gene(), prob) + prob.global_optimum();
		else if (prob.id() == 9)		indiv->fitness() = CEC2020::cf04(indiv->gene(), prob) + prob.global_optimum();
		else if (prob.id() == 10)		indiv->fitness() = CEC2020::cf05(indiv->gene(), prob) + prob.global_optimum();
		else							printf("Can't find the corresponding function!\n");
	}


	double Lunacek_bi_Rastrigin(const Individual::GeneVec& gene_vec, const CProblem& prob)
	{
		Individual::GeneVec shifted_gene_vec(gene_vec);

		double mu0 = 2.5, d = 1.0, s, mu1;
		s = 1.0 - 1.0 / (2.0 * pow(shifted_gene_vec.size() + 20.0, 0.5) - 8.2);
		mu1 = -pow((mu0 * mu0 - d) / s, 0.5);

		if (prob.shiftable())
			shifted_gene_vec = shiftfunc(shifted_gene_vec, prob.shift_val_vec());

		for (int i = 0; i < shifted_gene_vec.size(); i++)
		{
			shifted_gene_vec[i] *= 10.0 / 100.0;
		}


		Individual::GeneVec tmpx(gene_vec);

		for (int i = 0; i < tmpx.size(); i++)
		{
			tmpx[i] = 2 * shifted_gene_vec[i];
			if (prob.shift_val_vec()[i] < 0.0)
				tmpx[i] *= -1.;
		}

		Individual::GeneVec rotated_gene_vec(tmpx);

		for (int i = 0; i < tmpx.size(); i++)
		{
			tmpx[i] += mu0;
		}

		double tmp = 0.0, tmp1 = 0.0, tmp2 = 0.0;

		for (int i = 0; i < tmpx.size(); i++)
		{
			tmp = tmpx[i] - mu0;
			tmp1 += tmp * tmp;
			tmp = tmpx[i] - mu1;
			tmp2 += tmp * tmp;
		}
		tmp2 *= s;
		tmp2 += d * (int)tmpx.size();
		tmp = 0.0;

		if (prob.rotatable())
		{
			rotated_gene_vec = rotatefunc(rotated_gene_vec, prob.rotate_matrix());
		}

		double sum = 0;
		for (int i = 0; i < rotated_gene_vec.size(); i++)
		{
			tmp += cos(2.0 * PI * rotated_gene_vec[i]);
		}
		if (tmp1 < tmp2)
			sum = tmp1;
		else
			sum = tmp2;
		sum += 10.0 * ((int)rotated_gene_vec.size() - tmp);

		return sum;
	}

	double hf01(const Individual::GeneVec& gene_vec, const CProblem& prob) /* hybrid function 1 */
	{
		int tmp, cf_num = 3, G[3], G_nx[3];
		double fit[3], Gp[3] = { 0.3, 0.3, 0.4 };

		tmp = 0;
		for (int i = 1; i < cf_num; ++i)
		{
			G_nx[i] = ceil(Gp[i] * prob.dim());
			tmp += G_nx[i];
		}
		//G_nx[cf_num - 1] = prob.dim() - tmp;
		G_nx[0] = prob.dim() - tmp;

		G[0] = 0;
		for (int i = 1; i < cf_num; ++i)
		{
			G[i] = G[i - 1] + G_nx[i - 1];
		}

		Individual::GeneVec sr_gene_vec(gene_vec);
		sr_gene_vec = sr_func(sr_gene_vec, prob.shift_val_vec(), prob.rotate_matrix(), 1.0, prob.shiftable(), prob.rotatable());
		if (prob.shuffleable()) 	sr_gene_vec = shuffle_func(sr_gene_vec, prob.shuffle_val_vec());

		CProblem hybrid_prob = prob;
		hybrid_prob.set_shift_flag(false);
		hybrid_prob.set_rotate_flag(false);

		fit[0] = Modified_Schwefel(Individual::GeneVec(sr_gene_vec.begin() + G[0], sr_gene_vec.begin() + G[0] + G_nx[0]), hybrid_prob);
		fit[1] = Rastrihin(Individual::GeneVec(sr_gene_vec.begin() + G[1], sr_gene_vec.begin() + G[1] + G_nx[1]), hybrid_prob);
		fit[2] = Ellips(Individual::GeneVec(sr_gene_vec.begin() + G[2], sr_gene_vec.begin() + G[2] + G_nx[2]), hybrid_prob);

		double sum = 0.0;
		for (int i = 0; i < cf_num; ++i) sum += fit[i];
		return sum;
	}

	double hf06(const Individual::GeneVec& gene_vec, const CProblem& prob) /* hybrid function 2 */
	{
		int tmp, cf_num = 4, G[4], G_nx[4];
		double fit[4], Gp[4] = { 0.2, 0.2, 0.3, 0.3 };

		if (prob.dim() == 5)
		{
			G_nx[0] = 1;
			G_nx[1] = 1;
			G_nx[2] = 1;
			G_nx[3] = 2;
		}
		else
		{
			tmp = 0;
			for (int i = 1; i < cf_num ; ++i)
			{
				G_nx[i] = ceil(Gp[i] * prob.dim());
				tmp += G_nx[i];
			}
			G_nx[0] = prob.dim() - tmp;	
		}	

		G[0] = 0;
		for (int i = 1; i < cf_num; ++i)
		{
			G[i] = G[i - 1] + G_nx[i - 1];
		}

		Individual::GeneVec sr_gene_vec(gene_vec);
		sr_gene_vec = sr_func(sr_gene_vec, prob.shift_val_vec(), prob.rotate_matrix(), 1.0, prob.shiftable(), prob.rotatable());
		if (prob.shuffleable()) 	sr_gene_vec = shuffle_func(sr_gene_vec, prob.shuffle_val_vec());

		CProblem hybrid_prob = prob;
		hybrid_prob.set_shift_flag(false);
		hybrid_prob.set_rotate_flag(false);

		fit[0] = Expanded_Schaffer_F6(Individual::GeneVec(sr_gene_vec.begin() + G[0], sr_gene_vec.begin() + G[0] + G_nx[0]), hybrid_prob);
		fit[1] = Hgbat(Individual::GeneVec(sr_gene_vec.begin() + G[1], sr_gene_vec.begin() + G[1] + G_nx[1]), hybrid_prob);
		fit[2] = Rosenbrock(Individual::GeneVec(sr_gene_vec.begin() + G[2], sr_gene_vec.begin() + G[2] + G_nx[2]), hybrid_prob);
		fit[3] = Modified_Schwefel(Individual::GeneVec(sr_gene_vec.begin() + G[3], sr_gene_vec.begin() + G[3] + G_nx[3]), hybrid_prob);

		double sum = 0.0;
		for (int i = 0; i < cf_num; ++i) sum += fit[i];
		return sum;
	}

	double hf05(const Individual::GeneVec& gene_vec, const CProblem& prob) /* hybrid function 3 */
	{
		int tmp, cf_num = 5, G[5], G_nx[5];
		double fit[5], Gp[5] = { 0.1, 0.2, 0.2, 0.2, 0.3 };

		if (prob.dim() == 5)
		{
			for (int i = 0; i < 5; ++i)
			{
				G_nx[i] = 1;
			}
		}
		else
		{
			tmp = 0;
			for (int i = 1; i < cf_num; ++i)
			{
				G_nx[i] = ceil(Gp[i] * prob.dim());
				tmp += G_nx[i];
			}
			G_nx[0] = prob.dim() - tmp;
		}

		G[0] = 0;
		for (int i = 1; i < cf_num; ++i)
		{
			G[i] = G[i - 1] + G_nx[i - 1];
		}
		
		Individual::GeneVec sr_gene_vec(gene_vec);
		sr_gene_vec = sr_func(sr_gene_vec, prob.shift_val_vec(), prob.rotate_matrix(), 1.0, prob.shiftable(), prob.rotatable());
		if (prob.shuffleable()) 	sr_gene_vec = shuffle_func(sr_gene_vec, prob.shuffle_val_vec());

		CProblem hybrid_prob = prob;
		hybrid_prob.set_shift_flag(false);
		hybrid_prob.set_rotate_flag(false);
		
		//in D05 will cause some trouble in fit[0] and the program will crash here...
		fit[0] = Expanded_Schaffer_F6(Individual::GeneVec(sr_gene_vec.begin() + G[0], sr_gene_vec.begin() + G[0] + G_nx[0]), hybrid_prob);
		fit[1] = Hgbat(Individual::GeneVec(sr_gene_vec.begin() + G[1], sr_gene_vec.begin() + G[1] + G_nx[1]), hybrid_prob);
		fit[2] = Rosenbrock(Individual::GeneVec(sr_gene_vec.begin() + G[2], sr_gene_vec.begin() + G[2] + G_nx[2]), hybrid_prob);
		fit[3] = Modified_Schwefel(Individual::GeneVec(sr_gene_vec.begin() + G[3], sr_gene_vec.begin() + G[3] + G_nx[3]), hybrid_prob);
		fit[4] = Ellips(Individual::GeneVec(sr_gene_vec.begin() + G[4], sr_gene_vec.begin() + G[4] + G_nx[4]), hybrid_prob);
		
		double sum = 0.0;
		for (int i = 0; i < cf_num; ++i) sum += fit[i];
		return sum;
	}

	double cf02(const Individual::GeneVec& gene_vec, const CProblem& prob) /* composition function 1 */
	{
		int cf_num = 3;
		std::vector<double> fit(3);
		std::vector<double> delta = { 10, 20, 30 };
		std::vector<double> bias = { 0, 100, 200 };

		CProblem comp_prob = prob;
		comp_prob.set_rotate_flag(prob.rotatable());
		comp_prob.set_shift_flag(true);

		std::vector<double>::const_iterator sfit = prob.shift_val_vec().begin(), seit = prob.shift_val_vec().begin() + prob.dim();
		std::vector<CProblem::RNV>::const_iterator rfit = prob.rotate_matrix().begin(), reit = prob.rotate_matrix().begin() + prob.dim();

		comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
		comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
		fit[0] = Rastrihin(gene_vec, comp_prob);

		sfit = seit; seit += prob.dim();
		rfit = reit; reit += prob.dim();

		comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
		comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
		fit[1] = Griewangk(gene_vec, comp_prob);
		fit[1] = 1000 * fit[1] / 100;

		sfit = seit; seit += prob.dim();
		rfit = reit; reit += prob.dim();

		comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
		comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
		fit[2] = Modified_Schwefel(gene_vec, comp_prob);

		return cf_cal(gene_vec, prob.shift_val_vec(), fit, delta, bias);
	}

	double cf04(const Individual::GeneVec& gene_vec, const CProblem& prob) /* composition function 2 */
	{
		int cf_num = 4;
		std::vector<double> fit(4);
		std::vector<double> delta = { 10, 20, 30, 40 };
		std::vector<double> bias = { 0, 100, 200, 300 };

		CProblem comp_prob = prob;
		comp_prob.set_rotate_flag(prob.rotatable());
		comp_prob.set_shift_flag(true);

		std::vector<double>::const_iterator sfit = prob.shift_val_vec().begin(), seit = prob.shift_val_vec().begin() + prob.dim();
		std::vector<CProblem::RNV>::const_iterator rfit = prob.rotate_matrix().begin(), reit = prob.rotate_matrix().begin() + prob.dim();

		comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
		comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
		fit[0] = 1000 * Ackley(gene_vec, comp_prob) / 100;

		sfit = seit; seit += prob.dim();
		rfit = reit; reit += prob.dim();

		comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
		comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
		fit[1] = 10000 * Ellips(gene_vec, comp_prob) / 1e+10;

		sfit = seit; seit += prob.dim();
		rfit = reit; reit += prob.dim();

		comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
		comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
		fit[2] = 1000 * Griewangk(gene_vec, comp_prob) / 100;

		sfit = seit; seit += prob.dim();
		rfit = reit; reit += prob.dim();

		comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
		comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
		fit[3] = Rastrihin(gene_vec, comp_prob);

		return cf_cal(gene_vec, prob.shift_val_vec(), fit, delta, bias);
	}

	double cf05(const Individual::GeneVec& gene_vec, const CProblem& prob) /* composition function 3 */
	{
		int cf_num = 5;
		std::vector<double> fit(5);
		std::vector<double> delta = { 10, 20, 30, 40, 50 };
		std::vector<double> bias = { 0, 100, 200, 300, 400 };

		CProblem comp_prob = prob;
		comp_prob.set_rotate_flag(prob.rotatable());
		comp_prob.set_shift_flag(true);

		std::vector<double>::const_iterator sfit = prob.shift_val_vec().begin(), seit = prob.shift_val_vec().begin() + prob.dim();
		std::vector<CProblem::RNV>::const_iterator rfit = prob.rotate_matrix().begin(), reit = prob.rotate_matrix().begin() + prob.dim();

		comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
		comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
		fit[0] = 10000 * Rastrihin(gene_vec, comp_prob) / 1e+3;

		sfit = seit; seit += prob.dim();
		rfit = reit; reit += prob.dim();

		comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
		comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
		fit[1] = 1000 * Happy_Cat(gene_vec, comp_prob) / 1e+3;

		sfit = seit; seit += prob.dim();
		rfit = reit; reit += prob.dim();

		comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
		comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
		fit[2] = 1000 * Ackley(gene_vec, comp_prob) / 100;

		sfit = seit; seit += prob.dim();
		rfit = reit; reit += prob.dim();

		comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
		comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
		fit[3] = 10000 * Discus(gene_vec, comp_prob) / 1e+10;

		sfit = seit; seit += prob.dim();
		rfit = reit; reit += prob.dim();

		comp_prob.set_shift_vector(CProblem::RNV(sfit, seit));
		comp_prob.set_rotate_matrix(std::vector<CProblem::RNV>(rfit, reit));
		fit[4] = Rosenbrock(gene_vec, comp_prob);

		return cf_cal(gene_vec, prob.shift_val_vec(), fit, delta, bias);
	}
}