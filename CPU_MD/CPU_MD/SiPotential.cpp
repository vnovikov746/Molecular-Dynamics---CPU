

#include "SiPotential.h"

//-------------------- force between two Si particles ---------------------//
double f2_derivative_of_rij_tag(double r_ij_tag)
{
	double first = -4*B_Si*(1.0/(r_ij_tag*r_ij_tag*r_ij_tag*r_ij_tag*r_ij_tag));
	double second = ((B_Si*(1.0/(r_ij_tag*r_ij_tag*r_ij_tag*r_ij_tag)))-1)*(1.0/((r_ij_tag-a_Si)*(r_ij_tag-a_Si)));
	double print = first-second;

	double r_ij_tag_minus_a = r_ij_tag - a_Si;//r'ij-a
	double r_ij_tag_minus_a2 = r_ij_tag_minus_a*r_ij_tag_minus_a;

	double r_ij_tag_minus_a_in_mOne = (1.0/r_ij_tag_minus_a);//(r'ij-a)^(-1)
	double r_ij_tag_minus_a_in_mTwo = (1.0/(r_ij_tag_minus_a2));//(r'ij-a)^(-2)

	double r_ij_tag2 = r_ij_tag*r_ij_tag;
	double r_ij_tag4 = r_ij_tag2*r_ij_tag2;
	double r_ij_tag5 = r_ij_tag4*r_ij_tag;

	double r_ij_tag_in_mFive = (1.0/(r_ij_tag5));//r'ij^(-5)
	double r_ij_tag_in_mFour = (1.0/(r_ij_tag4));//r'ij^(-4)
	double expression = B_Si * r_ij_tag_in_mFour;
	expression = expression - 1.0;//(B*r'ij^(-4) - 1)

	double exponent = exp(r_ij_tag_minus_a_in_mOne);

	double f2_derivative_part_1 = -4.0 * B_Si * r_ij_tag_in_mFive;
	double f2_derivative_part_2 = expression * r_ij_tag_minus_a_in_mTwo;

	return A_Si*exponent*(f2_derivative_part_1 - f2_derivative_part_2);
}

double v2_derivative_of_rix(real3 i, real3 j, double r_ij)
{
	if(r_ij/sigma_Si == pow(2.0,1.0/6.0))
	{
		return 0;
	}

	double f2_derivative = f2_derivative_of_rij_tag(r_ij/sigma_Si);

	f2_derivative = f2_derivative * (epsilon_Si/sigma_Si);//v2 derivative of distance
	double dist_x = (i.x - j.x);
	dist_x = dist_x / (r_ij);

	double v2_derivative = f2_derivative * dist_x;

	return v2_derivative;
}

double v2_derivative_of_riy(real3 i, real3 j, double r_ij)
{
	if(r_ij/sigma_Si == pow(2.0,1.0/6.0))
	{
		return 0;
	}

	double f2_derivative = f2_derivative_of_rij_tag(r_ij/sigma_Si);

	f2_derivative = f2_derivative * (epsilon_Si/sigma_Si);
	double dist_y = i.y - j.y;
	dist_y = dist_y / (r_ij);

	double v2_derivative = f2_derivative * dist_y;

	return v2_derivative;
}

double v2_derivative_of_riz(real3 i, real3 j, double r_ij)
{
	if(r_ij/sigma_Si == pow(2.0,1.0/6.0))
	{
		return 0;
	}

	double f2_derivative = f2_derivative_of_rij_tag(r_ij/sigma_Si);

	f2_derivative = f2_derivative  * (epsilon_Si/sigma_Si);
	double dist_z = i.z - j.z;
	dist_z = dist_z / (r_ij);

	double v2_derivative = f2_derivative * dist_z;

	return v2_derivative;
}
//----------------------------------------------------------------------------//

//-------------------- potential between two Si particles ---------------------//
double f2(double r_ij_tag)
{
	if(r_ij_tag >= a_Si)
	{
		return 0;
	}
	double r_ij_tag_minus_a_in_mOne = r_ij_tag - a_Si;
	r_ij_tag_minus_a_in_mOne = (1.0/r_ij_tag_minus_a_in_mOne);

	double exponent = exp(r_ij_tag_minus_a_in_mOne);
	double r_ij_tag2 = r_ij_tag*r_ij_tag;
	double r_ij_tag4 = r_ij_tag2*r_ij_tag2;

	double expression = (1.0/(r_ij_tag4));
	expression *= B_Si;
	expression -= 1.0;

	return A_Si*expression*exponent;
}

double v2(double r_ij_tag)
{
	if(r_ij_tag == pow(2.0,1.0/6.0))
	{
		return -epsilon_Si;
	}
	return f2(r_ij_tag)*epsilon_Si;
}
//----------------------------------------------------------------------------//

//------------------------ force between three Si particles -------------------//
double hi_derivative_of_rij_tag(double r_ij_tag, double r_ik_tag, double r_jk_tag)
{
	double cosJik_plus_oneThird = ((r_ij_tag*r_ij_tag + r_ik_tag*r_ik_tag - r_jk_tag*r_jk_tag)/(2.0 * r_ij_tag * r_ik_tag)) + (1.0/3.0);

	double r_ij_tag_minus_a = r_ij_tag - a_Si;
	double r_ij_tag_minus_a_in_mOne = (1.0/r_ij_tag_minus_a) * gama_Si;
	double r_ij_tag_minus_a_in_mTwo = (1.0/(r_ij_tag_minus_a*r_ij_tag_minus_a)) * gama_Si;

	double r_ik_tag_minus_a_in_mOne = (1.0/r_ik_tag - a_Si) * gama_Si;

	double exponent = exp(r_ij_tag_minus_a_in_mOne+r_ik_tag_minus_a_in_mOne);

	double expression = (r_ij_tag*r_ij_tag - r_ik_tag*r_ik_tag + r_jk_tag*r_jk_tag) / (r_ij_tag*r_ij_tag * r_ik_tag);

	expression -= (r_ij_tag_minus_a_in_mTwo*cosJik_plus_oneThird);

	return lamda_Si*exponent*cosJik_plus_oneThird*expression;
}

double hi_derivative_of_rik_tag(double r_ij_tag, double r_ik_tag, double r_jk_tag)
{
	double cosJik_plus_oneThird = ((r_ij_tag*r_ij_tag + r_ik_tag*r_ik_tag - r_jk_tag*r_jk_tag)/(2.0 * r_ij_tag * r_ik_tag)) + (1.0/3.0);

	double r_ij_tag_minus_a_in_mOne = (1.0/(r_ij_tag - a_Si)) * gama_Si;

	double r_ik_tag_minus_a = r_ik_tag - a_Si;
	double r_ik_tag_minus_a_in_mOne = (1.0/r_ik_tag_minus_a) * gama_Si;
	double r_ik_tag_minus_a_in_mTwo = (1.0/(r_ik_tag_minus_a*r_ik_tag_minus_a)) * gama_Si;

	double exponent = exp(r_ij_tag_minus_a_in_mOne+r_ik_tag_minus_a_in_mOne);

	double expression = (r_ik_tag*r_ik_tag - r_ij_tag*r_ij_tag + r_jk_tag*r_jk_tag) / (r_ik_tag*r_ik_tag * r_ij_tag);

	expression -= (r_ik_tag_minus_a_in_mTwo*cosJik_plus_oneThird);

	return lamda_Si*exponent*cosJik_plus_oneThird*expression;
}

double hj_derivative_of_rij_tag(double r_ij_tag, double r_ik_tag, double r_jk_tag)
{
	double cosIjk_plus_oneThird = ((r_ij_tag*r_ij_tag + r_jk_tag*r_jk_tag - r_ik_tag*r_ik_tag)/(2.0 * r_ij_tag * r_jk_tag)) + (1.0/3.0);

	double r_ij_tag_minus_a = r_ij_tag - a_Si;
	double r_ij_tag_minus_a_in_mOne = (1.0/r_ij_tag_minus_a) * gama_Si;
	double r_ij_tag_minus_a_in_mTwo = (1.0/(r_ij_tag_minus_a*r_ij_tag_minus_a)) * gama_Si;

	double r_jk_tag_minus_a_in_mOne = (1.0/(r_jk_tag - a_Si)) * gama_Si;

	double exponent = exp(r_ij_tag_minus_a_in_mOne+r_jk_tag_minus_a_in_mOne);

	double expression = (r_ij_tag*r_ij_tag - r_jk_tag*r_jk_tag + r_ik_tag*r_ik_tag) / (r_ij_tag*r_ij_tag * r_jk_tag);

	expression -= (r_ij_tag_minus_a_in_mTwo*cosIjk_plus_oneThird);

	return lamda_Si*exponent*cosIjk_plus_oneThird*expression;
}

double hj_derivative_of_rik_tag(double r_ij_tag, double r_ik_tag, double r_jk_tag)
{
	double cosIjk_plus_oneThird = ((r_ij_tag*r_ij_tag + r_jk_tag*r_jk_tag - r_ik_tag*r_ik_tag)/(2.0 * r_ij_tag * r_jk_tag)) + (1.0/3.0);

	double r_ij_tag_minus_a_in_mOne = (1.0/(r_ij_tag - a_Si)) * gama_Si;
	double r_jk_tag_minus_a_in_mOne = (1.0/(r_jk_tag - a_Si)) * gama_Si;

	double exponent = exp(r_ij_tag_minus_a_in_mOne+r_jk_tag_minus_a_in_mOne);

	double expression = (-r_ik_tag) / (r_ij_tag * r_jk_tag);

	return lamda_Si*exponent*2.0*cosIjk_plus_oneThird*expression;
}

double hk_derivative_of_rij_tag(double r_ij_tag, double r_ik_tag, double r_jk_tag)
{
	double cosIkj_plus_oneThird = ((r_ik_tag*r_ik_tag + r_jk_tag*r_jk_tag - r_ij_tag*r_ij_tag)/(2 * r_ik_tag * r_jk_tag)) + (1.0/3.0);

	double r_ik_tag_minus_a_in_mOne = (1.0/(r_ik_tag - a_Si)) * gama_Si;

	double r_jk_tag_minus_a_in_mOne = (1.0/(r_jk_tag - a_Si)) * gama_Si;

	double exponent = exp(r_ik_tag_minus_a_in_mOne+r_jk_tag_minus_a_in_mOne);

	double expression = (-r_ij_tag) / (r_ik_tag * r_jk_tag);

	return lamda_Si*exponent*2.0*cosIkj_plus_oneThird*expression;
}

double hk_derivative_of_rik_tag(double r_ij_tag, double r_ik_tag, double r_jk_tag)
{
	double cosIkj_plus_oneThird = ((r_ik_tag*r_ik_tag + r_jk_tag*r_jk_tag - r_ij_tag*r_ij_tag)/(2.0 * r_ik_tag * r_jk_tag)) + (1.0/3.0);

	double r_jk_tag_minus_a_in_mOne = (1.0/(r_jk_tag - a_Si)) * gama_Si;

	double r_ik_tag_minus_a = r_ik_tag - a_Si;
	double r_ik_tag_minus_a_in_mOne = (1.0/r_ik_tag_minus_a) * gama_Si;
	double r_ik_tag_minus_a_in_mTwo = (1.0/(r_ik_tag_minus_a*r_ik_tag_minus_a)) * gama_Si;

	double exponent = exp(r_ik_tag_minus_a_in_mOne+r_jk_tag_minus_a_in_mOne);

	double expression = (r_ik_tag*r_ik_tag - r_jk_tag*r_jk_tag + r_ij_tag*r_ij_tag) / (r_ik_tag*r_ik_tag * r_jk_tag);

	expression -= (r_ik_tag_minus_a_in_mTwo*cosIkj_plus_oneThird);

	return lamda_Si*exponent*cosIkj_plus_oneThird*expression;
}

double f3_derivative_of_rij_tag(double r_ij_tag, double r_ik_tag, double r_jk_tag)
{
	double hi_derivative_of_rij = 0.0;
	double hj_derivative_of_rij = 0.0;
	double hk_derivative_of_rij = 0.0;
	if(r_ij_tag < a_Si && r_ik_tag < a_Si)
	{
		hi_derivative_of_rij = hi_derivative_of_rij_tag(r_ij_tag,r_ik_tag,r_jk_tag);
	}
	if(r_jk_tag < a_Si && r_ij_tag < a_Si)
	{
		hj_derivative_of_rij = hj_derivative_of_rij_tag(r_ij_tag,r_ik_tag,r_jk_tag);
	}
	if(r_jk_tag < a_Si && r_ik_tag < a_Si)
	{
		hk_derivative_of_rij = hk_derivative_of_rij_tag(r_ij_tag,r_ik_tag,r_jk_tag);
	}

	return hi_derivative_of_rij + hj_derivative_of_rij + hk_derivative_of_rij;
}

double f3_derivative_of_rik_tag(double r_ij_tag, double r_ik_tag, double r_jk_tag)
{
	double hi_derivative_of_rik = 0.0;
	double hj_derivative_of_rik = 0.0;
	double hk_derivative_of_rik = 0.0;
	if(r_ik_tag < a_Si && r_ij_tag < a_Si)
	{
		hi_derivative_of_rik = hi_derivative_of_rij_tag(r_ij_tag,r_ik_tag,r_jk_tag);
	}
	if(r_jk_tag < a_Si && r_ij_tag < a_Si)
	{
		hj_derivative_of_rik = hj_derivative_of_rij_tag(r_ij_tag,r_ik_tag,r_jk_tag);
	}
	if(r_jk_tag < a_Si && r_ik_tag < a_Si)
	{
		hk_derivative_of_rik = hk_derivative_of_rij_tag(r_ij_tag,r_ik_tag,r_jk_tag);
	}

	return hi_derivative_of_rik + hj_derivative_of_rik + hk_derivative_of_rik;
}

double v3_derivative_of_rix(real3 i, real3 j, real3 k, double r_ij, double r_ik, double r_jk)
{
	double v3_derived_by_rij = (f3_derivative_of_rij_tag(r_ij/sigma_Si, r_ik/sigma_Si, r_jk/sigma_Si))*(epsilon_Si/sigma_Si);
	double v3_derived_by_rik = (f3_derivative_of_rik_tag(r_ij/sigma_Si, r_ik/sigma_Si, r_jk/sigma_Si))*(epsilon_Si/sigma_Si);

	double dist_ijx = (i.x-j.x);
	double dist_ikx = (i.x-k.x);
	double expression1 = (dist_ijx/(r_ij));
	double expression2 = (dist_ikx/(r_ik));

	return v3_derived_by_rij*expression1 + v3_derived_by_rik*expression2;
}

double v3_derivative_of_riy(real3 i, real3 j, real3 k, double r_ij, double r_ik, double r_jk)
{
	double v3_derived_by_rij = (f3_derivative_of_rij_tag(r_ij/sigma_Si, r_ik/sigma_Si, r_jk/sigma_Si))*(epsilon_Si/sigma_Si);
	double v3_derived_by_rik = (f3_derivative_of_rik_tag(r_ij/sigma_Si, r_ik/sigma_Si, r_jk/sigma_Si))*(epsilon_Si/sigma_Si);

	double dist_ijy = (i.y-j.y);
	double dist_iky = (i.y-k.y);
	double expression1 = (dist_ijy/(r_ij));
	double expression2 = (dist_iky/(r_ik));

	return v3_derived_by_rij*expression1 + v3_derived_by_rik*expression2;
}

double v3_derivative_of_riz(real3 i, real3 j, real3 k, double r_ij, double r_ik, double r_jk)
{
	double v3_derived_by_rij = (f3_derivative_of_rij_tag(r_ij/sigma_Si, r_ik/sigma_Si, r_jk/sigma_Si))*(epsilon_Si/sigma_Si);
	double v3_derived_by_rik = (f3_derivative_of_rik_tag(r_ij/sigma_Si, r_ik/sigma_Si, r_jk/sigma_Si))*(epsilon_Si/sigma_Si);

	double dist_ijz = (i.z-j.z);
	double dist_ikz = (i.z-k.z);
	double expression1 = (dist_ijz/(r_ij));
	double expression2 = (dist_ikz/(r_ik));

	return v3_derived_by_rij*expression1 + v3_derived_by_rik*expression2;
}
//----------------------------------------------------------------------------//


//-------------------- potential between three Si particles -------------------//
double hi(double r_ij_tag, double r_ik_tag, double r_jk_tag)
{
	double cosJik_plus_oneThird = ((r_ij_tag*r_ij_tag + r_ik_tag*r_ik_tag - r_jk_tag*r_jk_tag)/(2.0 * r_ij_tag * r_ik_tag)) + (1.0/3.0);

	double r_ij_tag_minus_a_in_mOne = r_ij_tag - a_Si;
	r_ij_tag_minus_a_in_mOne = (1.0/r_ij_tag_minus_a_in_mOne);
	r_ij_tag_minus_a_in_mOne *= gama_Si;

	double r_ik_tag_minus_a_in_mOne = r_ik_tag - a_Si;
	r_ik_tag_minus_a_in_mOne = (1.0/r_ik_tag_minus_a_in_mOne);
	r_ik_tag_minus_a_in_mOne *= gama_Si;

	double exponent = exp(r_ij_tag_minus_a_in_mOne + r_ik_tag_minus_a_in_mOne);

	return lamda_Si*exponent*cosJik_plus_oneThird*cosJik_plus_oneThird;
}

double hj(double r_ij_tag, double r_ik_tag, double r_jk_tag)
{
	double cosIjk_plus_oneThird = ((r_ij_tag*r_ij_tag + r_jk_tag*r_jk_tag - r_ik_tag*r_ik_tag)/(2.0 * r_ij_tag * r_jk_tag)) + (1.0/3.0);

	double r_ij_tag_minus_a_in_mOne = r_ij_tag - a_Si;
	r_ij_tag_minus_a_in_mOne = (1.0/r_ij_tag_minus_a_in_mOne);
	r_ij_tag_minus_a_in_mOne *= gama_Si;

	double r_jk_tag_minus_a_in_mOne = r_ik_tag - a_Si;
	r_jk_tag_minus_a_in_mOne = (1/r_jk_tag_minus_a_in_mOne);
	r_jk_tag_minus_a_in_mOne *= gama_Si;

	double exponent = exp(r_ij_tag_minus_a_in_mOne + r_jk_tag_minus_a_in_mOne);

	return lamda_Si*exponent*cosIjk_plus_oneThird*cosIjk_plus_oneThird;
}

double hk(double r_ij_tag, double r_ik_tag, double r_jk_tag)
{
	double cosIkj_plus_oneThird = ((r_ik_tag*r_ik_tag + r_jk_tag*r_jk_tag - r_ij_tag*r_ij_tag)/(2.0 * r_ik_tag * r_jk_tag)) + (1.0/3.0);

	double r_ik_tag_minus_a_in_mOne = r_ik_tag - a_Si;
	r_ik_tag_minus_a_in_mOne = (1.0/r_ik_tag_minus_a_in_mOne);
	r_ik_tag_minus_a_in_mOne *= gama_Si;

	double r_jk_tag_minus_a_in_mOne = r_jk_tag - a_Si;
	r_jk_tag_minus_a_in_mOne = (1.0/r_jk_tag_minus_a_in_mOne);
	r_jk_tag_minus_a_in_mOne *= gama_Si;

	double exponent = exp(r_ik_tag_minus_a_in_mOne + r_jk_tag_minus_a_in_mOne);

	return lamda_Si*exponent*cosIkj_plus_oneThird*cosIkj_plus_oneThird;
}

double f3(double r_ij_tag, double r_ik_tag, double r_jk_tag)
{
	double h_i = 0.0;
	double h_j = 0.0;
	double h_k = 0.0;
	if(r_ij_tag < a_Si && r_ik_tag < a_Si)
	{
		h_i = hi(r_ij_tag,r_ik_tag,r_jk_tag);
	}
	if(r_ij_tag < a_Si && r_ik_tag < a_Si)
	{
		h_j = hj(r_ij_tag,r_ik_tag,r_jk_tag);
	}
	if(r_jk_tag < a_Si && r_ik_tag < a_Si)
	{
		h_k = hk(r_ij_tag,r_ik_tag,r_jk_tag);
	}
	return h_i + h_j + h_k;
}

double v3(double r_ij_tag, double r_ik_tag, double r_jk_tag)
{
	return f3(r_ij_tag,r_ik_tag,r_jk_tag)*epsilon_Si;
}
//----------------------------------------------------------------------------//
