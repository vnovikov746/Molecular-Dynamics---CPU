/*
* Molecular Dynamics Project.
* Outhor: Vladimir Novikov.
*/

#include "Equi.h"

int writeParticlesInput(configurations *config)
{
	string line;
	config->SI_PARTICLES = 0;
	config->XE_PARTICLES = 0;
	int pos = 0;
	int posInLine_1 = 0;
	int posInLine_2 = 0;
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;

	ifstream conditionsFile;
	conditionsFile.open(*config->CONDITIONS_FILE);
	if (!conditionsFile.good())
	{
		cout << "Can not open conditions file" << endl;
		return EXIT_FAILURE;
	}

	ofstream MDinputFile;
	MDinputFile.open(*config->INPUT_FILE);
	if (!MDinputFile.good())
	{
		cout << "Can not open simulations input file" << endl;
		return EXIT_FAILURE;
	}

	MDinputFile.precision(19);

	while (getline(conditionsFile, line))
	{
		if (line.find("//") && !line.empty())// If comment is at the start of the line, find will return 0.
		{
			if (line.find("Si") == 0)
			{
				pos = line.find_first_of('=');
				posInLine_1 = line.find_first_of('x');
				posInLine_2 = line.find_first_of('y');
				x = atof(line.substr(posInLine_1 + 2, posInLine_2).c_str());
				posInLine_1 = line.find_first_of('z');
				y = atof(line.substr(posInLine_2 + 2, posInLine_1).c_str());
				z = atof(line.substr(posInLine_1 + 2, line.size()).c_str());

				int x_size = (config->SI_LENGTH);
				int y_size = (config->SI_WIDTH);
				int z_size = (config->SI_HEIGHT);

				double movedX;
				double movedY;
				double movedZ;

				MDinputFile << "\n//------------------------------------------------\n";
				MDinputFile << "//" << line << "\n";//0 0 0
				MDinputFile << "//------------------------------------------------\n";

				double check1 = 0.0;
				double check2 = 0.0;
				double check3 = 0.0;

				for (int i = 0; i <= (x_size); i++)
				{
					movedX = (x + i);
					if (x == 0.0)
					{
						check1 = 0.75;
					}
					else if (x == 0.25)
					{
						check1 = 0.5;
					}
					else if (x == 0.5)
					{
						check1 = 0.25;
					}
					else if (x == 0.75)
					{
						check1 = 0.0;
					}
					for (int j = 0; j <= (y_size); j++)
					{
						movedY = (y + j);
						if (y == 0.0)
						{
							check2 = 0.75;
						}
						else if (y == 0.25)
						{
							check2 = 0.5;
						}
						else if (y == 0.5)
						{
							check2 = 0.25;
						}
						else if (y == 0.75)
						{
							check2 = 0.0;
						}
						for (int k = 0; k <= (z_size); k++)
						{
							movedZ = (z + k);
							if (z == 0.0)
							{
								check3 = 0.75;
							}
							else if (z == 0.25)
							{
								check3 = 0.5;
							}
							else if (z == 0.5)
							{
								check3 = 0.25;
							}
							else if (z == 0.75)
							{
								check3 = 0.0;
							}
							if ((movedZ + check3) <= z_size && (movedX + check1) <= x_size && (movedY + check2) <= y_size)
							{
								MDinputFile << "Si   x=" << movedX*(space_Si) << "  y=" << movedY*(space_Si) << "  z=" << movedZ*(space_Si) << "\n";
								config->SI_PARTICLES++;
							}
						}
					}
				}
			}
			else if (line.find("Xe") == 0)
			{
				pos = line.find_first_of('=');
				posInLine_1 = line.find_first_of('x');
				posInLine_2 = line.find_first_of('y');
				x = atof(line.substr(posInLine_1 + 2, posInLine_2).c_str());
				posInLine_1 = line.find_first_of('z');
				y = atof(line.substr(posInLine_2 + 2, posInLine_1).c_str());
				z = atof(line.substr(posInLine_1 + 2, line.size()).c_str());

				int x_size = (config->XE_LENGTH);
				int y_size = (config->XE_WIDTH);
				int z_size = (config->XE_HEIGHT);

				double movedX;
				double movedY;
				double movedZ;

				MDinputFile << "\n//------------------------------------------------\n";
				MDinputFile << "//" << line << "\n";//0 0 0
				MDinputFile << "//------------------------------------------------\n";

				double check1 = 0.0;
				double check2 = 0.0;
				double check3 = 0.0;
				for (int i = 0; i <= (x_size); i++)
				{
					movedX = (x + i);
					if (x == 0.0)
					{
						check1 = 0.75;
					}
					else if (x == 0.5)
					{
						check1 = 0.25;
					}
					for (int j = 0; j <= (y_size); j++)
					{
						movedY = (y + j);
						if (y == 0.0)
						{
							check2 = 0.75;
						}
						else if (y == 0.5)
						{
							check2 = 0.25;
						}
						for (int k = 0; k <= (z_size); k++)
						{
							movedZ = (z + k);
							if (z == 0.0)
							{
								check3 = 0.75;
							}
							else if (z == 0.5)
							{
								check3 = 0.25;
							}
							if ((movedZ + check3) <= z_size && (movedX + check1) <= x_size && (movedY + check2) <= y_size)
							{
								MDinputFile << "Xe   x=" << movedX*space_Xe + 0.25*config->SI_LENGTH*space_Si << "  y=" << movedY*space_Xe + 0.25*config->SI_WIDTH*space_Si << "  z=" << movedZ*space_Xe + config->SI_HEIGHT*space_Si + config->LA_SPACE << "\n";
								config->XE_PARTICLES++;
							}
						}
					}
				}
			}
		}
	}
	config->SI_HEIGHT *= (int)space_Xe;
	config->SI_WIDTH *= (int)space_Xe;
	config->SI_LENGTH *= (int)space_Xe;
	config->XE_HEIGHT *= (int)space_Xe;
	config->XE_WIDTH *= (int)space_Xe;
	config->XE_LENGTH *= (int)space_Xe;
	MDinputFile.close();
	conditionsFile.close();
	return EXIT_SUCCESS;
}

/***********************************************************************

Inverse of the cumulative normal distribution

Implementation:

***********************************************************************/

double Erf(double x)
//Erf(x) = 2/sqrt(pi)*integral(exp(-t^2),t,0,x)
//       = 2/sqrt(pi)*[x - x^3/3 + x^5/5*2! - x^7/7*3! + ...]
//       = 1-Erfc(x)
{
	static const double two_sqrtpi = 1.128379167095512574;        // 2/sqrt(pi)
	if (fabs(x) > 2.2) {
		return 1.0 - Erfc(x);        //use continued fraction when fabs(x) > 2.2
	}
	double sum = x, term = x, xsqr = x*x;
	int j = 1;
	do {
		term *= xsqr / j;
		sum -= term / (2 * j + 1);
		++j;
		term *= xsqr / j;
		sum += term / (2 * j + 1);
		++j;
	} while (fabs(term) / sum > rel_error);
	return two_sqrtpi*sum;
}

double Erfc(double x)
//Erfc(x) = 2/sqrt(pi)*integral(exp(-t^2),t,x,inf)
//        = exp(-x^2)/sqrt(pi) * [1/x+ (1/2)/x+ (2/2)/x+ (3/2)/x+ (4/2)/x+ ...]
//        = 1-Erf(x)
//expression inside [] is a continued fraction so '+' means add to denominator only
{
	static const double one_sqrtpi = 0.564189583547756287;        // 1/sqrt(pi)
	if (fabs(x) < 2.2) {
		return 1.0 - Erf(x);        //use series when fabs(x) < 2.2
	}
	if (x<0) {               //continued fraction only valid for x>0
		return 2.0 - Erfc(-x);
	}
	double a = 1, b = x;                //last two convergent numerators
	double c = x, d = x*x + 0.5;          //last two convergent denominators
	double q1, q2 = b / d;             //last two convergents (a/c and b/d)
	double n = 1.0, t;
	do {
		t = a*n + b*x;
		a = b;
		b = t;
		t = c*n + d*x;
		c = d;
		d = t;
		n += 0.5;
		q1 = q2;
		q2 = b / d;
	} while (fabs(q1 - q2) / q2 > rel_error);
	return one_sqrtpi*exp(-x*x)*q2;
}

double Phi(double z)
{
	return 0.5 * Erfc(-C1sqrt2*z);
}

/***********************************************************************

Inverse of the error function Erf.

Implementation: Inversion by Newton iteration of Erf(x).
The initial value x0 = 0.
For |z| <= 0.84 (=Erf(1)) at most 4 iterations are necessary.

***********************************************************************/

double InvErfSmall(const double z)
{
	/* f(x)   = Erf(x) - z   */
	/* f'(x)  = c*exp(-x*x)  */
	/* f''(x) = -2 f'(x)     */
	double c = C2sqrtPi;
	double f = -z, f1 = c;
	double q = f / f1, x = -q, x0 = 0;

	while (fabs(x - x0) > 1e-12 && fabs(f) > 1e-14) {
		/* Newton 2nd order: x <- x - f/f'(1 + f*f''/(2 f'^2)) */
		x0 = x;
		f = Erf(x) - z;
		f1 = c*exp(-x*x);
		q = f / f1;
		x -= q*(1 - x*q);  /* Newton Step 2nd order */
	}

	return x;
}

/***********************************************************************

lambert_w2 is that "Dilbert lambda" (inverse x e^x^2)
If you don't have W functions, you can probably get away with

log(log(x))
lambert_w2(x) ~ sqrt(log(x)) - -------------- .
4 sqrt(log(x))

***********************************************************************/

double lambertW2(const double z)
{
	/* "Dilbert lambda" (inverse x e^x^2) approximation */
	double logz = log(z);
	double slz = sqrt(logz);
	return slz - 0.25*log(logz) / slz;
}


/***********************************************************************

Inverse of the error function Erfc.

Implementation: Inversion by Newton iteration of Erfc(sqrt(log(x))).
The initial value is computed via lambertW2.
For z < 0.25 at most 4 iterations are necessary.

***********************************************************************/

double InvErfcSmall(const double z)
{
	/* f(x)   = Erfc(sqrt(log(x))) - z   */
	/* f'(x)  = 1/(c x^2 sqrt(log(x)))   */
	/* f''(x) = -c*x*f'(x)^2*(2*sqrt(log(x))+1/(2*sqrt(log(x))))   */
	double c = CsqrtPi;
	double f = 1, f1i = 0;
	double a = lambertW2(1 / (c*z));
	double x = Erfc(a) / (c * a * z*z), x0 = 0;

	while (fabs(x - x0) > 1e-12*x && fabs(f) > 1e-15*z) {
		/* Newton 2nd order: x <- x - f/f'(1 + f*f''/(2 f'^2)) */
		double slx;
		x0 = x;
		slx = sqrt(log(x));
		f = z - Erfc(slx);
		f1i = c * x * x * slx;
		x -= f*f1i*(1 - c*f*x*(slx + 0.25 / slx));       /* Newton Step */
	}

	return  sqrt(log(x));
}

/***********************************************************************

Inverse of the error function Erf.

Implementation:
For small and big values two differnt approximations are used.

Return values:
If x is +1, InvErf() returns +INFINITY.
If x is -1, InvErf() returns -INFINITY.
If x is not in the range [-1, 1],
InvErf() returns NaN and sets errno to EDOM.

***********************************************************************/

double inverf(const double z)
{
	double az = fabs(z);
	double x = 0;
	if (az >= 1)
	{
		double huge = DBL_MAX;
		double inf = 2 * huge;
		/* Infinity of NaN */
		if (z == 1)
			return  inf; /* +infinity */
		if (z == -1)
			return -inf; /* -infinity */
		/* out of range, set errno=EDOM, return NaN */
		errno = EDOM;
		return inf - inf;
	}

	if (az < 0.8125)          /* -13/16 < z < 13/16 */
		x = InvErfSmall(z);
	else if (z > 0)           /* z >=  13/16 */
		x = InvErfcSmall(1 - z);
	else                      /* z <= -13/16 */
		x = -InvErfcSmall(z + 1);
	return x;
}

double inverfc(const double z)
{
	double x = 0;
	if (z <= 0 || z >= 2)
	{
		double huge = DBL_MAX;
		double inf = 2 * huge;
		/* Infinity of NaN */
		if (z == 0)
			return  inf; /* +infinity */
		if (z == 2)
			return -inf; /* -infinity */
		/* out of range, set errno=EDOM, return NaN */
		errno = EDOM;
		return inf - inf;
	}

	if (z <= 0.1875)          /* z <=   3/16 */
		x = InvErfcSmall(z);
	else if (z < 1.8125)      /* 3/16 < z < 29/16 */
		x = InvErfSmall(1 - z);
	else                      /* z >=  29/16 */
		x = -InvErfcSmall(2 - z);
	return x;
}

void chooseVelocities(particleStruct *particles, int listSize, double temp, configurations *config, double mass)
{
	double cmvx = 0;
	double cmvy = 0;
	double cmvz = 0;
	double totalm = 0;
	double t = 0;

	for (int i = 0; i < listSize; i++)
	{
		double y = ((double)rand() + 1) / (double)(RAND_MAX + 2.0);

		y = 2.0*y - 1.0;

		y = inverf(y);

		particles[i].velocity.x = y*sqrt(2.0*temp / mass);

		y = ((double)rand() + 1) / (double)(RAND_MAX + 2.0);

		y = 2.0*y - 1.0;

		y = inverf(y);

		particles[i].velocity.y = y*sqrt(2.0*temp / mass);

		y = ((double)rand() + 1) / (double)(RAND_MAX + 2.0);

		y = 2.0*y - 1.0;

		y = inverf(y);

		particles[i].velocity.z = y*sqrt(2.0*temp / mass);

		cmvx = cmvx + mass*particles[i].velocity.x;
		cmvy = cmvy + mass*particles[i].velocity.y;

		cmvz = cmvz + mass*particles[i].velocity.z;
		totalm += mass;
	}

	cmvx = cmvx / totalm;
	cmvy = cmvy / totalm;
	cmvz = cmvz / totalm;

	while ((fabs(cmvx) > 1.0e-15) || (fabs(cmvy) > 1.0e-15) || (fabs(cmvz) > 1.0e-15))
	{
		for (int i = 0; i < listSize; i++)
		{
			particles[i].velocity.x = particles[i].velocity.x - cmvx;
			particles[i].velocity.y = particles[i].velocity.y - cmvy;
			particles[i].velocity.z = particles[i].velocity.z - cmvz;
		}

		cmvx = 0.0;
		cmvy = 0.0;
		cmvz = 0.0;

		for (int i = 0; i < listSize; i++)
		{
			cmvx = cmvx + particles[i].velocity.x*mass;
			cmvy = cmvy + particles[i].velocity.y*mass;
			cmvz = cmvz + particles[i].velocity.z*mass;
		}

		cmvx = cmvx / totalm;
		cmvy = cmvy / totalm;
		cmvz = cmvz / totalm;
	}

	for (int i = 0; i < listSize; i++)
	{
		t += mass*(particles[i].velocity.x*particles[i].velocity.x
			+ particles[i].velocity.y*particles[i].velocity.y
			+ particles[i].velocity.z*particles[i].velocity.z);
	}

	t = t / (3.0*(double)(listSize));
	double scalev = sqrt(temp / t);

	for (int i = 0; i < listSize; i++)
	{
		particles[i].velocity.x = particles[i].velocity.x*scalev;
		particles[i].velocity.y = particles[i].velocity.y*scalev;
		particles[i].velocity.z = particles[i].velocity.z*scalev;
	}
}
