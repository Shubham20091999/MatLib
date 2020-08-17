#include"Polynomial.h"

//Object Operations
void Polynomial::Differentiate(unsigned int n)
{
	int sz = m_coeffs.size() - n - 1;
	if (sz < 0)
	{
		m_coeffs.resize(1);
		m_coeffs[0] = 0;
		return;
	}
	m_coeffs.resize(sz + 1);

	for (int i = 0; i < sz + 1; i++)
	{
		m_coeffs[i] = double(Permutations(sz - i + n, sz - i)) * m_coeffs[i];
	}
}

void Polynomial::Integrate(double x, double y)
{
	unsigned int sz = m_coeffs.size();

	if (sz != 1 || m_coeffs[0] != 0.0)
	{
		m_coeffs.resize(sz + 1);
		for (unsigned int i = 0; i < sz; i++)
			m_coeffs[i] = m_coeffs[i] / double(sz - i);
		m_coeffs[sz] = y - (*this)(x);
		return;
	}
	m_coeffs[0] = y;
}


std::complex<double> Polynomial::FindSol_N(const double& StopCriterion) const
{
	Polynomial derivative = Polynomial::Differentiate(*this, 1);
	std::complex<double> c(0.75812488723559, 0.95476369724224);
	std::complex<double> a=c;
	std::complex<double> a_previous;
	unsigned int change=0;
	while (abs(a.real() - a_previous.real()) > StopCriterion && abs(a.imag() - a_previous.imag()) > StopCriterion)
	{
	
		a_previous = a;
		a -= (*this)(a) / derivative(a);
		if (isnan(a.real()) || isnan(a.imag()))
		{
			change += 1;
			a = pow(c, change);
		}
	}
	return a;
}


std::vector<std::complex<double>> Polynomial::FindAllSol_A(bool safe, const double& StopCriterion) const
{
	unsigned int deg = this->Get_Degree();
	Polynomial derivative = Polynomial::Differentiate(*this, 1);
	std::vector<std::complex<double>> a=FindAllSol_N(StopCriterion);
	std::vector<std::complex<double>> a_pre(deg);

	for(unsigned int i=0;i<100;i++)
	{
		bool done = true;

		for(unsigned int i=0;i<deg;i++)
		{
			if(safe)
				if (isnan(a[i].real()) || isnan(a[i].imag()))
					a[i] = std::complex<double>(1, 1);
				  
			auto deno=(*this)(a[i]);
			
			if (deno == 0.0)
			{
				continue;
			}
			std::complex<double> tmp = derivative(a[i]) / deno;

			bool Broken = false;
			for (unsigned int j = 0; j < deg; j++)
			{
				if(j != i)
				{
					if(safe)
						if (a[i] - a[j] == 0.0)
						{
							Broken = true;
							break;
						}
					tmp -= 1.0 / (a[i] - a[j]);
				}
			}
			if (!Broken)
			{
				if (tmp != 0.0)
					a[i] -= 1.0 / tmp;
			}
			
			if (abs(a[i] - a_pre[i]) > 10e-15)
				done = false;
		}
		if(done)
			break;

		a_pre=a;
	}

	return a;
}

std::vector<std::complex<double>> Polynomial::FindAllSol_N(const double& StopCriterion) const
{
	unsigned int deg = this->Get_Degree();
	Polynomial p = *this;
	std::vector<std::complex<double>> a(deg);

	for(unsigned int i=0;i<deg;i++)
	{
		std::complex<double> s = p.FindSol_N(StopCriterion);
		a[i] = s;
		p.DivideBySolution(s);
	}

	return a;
}

//Static Operations
Polynomial Polynomial::Differentiate(const Polynomial& alt, unsigned int n = 1)
{
	int sz = alt.m_coeffs.size() - n - 1;
	if (sz < 0)
	{
		Polynomial ret(1);
		ret.m_coeffs[0] = 0;
		return ret;
	}

	Polynomial ret(sz);

	for (int i = 0; i < sz + 1; i++)
	{
		ret.m_coeffs[i] = double(Permutations(sz - i + n, sz - i)) * alt.m_coeffs[i];
	}
	return ret;
}

Polynomial Polynomial::Integrate(const Polynomial& alt, double x, double y)
{
	unsigned int sz = alt.m_coeffs.size();

	Polynomial a(sz);
	if (sz != 1 || alt.m_coeffs[0] != 0.0)
	{
		for (unsigned int i = 0; i < sz; i++)
			a.m_coeffs[i] = alt.m_coeffs[i] / double(sz - i);
		a.m_coeffs[sz] = y - a(x);
		return a;
	}
	a.m_coeffs[0] = y;
	return a;
}

std::vector<std::complex<double>> Polynomial::FindAllSol_N(const Polynomial& A, const double& StopCriterion)
{
	return A.FindAllSol_N(StopCriterion);
}

std::vector<std::complex<double>> Polynomial::FindAllSol_A(const Polynomial& A,bool safe, const double& StopCriterion)
{
	return A.FindAllSol_A(safe,StopCriterion);
}


