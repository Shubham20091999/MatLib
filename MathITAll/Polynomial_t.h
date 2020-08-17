#ifndef Poly_t
#define Poly_t

#include"Polynomial.h"

template<typename T>
Polynomial::Polynomial(std::initializer_list<T> args) :
	m_coeffs(std::vector<std::complex<double>>(args.size()))
{
	debug("Construction of Polynomial--(3) ");
	unsigned int i = 0;
	for (T x : args)
	{
		m_coeffs[i] = (std::complex<double>)x;
		i++;
	}
}

template<typename T>
Polynomial::Polynomial(std::initializer_list<std::initializer_list<T>> args) :
	m_coeffs(std::vector<std::complex<double>>(args.size(), 0))
{
	debug("Construction of Polynomial--(3) ");
	unsigned int i = 0;
	for (auto x : args)
	{
		if (x.size() > 2)
			error("Invalid Conversion/There should be only two values one corresponding to real and other to imaginary");
		else if (x.size() == 2)
			m_coeffs[i] = std::complex<double>(*x.begin(), *(x.end() - 1));
		else if (x.size() == 1)
			m_coeffs[i] = std::complex<double>(*x.begin());
		i++;
	}
}

template<typename T>
std::complex<double> Polynomial::operator()(const T& x) const
{
	std::complex<double> a = 0;
	unsigned int sz = m_coeffs.size();
	for (unsigned int i = 0; i < sz; i++)
		a = a * std::complex<double>(x) + m_coeffs[i];
	return a;
}

template<typename T>
Polynomial::Polynomial(unsigned int deg, T* ls) :
	m_coeffs(std::vector < std::complex<double>>(deg + 1, 0))
{
	debug("Construction of Polynomial--(2) ");

	for (unsigned int i = 0; i < deg + 1; i++)
		m_coeffs[i] = (double)ls[i];
}

template<typename T>
Polynomial Polynomial::operator*(const T& c) const
{
	Polynomial n(m_coeffs.size() - 1);
	for (unsigned int i = 0; i < m_coeffs.size(); i++)
		n.m_coeffs[i] = m_coeffs[i] * std::complex<double>(c);
	return n;
}

template<typename T>
void Polynomial::operator*=(const T& c)
{
	for (unsigned int i = 0; i < m_coeffs.size(); i++)
		m_coeffs[i] *= std::complex<double>(c);
}

template<typename T>
void Polynomial::operator/=(const T& c)
{
	*this *= 1.0 / c;
}

template<typename T>
Polynomial Polynomial::operator/(const T& c) const
{
	Polynomial n(Get_Degree());
	for (unsigned int i = 0; i < m_coeffs.size(); i++)
	{
		n.m_coeffs[i] = m_coeffs[i] / std::complex<double>(c);
	}
	return n;
}



template<typename T>
void Polynomial::DivideBySolution(const T& c)
{
	unsigned int deg = m_coeffs.size() - 1;

	Polynomial t(deg - 1);
	std::complex<double> sum = 0;
	for (unsigned int i = 0; i < deg; i++)
	{
		sum = std::complex<double>(c) * sum + m_coeffs[i];
		t.m_coeffs[i] = sum;

	}
	m_coeffs = std::move(t.m_coeffs);
}


template<typename T>
Polynomial operator*(const T& other, const Polynomial& A)
{
	return A * other;
}

#endif // !Poly_t
