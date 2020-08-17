#include"Polynomial.h"

//Constructors
Polynomial::Polynomial(unsigned int deg) :
	m_coeffs(std::vector<std::complex<double>>(deg + 1, 0))
{
	debug("Construction of Polynomial--(1) ");
}

//Destructor
Polynomial::~Polynomial()
{
	debug("Polynomial Destroyed");
}

//Copy Constructor
Polynomial::Polynomial(const Polynomial& alt)
{
	debug("Copying Polynomial");
	m_coeffs = alt.m_coeffs;
}

//Move Constructor
Polynomial::Polynomial(Polynomial&& alt) noexcept
{
	debug("Move Constructor");
	m_coeffs = std::move(alt.m_coeffs);
}


//Visualization
std::ostream& operator <<(std::ostream& os, const Polynomial& a)
{
	unsigned int sz = a.m_coeffs.size();
	for (unsigned int i = 0; i < sz; i++)
	{
		os << a.m_coeffs[i];
		if (i == sz - 2)
			os << "x +";
		else if (i != sz - 1)
			os << "x^" << sz - i - 1 << " +";
	}
	return os;
}

//Getters
unsigned int Polynomial::Get_Degree() const
{
	return m_coeffs.size() - 1;
}

std::vector<std::complex<double>> Polynomial::Get_Coeffs() const
{
	return m_coeffs;
}


//Operator Overloading
void Polynomial::operator=(const Polynomial& other)
{
	debug("Copying Polynomial");
	m_coeffs = other.m_coeffs;
}



Polynomial Polynomial::operator*(const Polynomial& alt) const
{
	Polynomial n(m_coeffs.size() + alt.m_coeffs.size() - 2);
	unsigned int szr = alt.m_coeffs.size() - 1;
	unsigned int szl = m_coeffs.size() - 1;
	for (unsigned int i = 0; i < m_coeffs.size(); i++)
	{
		for (unsigned int j = 0; j < alt.m_coeffs.size(); j++)
		{
			n[i + j] += (*this).m_coeffs[szl - i] * alt.m_coeffs[szr - j];
		}
	}
	return n;
}

void Polynomial::operator*=(Polynomial& alt)
{
	*this = (*this) * alt;
}



std::complex<double> Polynomial::operator [](unsigned int degree) const
{
	return m_coeffs[Get_Degree() - degree];
}

std::complex<double>& Polynomial::operator [](unsigned int degree)
{
	return m_coeffs[Get_Degree() - degree];
}


Polynomial Polynomial::operator+(const Polynomial& other) const
{
	unsigned int small = 0;
	unsigned int large = 0;
	bool R = false;

	unsigned int r = other.m_coeffs.size();
	unsigned int l = m_coeffs.size();

	if (r > l)
	{
		small = l;
		large = r;
		R = true;
	}
	else
	{
		small = r;
		large = l;
	}

	Polynomial n(large - 1);

	for (unsigned int i = 0; i < large; i++)
	{
		if (i > small - 1)
		{
			if (R)
				n[i] = other.m_coeffs[r - i - 1];
			else
				n[i] = (*this).m_coeffs[l - i - 1];
			continue;
		}
		n[i] = (*this).m_coeffs[l - i - 1] + other.m_coeffs[r - i - 1];
	}

	return n;
}

void Polynomial::operator+=(const Polynomial& other)
{
	*this = std::move(*this + other);
}



Polynomial Polynomial::operator-(const Polynomial& other) const
{
	unsigned int small = 0;
	unsigned int large = 0;
	bool R = false;

	unsigned int r = other.m_coeffs.size();
	unsigned int l = m_coeffs.size();

	if (r > l)
	{
		small = l;
		large = r;
		R = true;
	}
	else
	{
		small = r;
		large = l;
	}

	Polynomial n(large - 1);

	for (unsigned int i = 0; i < large; i++)
	{
		if (i > small - 1)
		{
			if (R)
				n[i] = -other.m_coeffs[r - i - 1];
			else
				n[i] = (*this).m_coeffs[l - i - 1];
			continue;
		}
		n[i] = (*this).m_coeffs[l - i - 1] - other.m_coeffs[r - i - 1];
	}

	return n;
}

void Polynomial::operator-=(const Polynomial& other)
{
	*this = std::move(*this - other);
}



bool Polynomial::operator==(const Polynomial& other) const
{
	return m_coeffs == other.m_coeffs;
}

bool Polynomial::operator!=(const Polynomial& other) const
{
	return !(*this == other);
}