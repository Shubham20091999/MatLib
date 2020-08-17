#ifndef Mat_t
#define Mat_t
#include"Matrix.h"

template<typename T>
Matrix::Matrix(std::initializer_list<std::initializer_list<T>> e) :
	m_nrow(e.size()), m_ncol((*e.begin()).size()), m_mat(std::vector<std::complex<double>>(m_nrow* m_ncol))
{
	debug("Construction of Matrix--(1)");
	unsigned int i = 0;
	for (auto row : e)
	{
		for (auto col : row)
		{
			m_mat[i] =std::complex<double>(col);
			i += 1;
		}
		if (i % m_ncol != 0)
			error("Mismatch Column Size");
	}
}

template<typename T>
Matrix::Matrix(std::initializer_list<std::initializer_list<std::initializer_list<T>>> e):
	m_nrow(e.size()), m_ncol((*e.begin()).size()), m_mat(std::vector<std::complex<double>>(m_nrow* m_ncol))
{
	debug("Construction of Matrix--(2)");
	unsigned int i = 0;
	for (auto row : e)
	{
		for (auto col : row)
		{
			if(col.size()>2)
				error("Cannot Be Converted to Complex Number");
			else if(col.size()==2)
				m_mat[i] = std::complex<double>(*col.begin(),*(col.end()-1));
			else if(col.size()==1)
				m_mat[i]=std::complex<double>(*col.begin());
			else
				m_mat[i]=std::complex<double>(0);
			i += 1;
		}
		if (i % m_ncol != 0)
			error("Mismatch Column Size");
	}
}

template<typename T>
void Matrix::operator/=(const T& c)
{
	*this = std::move(*this / c);
}

template<typename T>
Matrix Matrix::operator/ (const T& c) const
{
	Matrix ret(m_nrow, m_ncol);
	for (unsigned int i = 0; i < m_ncol * m_nrow; i++)
		ret.m_mat[i] = m_mat[i] / std::complex<double>(c);
	return ret;
}


template<typename T>
Matrix Matrix::operator*(const T& c) const
{
	Matrix ret(m_nrow, m_ncol);
	for (unsigned int i = 0; i < m_ncol * m_nrow; i++)
		ret.m_mat[i] = std::complex<double>(c) * m_mat[i];
	return ret;
}

template<typename T>
void Matrix::operator*=(const T& c)
{
	*this=std::move(*this*c);
}

template<typename T>
Matrix operator*(const T& c, const Matrix& alt)
{
	return alt*c;
}
	
#endif // !Mat_t
