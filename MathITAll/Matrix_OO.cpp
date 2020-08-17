#include"Matrix.h"
//Constructor
Matrix::Matrix(const unsigned int& _m, const unsigned int& _n) :
	m_nrow(_m), m_ncol(_n), m_mat(std::vector<std::complex<double>>(_m* _n, 0))
{
	debug("Construction of Matrix--(0)");
}

Matrix Matrix::Identity(const unsigned int& _n)
{
	Matrix ret(_n, _n);

	for (unsigned int i = 0; i < _n; i++)
	{
		ret.m_mat[i + i * _n] = 1;
	}
	return ret;
}

Matrix Matrix::AxisVector(const unsigned int& _dim, const unsigned int& _m)
{
	Matrix ret(_dim, 1);
	ret.Get(_m, 0) = 1;
	return ret;
}

Matrix Matrix::AxisVector(const unsigned int& _m) const
{
	Matrix ret(m_nrow, 1);
	ret.Get(_m, 0) = 1;
	return ret;
}

Matrix::Matrix(std::vector<Matrix> e) :
	m_nrow(e[0].m_nrow), m_ncol(e.size()), m_mat(std::vector<std::complex<double>>(m_nrow* m_ncol))
{
	for (unsigned int i = 0; i < e.size(); i++)
	{
		SetCol(e[i], i);
	}

}


//Destructor
Matrix::~Matrix()
{
	debug("Destruction of Matrix");
}

//Copy Constructor
Matrix::Matrix(const Matrix& alt) :
	m_nrow(alt.m_nrow), m_ncol(alt.m_ncol), m_mat(alt.m_mat)
{
	debug("Matrix Copied");
}

//Move Constructor
Matrix::Matrix(Matrix&& alt) noexcept :
	m_nrow(alt.m_nrow), m_ncol(alt.m_ncol)
{
	m_mat = std::move(alt.m_mat);
	debug("Matrix Moved");
}

//Visualization
std::ostream& operator <<(std::ostream& os, const Matrix& alt)
{
	os << "---------------\n";
	for (unsigned int i = 0; i < alt.m_nrow; i++)
	{
		for (unsigned int j = 0; j < alt.m_ncol; j++)
		{
			os << alt.m_mat[j + i * alt.m_ncol];
			os << " | ";
		}
		os << "\n";
	}
	os << "---------------\n";
	return os;
}

//Operator Overloading
void Matrix::operator=(const Matrix& alt)
{
	debug("Matrix Copied");
	m_ncol = alt.m_ncol;
	m_nrow = alt.m_nrow;
	m_mat = alt.m_mat;
}

void Matrix::operator=(Matrix&& alt) noexcept
{
	debug("Matrix Moved");
	m_ncol = std::move(alt.m_ncol);
	m_nrow = std::move(alt.m_nrow);
	m_mat = std::move(alt.m_mat);
}



Matrix Matrix::operator*(const Matrix& alt) const
{
	if (m_ncol != alt.m_nrow)
		error("Dimesnion Mismatch/Matrix cannot be Multiplied");
	Matrix ret(m_nrow, alt.m_ncol);
	for (unsigned int i = 0; i < m_nrow; i++)
	{
		for (unsigned int j = 0; j < alt.m_ncol; j++)
		{
			for (unsigned int k = 0; k < m_ncol; k++)
			{
				ret.Get(i, j) += Get(i, k) * alt.Get(k, j);
			}
		}
	}
	return ret;
}

void Matrix::operator*=(const Matrix& alt)
{
	*this = std::move(*this * alt);
}


Matrix Matrix::operator+(const Matrix& alt) const
{
	if (!isSameDim(alt))
		error("Dimension mismatch/Matrix cannot be added");

	Matrix ret(m_nrow, m_ncol);
	for (unsigned int i = 0; i < m_ncol * m_nrow; i++)
		ret.m_mat[i] = m_mat[i] + alt.m_mat[i];
	return ret;
}

void Matrix::operator+=(const Matrix& alt)
{
	*this = std::move(*this + alt);
}


Matrix Matrix::operator-(const Matrix& alt) const
{
	if (!isSameDim(alt))
		error("Dimension mismatch/Matrix cannot be added");

	Matrix ret(m_nrow, m_ncol);
	for (unsigned int i = 0; i < m_ncol * m_nrow; i++)
		ret.m_mat[i] = m_mat[i] - alt.m_mat[i];
	return ret;
}

void Matrix::operator-=(const Matrix& alt)
{
	*this = std::move(*this - alt);
}

Matrix Matrix::operator-() const
{
	return *this * -1;
}


bool Matrix::operator==(const Matrix& alt) const
{
	if (!isSameDim(alt))
		return false;
	for (unsigned int i = 0; i < m_ncol * m_nrow; i++)
		if (m_mat[i] != alt.m_mat[i])
			return false;
	return true;
}

bool Matrix::operator!=(const Matrix& alt) const
{
	return !(*this == alt);
}

std::complex<double>& Matrix::operator()(const unsigned int& _m, const unsigned int& _n)
{
	if (_m >= m_nrow || _n >= m_ncol)
		error("Index out of range");

	return m_mat[m_ncol * _m + _n];
}

std::complex<double> Matrix::operator()(const unsigned int& _m, const unsigned int& _n) const
{
	if (_m >= m_nrow || _n >= m_ncol)
		error("Index out of range");

	return m_mat[m_ncol * _m + _n];
}



//Getters & setters
unsigned int Matrix::Get_size_rows() const
{
	return m_nrow;
}

unsigned int Matrix::Get_size_cols() const
{
	return m_ncol;
}



std::complex<double>& Matrix::Get(const unsigned int& _m, const unsigned int& _n)
{
	if (_m >= m_nrow || _n >= m_ncol)
		error("Index out of range");

	return m_mat[m_ncol * _m + _n];
}

std::complex<double> Matrix::Get(const unsigned int& _m, const unsigned int& _n) const
{
	if (_m >= m_nrow || _n >= m_ncol)
		error("Index out of range");

	return m_mat[m_ncol * _m + _n];
}



Matrix Matrix::GetRow(const unsigned int& _m) const
{
	Matrix ret(1, m_ncol);

	for (unsigned int i = 0; i < m_ncol; i++)
		ret.Get(0, i) = Get(_m, i);

	return ret;
}

Matrix Matrix::GetRow(const Matrix& alt, const unsigned int& _m)
{
	return alt.GetRow(_m);
}



Matrix Matrix::GetCol(const unsigned int& _n) const
{
	Matrix ret(m_nrow, 1);

	for (unsigned int i = 0; i < m_nrow; i++)
		ret.Get(i, 0) = Get(i, _n);

	return ret;
}

Matrix Matrix::GetCol(const Matrix& alt, const unsigned int& _n)
{
	return alt.GetCol(_n);
}


Matrix Matrix::GetPart(unsigned int m1, unsigned int m2, unsigned int n1, unsigned int n2) const
{
	if (m1 == NULL)
		m1 = 0;
	if (m2 == NULL)
		m2 = m_nrow;
	if (n1 == NULL)
		n1 = 0;
	if (n2 == NULL)
		n2 = m_ncol;
	if (m1 >= m2 || n1 >= n2)
		error("Matrix Cannot have one dimension zero or negative");


	Matrix ret(m2 - m1, n2 - n1);

	for (unsigned int i = m1; i < m2; i++)
		for (unsigned int j = n1; j < n2; j++)
			ret.Get(i - m1, j - n1) = Get(i, j);

	return ret;
}

void Matrix::SetPart(const Matrix& alt, unsigned int m1, unsigned int m2, unsigned int n1, unsigned int n2)
{
	if (m1 == NULL)
		m1 = 0;
	if (m2 == NULL)
		m2 = m_nrow;
	if (n1 == NULL)
		n1 = 0;
	if (n2 == NULL)
		n2 = m_ncol;

	if (m1 >= m2 || n1 >= n2)
		error("SetPart: Matrix Cannot have one dimension zero or negative");
	if (m2 - m1 != alt.m_nrow || n2 - n1 != alt.m_ncol)
		error("SetPart: Matrix Size Mismatch");

	for (unsigned int i = m1; i < m2; i++)
		for (unsigned int j = n1; j < n2; j++)
			Get(i, j) = alt.Get(i - m1, j - n1);

}


void Matrix::SetRow(const Matrix& alt, const unsigned int& _m)
{
	if (alt.m_ncol != m_ncol)
		error("Size Missmatch");
	if (_m + alt.m_nrow > m_nrow)
		error("Matrix Size Too Big Or Index out of range");

	for (unsigned int i = 0; i < alt.m_nrow; i++)
		for (unsigned int j = 0; j < m_ncol; j++)
			Get(i + _m, j) = alt.Get(i, j);
}

Matrix Matrix::SetRow(Matrix A, const Matrix& alt, const unsigned int& _m)
{
	A.SetRow(alt, _m);
	return A;
}



void Matrix::SetCol(const Matrix& alt, const unsigned int& _n)
{
	if (alt.m_nrow != m_nrow)
		error("Size Missmatch");
	if (_n + alt.m_ncol > m_ncol)
		error("Matrix Size Too Big Or Index out of range");

	for (unsigned int i = 0; i < m_nrow; i++)
		for (unsigned int j = 0; j < alt.m_ncol; j++)
			Get(i, j + _n) = alt.Get(i, j);
}

Matrix Matrix::SetCol(Matrix A, const Matrix& alt, const unsigned int& _n)
{
	A.SetCol(alt, _n);
	return A;
}



void Matrix::DeleteRow(const unsigned int& _m)
{
	if (_m >= m_nrow)
		error("Index Out of range");
	unsigned int change = 0;

	for (unsigned int i = 0; i < m_nrow * m_ncol; i++)
	{
		if (std::floor(i / (m_ncol)) != _m)
			m_mat[i - change] = m_mat[i];
		else
			change += 1;
	}
	m_mat.erase(m_mat.end() - m_ncol, m_mat.end());
	m_nrow -= 1;
}

Matrix Matrix::DeleteRow(const Matrix& alt, const unsigned int& _m)
{
	if (_m >= alt.m_nrow)
		error("Index Out of range");

	Matrix ret(alt.m_nrow - 1, alt.m_ncol);
	unsigned int change = 0;

	for (unsigned int i = 0; i < alt.m_nrow * alt.m_ncol; i++)
	{
		if (std::floor(i / (alt.m_ncol)) != _m)
			ret.m_mat[i - change] = alt.m_mat[i];
		else
			change += 1;
	}
	return ret;
}



void Matrix::DeleteCol(const unsigned int& _n)
{
	if (_n >= m_ncol)
		error("Index Out of range");

	unsigned int change = 0;
	for (unsigned int i = 0; i < m_nrow * m_ncol; i++)
	{
		if (i % m_ncol != _n)
			m_mat[i - change] = m_mat[i];
		else
			change += 1;
	}
	m_mat.erase(m_mat.end() - m_nrow, m_mat.end());
	m_ncol -= 1;
}

Matrix Matrix::DeleteCol(const Matrix& alt, const unsigned int& _n)
{
	if (_n >= alt.m_ncol)
		error("Index Out of range");

	Matrix ret(alt.m_nrow, alt.m_ncol - 1);
	unsigned int change = 0;
	for (unsigned int i = 0; i < alt.m_nrow * alt.m_ncol; i++)
	{
		if (i % alt.m_ncol != _n)
			ret.m_mat[i - change] = alt.m_mat[i];
		else
			change += 1;
	}
	return ret;
}



Matrix Matrix::InsertRow(const Matrix& A, const Matrix& alt, unsigned int _m)
{
	if (alt.m_ncol != A.m_ncol)
		error("Column Size Mismatch");

	if (_m == NULL)
		_m = A.m_nrow;

	Matrix ret(A.m_nrow + alt.m_nrow, A.m_ncol);

	unsigned int inserted = 0;
	for (unsigned int i = 0; i < A.m_nrow + 1; i++)
	{
		if (i == _m && inserted < alt.m_nrow)
		{
			for (unsigned int j = 0; j < A.m_ncol; j++)
				ret.Get(i + inserted, j) = alt.Get(inserted, j);
			inserted += 1;
			i -= 1;
		}
		else if (i < A.m_nrow)
			for (unsigned int j = 0; j < A.m_ncol; j++)
				ret.Get(i + inserted, j) = A.Get(i, j);
	}
	return ret;
}

void Matrix::InsertRow(const Matrix& alt, unsigned int _m)
{
	*this = std::move(InsertRow(*this, alt, _m));
}


Matrix Matrix::InsertCol(const Matrix& A, const Matrix& alt, unsigned int _n)
{
	if (alt.m_nrow != A.m_nrow)
		error("Row Size Mismatch");

	if (_n == NULL)
		_n = A.m_ncol;

	Matrix ret(A.m_nrow, A.m_ncol + alt.m_ncol);

	unsigned int inserted = 0;

	for (unsigned int i = 0; i < A.m_ncol + 1; i++)
	{
		if (i == _n && inserted < alt.m_ncol)
		{
			for (unsigned int j = 0; j < A.m_nrow; j++)
				ret.Get(j, i + inserted) = alt.Get(j, inserted);

			i -= 1;
			inserted += 1;
		}
		else if (i < A.m_ncol)
			for (unsigned int j = 0; j < A.m_nrow; j++)
				ret.Get(j, i + inserted) = A.Get(j, i);
	}

	return ret;
}

void Matrix::InsertCol(const Matrix& alt, unsigned int _n)
{
	*this = std::move(InsertCol(*this, alt, _n));
}




//Checks
bool Matrix::isSquare() const
{
	if (m_ncol == m_nrow)
		return true;
	return false;
}

bool Matrix::isSquare(const Matrix& alt)
{
	return alt.isSquare();
}


bool Matrix::isSymm() const
{
	if (!isSquare())
		return false;

	for (unsigned int i = 0; i < m_ncol; i++)
		for (unsigned int j = 0; j < m_nrow; j++)
			if (Get(i, j) != Get(j, i))
				return false;

	return true;
}

bool Matrix::isSymm(const Matrix& alt)
{
	return alt.isSymm();
}


bool Matrix::isSkewSymm() const
{
	if (!isSquare())
		return false;

	for (unsigned int i = 0; i < m_ncol; i++)
		for (unsigned int j = 0; j < m_nrow; j++)
			if (Get(i, j) != -Get(j, i))
				return false;

	return true;
}

bool Matrix::isSkewSymm(const Matrix& alt)
{
	return alt.isSkewSymm();
}


bool Matrix::isSameDim(const Matrix& B) const
{
	if (m_ncol == m_ncol && m_nrow == m_nrow)
		return true;
	return false;
}

bool Matrix::isSameDim(const Matrix& A, const Matrix& B)
{
	return A.isSameDim(B);
}