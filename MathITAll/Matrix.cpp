#include"Matrix.h"


//Operations
Matrix Matrix::Transpose() const
{
	Matrix ret(m_ncol, m_nrow);
	for (unsigned int i = 0; i < m_nrow; i++)
		for (unsigned int j = 0; j < m_ncol; j++)
			ret.Get(j, i) = Get(i, j);
	return ret;
}

Matrix Matrix::Transpose(const Matrix& A)
{
	return A.Transpose();
}



Matrix Matrix::Hermitian(const Matrix& A)
{
	return A.Hermitian();
}

Matrix Matrix::Hermitian() const
{
	Matrix ret(m_ncol, m_nrow);
	for (unsigned int i = 0; i < m_nrow; i++)
		for (unsigned int j = 0; j < m_ncol; j++)
			ret.Get(j, i) = std::conj(Get(i, j));
	return ret;

}



Matrix Matrix::Anti_Transpose() const
{
	return Anti_Transpose(*this);
}

Matrix Matrix::Anti_Transpose(const Matrix& A)
{
	Matrix ret(A.m_ncol, A.m_nrow);
	for (unsigned int i = 0; i < A.m_nrow; i++)
	{
		for (unsigned int j = 0; j < A.m_ncol; j++)
			ret.Get(j, i) = A.Get(A.m_nrow - i - 1, A.m_ncol - j - 1);
	}
	return ret;
}


void Matrix::InterchangeRow(unsigned int p, unsigned int q)
{
	for (unsigned int i = 0; i < m_ncol; i++)
	{
		std::complex<double> temp = Get(p, i);
		Get(p, i) = Get(q, i);
		Get(q, i) = temp;
	}
}


Matrix Matrix::UTform(Matrix A, Matrix* B)
{
	A.UTform(B);
	return A;
}

void Matrix::UTform(Matrix* B)
{
	unsigned int smll = (m_nrow < m_ncol) ? m_nrow : m_ncol;

	for (unsigned int i = 0; i < smll; i++)
	{
		double max = abs(Get(i, i));
		unsigned int pos = i;
		for (unsigned int j = i + 1; j < m_nrow; j++)
			if (max < abs(Get(j, i)))
			{
				max = abs(Get(j, i));
				pos = j;
			}

		if (i != pos)
		{
			InterchangeRow(i, pos);
			if (B)
				B->InterchangeRow(i, pos);
		}

		if (isZero(Get(i, i)))
			for (unsigned int j = i + 1; j < m_nrow; j++)
				if (!isZero(Get(j, i)))
				{
					InterchangeRow(i, j);
					if (B)
						B->InterchangeRow(i, j);
				}

		if (!isZero(Get(i, i)))
			for (unsigned int j = i + 1; j < m_nrow; j++)
				if (!isZero(Get(j, i)))
				{
					std::complex<double> change = Get(j, i) / Get(i, i);
					for (unsigned int k = i; k < m_ncol; k++)
						Get(j, k) -= change * Get(i, k);
					if (B)
						for (unsigned int k = 0; k < B->m_ncol; k++)
							B->Get(j, k) -= change * B->Get(i, k);
				}
	}
}


lu Matrix::LU(const Matrix& input, unsigned int* cnt)
{
	Matrix U = input;
	Matrix P = Matrix::Identity(U.m_nrow);
	Matrix L = Identity(U.m_nrow);
	unsigned int count = 0;

	unsigned int smll = (U.m_nrow < U.m_ncol) ? U.m_nrow : U.m_ncol;

	for (unsigned int i = 0; i < smll; i++)
	{
		double max = abs(U.Get(i, i));
		unsigned int pos = i;

		for (unsigned int j = i + 1; j < U.m_nrow; j++)
		{
			if (max < abs(U.Get(j, i)))
			{
				max = abs(U.Get(j, i));
				pos = j;
			}
		}
		if (i != pos)
		{
			U.InterchangeRow(i, pos);
			P.InterchangeRow(i, pos);
			count += 1;
		}

		if (isZero(U.Get(i, i)))
			for (unsigned int j = i + 1; j < U.m_nrow; j++)
				if (!isZero(U.Get(j, i)))
				{
					U.InterchangeRow(i, j);
					P.InterchangeRow(i, j);
					count += 1;
				}

		if (!isZero(U.Get(i, i)))
			for (unsigned int j = i + 1; j < U.m_nrow; j++)
				if (!isZero(U.Get(j, i)))
				{
					std::complex<double> change = U.Get(j, i) / U.Get(i, i);
					for (unsigned int k = i; k < U.m_ncol; k++)
						U.Get(j, k) -= change * U.Get(i, k);
				}
	}

	Matrix tmp = P * input;

	for (unsigned int i = 0; i < L.m_ncol - 1; i++)
	{
		for (unsigned int j = i + 1; j < L.m_nrow; j++)
		{
			L.Get(j, i) = tmp.Get(j, i);
			for (unsigned int k = 0; k < i; k++)
			{
				L.Get(j, i) -= L.Get(j, k) * U.Get(k, i);
			}
			L.Get(j, i) /= U.Get(i, i);
		}
	}
	if (cnt)
		*cnt = count;
	return lu{ std::move(P), std::move(L), std::move(U) };
}



std::complex<double> Matrix::Determinent(const Matrix& A)
{
	if (!A.isSquare())
		error("Determinent of Non Square Matrix Does not Exist ");

	unsigned int count = 0;
	auto tmp = LU(A, &count);
	std::complex<double> ret = (count % 2 == 0) ? 1 : -1;
	for (unsigned int i = 0; i < A.m_ncol; i++)
		ret *= tmp.L.Get(i, i) * tmp.U.Get(i, i);

	return ret;
}

std::complex<double> Matrix::Determinent() const
{
	return Determinent(*this);
}

std::complex<double> Matrix::Trace(const Matrix& A)
{
	if (!A.isSquare())
		error("Trace of rectangular Matrix Doesnot exist");

	std::complex<double> ret = 0;
	for (unsigned int i = 0; i < A.m_ncol; i++)
		ret += A.Get(i, i);

	return ret;
}

std::complex<double> Matrix::Trace() const
{
	return Trace(*this);
}

Polynomial Matrix::CharPoly(const Matrix& A)
{
	if (!isSquare(A))
		error("Characteristic polynomial of non square matrix doesnot exixt");

	Matrix tmp = A;
	Matrix I = Identity(A.m_ncol);

	Polynomial ret(A.m_ncol);

	ret[A.m_ncol] = 1;

	ret[A.m_ncol - 1] = -tmp.Trace();

	for (unsigned int i = 2; i < A.m_ncol + 1; i++)
	{
		tmp = A * (tmp + I * ret[A.m_ncol - i + 1]);
		ret[A.m_ncol - i] = -tmp.Trace() / std::complex<double>(i);
	}

	return ret;
}


std::vector<std::complex<double>> Matrix::EigenValues(const Matrix& A,unsigned int* zeros)
{
	auto lams= CharPoly(A).FindAllSol_A();
	SameTogether(lams);
	if(zeros)
		*zeros=MoveZeroAtEnd(lams);

	return lams;
}



void Matrix::EF(Matrix* B)
{
	this->UTform(B);
	unsigned int smll = (m_ncol < m_nrow) ? m_ncol : m_nrow;

	for (unsigned int i = 0; i < smll; i++)
	{
		if (isZero(Get(i, i)))
			for (unsigned int j = i + 1; j < m_ncol; j++)
			{
				bool done = false;
				if (isZero(Get(i, j)))
					for (unsigned int k = i + 1; k < m_nrow; k++)
					{
						if (!isZero(Get(k, j)))
						{
							InterchangeRow(k, i);
							if (B)
								B->InterchangeRow(k, i);
							done = true;
							break;
						}
					}
				if (done)
					break;
			}
	}

	for (unsigned int i = 0; i < m_nrow; i++)
	{
		bool found = false;
		std::complex<double> change = 1.0;
		for (unsigned int j = i; j < m_ncol; j++)
		{
			if (!isZero(Get(i, j)) && !found)
			{
				found = true;
				change = Get(i, j);
			}
			Get(i, j) /= change;
		}
		if (B && change != 1.0)
			for (unsigned int j = 0; j < B->m_ncol; j++)
				B->Get(i, j) /= change;
	}

}

Matrix Matrix::EF(Matrix A, Matrix* B)
{
	A.EF(B);
	return A;
}


void Matrix::RREF(Matrix* B)
{
	EF(B);

	unsigned int smll = (m_ncol < m_nrow) ? m_ncol : m_nrow;

	for (unsigned int i = 0; i < smll; i++)
	{
		for (unsigned int j = i; j < m_ncol; j++)
			if (!isZero(Get(i, j)))
			{
				for (unsigned int k = 0; k < m_nrow; k++)
					if (!isZero(Get(k, j)) && k != i)
					{
						std::complex<double> change = Get(k, j) / Get(i, j);
						for (unsigned int l = j; l < m_ncol; l++)
							Get(k, l) -= change * Get(i, l);
						if (B)
							for (unsigned int l = 0; l < B->m_ncol; l++)
								B->Get(k, l) -= change * B->Get(i, l);
					}
				break;
			}
	}
}

Matrix Matrix::RREF(Matrix A, Matrix* B)
{
	A.RREF(B);
	return A;
}


std::vector<Matrix> Matrix::NullSpace(Matrix A)
{
	std::vector<Matrix> ret;
	A.RREF();
	unsigned int smll = (A.m_ncol < A.m_nrow) ? A.m_ncol : A.m_nrow;
	for (unsigned int i = 0; i < smll; i++)
		for (unsigned int j = i; j < A.m_ncol; j++)
			if (!isZero(A.Get(i, j)))
			{
				if(i==j)
					break;
				A.InterchangeRow(i, j);
				i--;
				break;
			}
	for (unsigned int i = 0; i < smll; i++)
	{
		bool nopivot = true;

		for (unsigned int j = i; j < A.m_ncol; j++)
			if (!isZero(A.Get(i, j)))
			{
				nopivot = false;
				break;
			}

		if (nopivot)
		{
			Matrix tmp(A.m_ncol, 1);
			for (unsigned int j = 0; j < A.m_nrow; j++)
				tmp.Get(j, 0) = A.Get(j, i);

			tmp.Get(i, 0) = -1;
			ret.emplace_back(tmp);
		}
	}
	if (A.m_nrow < A.m_ncol)
		for (unsigned int i = 0; i < A.m_ncol - A.m_nrow; i++)
		{
			Matrix tmp(A.m_ncol, 1);
			for (unsigned int j = 0; j < A.m_nrow; j++)
				tmp.Get(j, 0) = A.Get(j, A.m_nrow + i);

			tmp.Get(A.m_nrow + i, 0) = -1;
			ret.emplace_back(tmp);
		}

	if (ret.size() == 0)
		ret.emplace_back(Matrix(A.m_ncol, 1));

	return ret;
}

std::vector<Matrix> Matrix::NullSpace()
{
	return NullSpace(*this);
}



std::complex<double> Matrix::Cofactor(unsigned int _m, unsigned int _n)
{
	Matrix tmp = DeleteCol(DeleteRow(*this,_m),_n);
	return std::pow(-1, _m + _n) * Matrix::Determinent(tmp);
}


std::vector<Matrix> Matrix::EigenVector(const Matrix& A, const std::complex<double>& lam)
{	
	return NullSpace(A - Identity(A.m_ncol) * lam);
}

std::vector<Matrix> Matrix::EigenVector(const std::complex<double>& lam) const
{
	return EigenVector(*this,lam);
}


std::vector<Matrix> Matrix::AllEigenVectors() const
{
	std::vector<std::complex<double>> lams = EigenValues(*this);

	RemoveDuplicates_Complex(lams);
	std::vector<Matrix> ret;
	for (auto sol : lams)
	{
			auto tmp = EigenVector(sol);
			ret.insert(ret.end(), tmp.begin(), tmp.end());
	}

	return ret;
}

std::vector<Matrix> Matrix::AllEigenVectors(const Matrix& A)
{
	return A.AllEigenVectors();
}


std::complex<double> Matrix::Dot(const Matrix& B) const
{
	return Trace(Transpose(*this) * Conjugate(B));
}

std::complex<double> Matrix::Dot(const Matrix& A, const Matrix& B)
{
	return A.Dot(B);
}



void Matrix::Orthonormalize()
{
	*this=Orthonormalize(*this);
}

Matrix Matrix::Orthonormalize(const Matrix& A)
{
	Matrix QT = Identity(A.m_nrow);
	for (unsigned int i = 0; i < A.m_ncol; i++)
	{
		auto y = A.GetPart(i, NULL, i, i + 1);
		auto w = y + y.AxisVector(0) * y(0, 0) / abs(y(0, 0)) * Magnitude(y);
		auto v = w / Magnitude(w);

		Matrix H = Matrix::Identity(A.m_nrow);
		H.SetPart(Matrix::Identity(v.m_nrow) - v * Hermitian(v) * 2, i, NULL, i, NULL);

		QT = H * QT;
	}
	return Hermitian(QT);
}

void Matrix::Orthonormalize(std::vector<Matrix>& A)
{
	Matrix rett(A[0].m_nrow,A.size());
	for (unsigned int i = 0; i < A.size(); i++)
	{
		 rett.SetCol(A[i],i);
	}

	rett.Orthonormalize();

	for (unsigned int i = 0; i < A.size(); i++)
	{
		A[i]=rett.GetCol(i);
	}
}


qr Matrix::QR() const
{
	Matrix R=*this;
	Matrix QT=Identity(m_nrow);
	for (unsigned int i = 0; i < m_ncol; i++)
	{
		auto y=R.GetPart(i,NULL,i,i+1);
		auto w=y+y.AxisVector(0)*y(0,0)/abs(y(0,0))*Magnitude(y);
		auto v=w/Magnitude(w);
		
		Matrix H=Matrix::Identity(m_nrow);
		H.SetPart(Matrix::Identity(v.m_nrow) - v * Hermitian(v) * 2,i,NULL,i,NULL);

		QT=H*QT;
		R=H*R;
	}
	return qr{Hermitian(QT),R};
}

qr Matrix::QR(const Matrix& A)
{
	return A.QR();
}


std::complex<double> Matrix::Magnitude() const
{
	return Magnitude(*this);
}

std::complex<double> Matrix::Magnitude(const Matrix& A)
{
	if (A.m_ncol != 1)
		error("Matrix is Not a Column Vector");

	std::complex<double> ret = 0;
	for (unsigned int i = 0; i < A.m_nrow; i++)
		ret += norm(A.Get(i, 0));
	return sqrt(ret);
}


Matrix Matrix::Cholesky() const
{
	return Cholesky(*this);
}

Matrix Matrix::Cholesky(const Matrix& A)
{
	if (!isSymm(A))
		error("Cholesky Decomposition of Non Symm. matrix does not exist");
	
	for (unsigned int i = 0; i < A.m_ncol; i++)
	{
		auto det=Determinent(A.GetPart(NULL, i, NULL, i));

		if(isComplex(det) || (!isZero(det) && det.real()<0))
			error("Cholesky Decomposition only on Semi/- Positive Definite Matrix");
	}

	Matrix L(A.m_nrow, A.m_ncol);

	for (unsigned int i = 0; i < A.m_ncol; i++)
	{
		L.Get(i, i) = A.Get(i, i);
		for (unsigned int j = 0; j < i; j++)
			L.Get(i, i) -= L(i, j) * std::conj(L(i, j));
		L.Get(i, i) = sqrt(L.Get(i, i));

		for (unsigned int j = i + 1; j < A.m_nrow; j++)
		{
			L.Get(j, i) = A.Get(j, i);
			for (unsigned int k = 0; k < i; k++)
			{
				L.Get(j, i) -= L.Get(j, k) * conj(L.Get(i, k));
			}
			L.Get(j, i) /= L(i, i);
		}
	}
	return L;
}



Matrix Matrix::Invert_U(Matrix A)
{
	if(!isSquare(A))
		error("Rectangular Matrix Inverse Doesnot Exist");

	Matrix Inv=Identity(A.m_ncol);

	for (unsigned int i = 0; i < A.m_ncol; i++)
	{
		std::complex<double> change = A.Get(i, i);
		if(isZero(A.Get(i, i)))
			error("Singular Matrix cannot be inverted");
		
		for(unsigned int j=i;j<A.m_ncol;j++)
			A.Get(i,j)/=change;
		for (unsigned int j = 0; j < Inv.m_ncol; j++)
			Inv.Get(i, j) /= change;
	}

	for (unsigned int i = 1; i < A.m_ncol; i++)
		for (unsigned int j = 0; j < i; j++)
		{
			std::complex<double> change= A.Get(j,i);
			if(!isZero(change))
			{
				for (unsigned int k = i; k < A.m_ncol; k++)
					A.Get(j,k)-=change* A.Get(i,k);
				for (unsigned int k = 0; k < Inv.m_ncol; k++)
					Inv.Get(j, k) -= change * Inv.Get(i, k);
			}
		}
	
	return Inv;
}

Matrix Matrix::Invert_U() const
{
	return Invert_U(*this);
}



Matrix Matrix::Invert_L(Matrix A)
{
	if (!isSquare(A))
		error("Rectangular Matrix Inverse Doesnot Exist");

	Matrix Inv = Identity(A.m_ncol);

	for (unsigned int i = 0; i < A.m_nrow; i++)
	{
		std::complex<double> change = A.Get(i, i);
		if (isZero(A.Get(i, i)))
			error("Singular Matrix cannot be inverted");

		for (unsigned int j = 0; j < i+1; j++)
			A.Get(i, j) /= change;
		for (unsigned int j = 0; j < Inv.m_ncol; j++)
			Inv.Get(i, j) /= change;
	}

	for (unsigned int i = 0; i < A.m_ncol - 1; i++)
	{
		for (unsigned int j = i + 1; j < A.m_nrow; j++)
		{
			std::complex<double> change=A.Get(j,i);

			for (unsigned int k = 0; k < i+1; k++)
				A.Get(j,k)-=change*A.Get(i,k);
			for(unsigned int k=0;k<Inv.m_ncol;k++)
				Inv.Get(j, k) -= change * Inv.Get(i, k);
		}
	}

	return Inv;
}

Matrix Matrix::Invert_L() const
{
	return Invert_L(*this);
}


Matrix Matrix::Invert(const Matrix& A)
{
	auto plu=LU(A);
	return Invert_U(plu.U)*Invert_L(plu.L)*Transpose(plu.P);
}



void Matrix::Normalize()
{
	for (unsigned int i = 0; i < m_ncol; i++)
	{
		auto change=Magnitude(GetCol(i));
		for (unsigned int j = 0; j < m_nrow; j++)
		{	
			Get(j,i)/=change;
		}
	}
}

Matrix Matrix::Normalize(Matrix A)
{
	A.Normalize();
	return A;
}



svd Matrix::SVD(const Matrix& A)
{

	Matrix AAT=A*A.Hermitian();
	unsigned int zeros=0;
	auto eVal=EigenValues(AAT,&zeros);
	Matrix S(AAT.m_nrow,AAT.m_ncol-zeros);

	for (unsigned int i = 0; i < eVal.size()-zeros; i++)
	{
		if(isZero(eVal[i]))
			eVal[i]=0;
		else
			eVal[i]=eVal[i].real();
	}
	for (unsigned int i = 0; i < eVal.size()-zeros; i++)
		S(i,i)=sqrt(eVal[i].real());

	RemoveDuplicates_Complex(eVal);

	Matrix U=Matrix::Identity(AAT.m_nrow);
	unsigned int inserted=0;

	for (unsigned int i = 0; i < eVal.size(); i++)
	{
		auto tmp=std::move(AAT.EigenVector(eVal[i]));
		for (auto s : tmp)
		{
			U.SetCol(s,inserted);
			inserted+=1;
		}
	}
	
	U.GramSchmidth();

	
	Matrix Sinv(S.m_ncol,S.m_nrow);
	for (unsigned int i = 0; i < S.m_ncol; i++)
		Sinv(i,i)=1.0/S(i,i);


	return svd({U,S,Hermitian(Sinv * Hermitian(U) * A)});
}



void Matrix::GramSchmidth()
{
	for (unsigned int i = 0; i < m_ncol; i++)
	{
		Matrix temp = GetCol(i);
		for (unsigned int j = 0; j < i; j++)
		{
			auto tmp = GetCol(j);
			temp -= tmp * Dot(GetCol(i), tmp);
		}
		SetCol(temp / Magnitude(temp), i);
	}
}

Matrix Matrix::GramSchmidth(Matrix A)
{
	A.GramSchmidth();
	return A;
}


void Matrix::Conjugate()
{
	*this = Conjugate(*this);
}

Matrix Matrix::Conjugate(const Matrix& A)
{
	Matrix ret(A.m_nrow, A.m_ncol);
	for (unsigned int i = 0; i < A.m_ncol * A.m_nrow; i++)
		ret.m_mat[i] = conj(A.m_mat[i]);
	return ret;
}