#ifndef Mat_h
#define Mat_h

#include"MathIT.h"

class Polynomial;
struct lu;
struct qr;
struct svd;

class Matrix
{
private:
	unsigned int m_nrow;
	unsigned int m_ncol;
	std::vector<std::complex<double>> m_mat;

public:
	//Constructors
	Matrix()=delete;

	Matrix(const unsigned int& _m, const unsigned int& _n);

	template<typename T>
	Matrix(std::initializer_list<std::initializer_list<T>> e);

	template<typename T>
	Matrix(std::initializer_list<std::initializer_list<std::initializer_list<T>>> e);

	explicit Matrix( std::vector<Matrix> e) ;

	static Matrix Identity(const unsigned int& _n);

	static Matrix AxisVector(const unsigned int& _dim, const unsigned int& _m);
	Matrix AxisVector(const unsigned int& _m) const;
	//Destructors
	~Matrix();

	//Copy Constructor
	Matrix(const Matrix& alt);

	//Move Constructor
	Matrix(Matrix&& alt) noexcept;

	//Visualization
	friend std::ostream& operator <<(std::ostream& os, const Matrix& a);

	//Getters & Setters
	unsigned int Get_size_rows() const;

	unsigned int Get_size_cols() const;

	std::complex<double>& Get(const unsigned int& _m, const unsigned int& _n);

	std::complex<double> Get(const unsigned int& _m, const unsigned int& _n) const;

	Matrix GetRow(const unsigned int& _m) const;
	static Matrix GetRow(const Matrix& alt,const unsigned int& _m);

	Matrix GetCol(const unsigned int& _n) const;
	static Matrix GetCol(const Matrix& alt, const unsigned int& _n);

	void SetRow(const Matrix& alt, const unsigned int& _m);
	static Matrix SetRow(Matrix A,const Matrix& alt, const unsigned int& _m);

	void SetCol(const Matrix& alt,const unsigned int& _n);
	static Matrix SetCol(Matrix A,const Matrix& alt, const unsigned int& _n);

	void DeleteRow(const unsigned int& _m);
	static Matrix DeleteRow(const Matrix& alt, const unsigned int& _m);

	void DeleteCol(const unsigned int& _n);
	static Matrix DeleteCol(const Matrix& alt, const unsigned int& _m);

	Matrix GetPart(unsigned int m1 = NULL, unsigned int m2 = NULL, unsigned int n1 = NULL, unsigned int n2 = NULL) const;

	void SetPart(const Matrix& alt, unsigned int m1 = NULL, unsigned int m2 = NULL, unsigned int n1 = NULL, unsigned int n2 = NULL);

	static Matrix InsertRow(const Matrix& A,const Matrix& alt,unsigned int _m=NULL);
	void InsertRow(const Matrix& alt, unsigned int _m = NULL);

	static Matrix InsertCol(const Matrix& A, const Matrix& alt,unsigned int _n=NULL);
	
	void InsertCol(const Matrix& alt, unsigned int _n = NULL);

	//checks
	bool isSquare() const;
	static bool isSquare(const Matrix& alt);

	bool isSymm() const;
	static bool isSymm(const Matrix& alt);

	bool isSkewSymm() const;
	static bool isSkewSymm(const Matrix& alt);

	static bool isSameDim(const Matrix& A, const Matrix& B);
	bool isSameDim(const Matrix& alt) const;



	//Operator Overloading
	//copy
	void operator=(const Matrix& alt);

	//move
	void operator=(Matrix&& alt) noexcept;

	template<typename T>
	Matrix operator/ (const T& c) const;
	template<typename T>
	void operator/=(const T& c);

	template<typename T>
	Matrix operator* (const T& c) const;
	template<typename T>
	void operator*=(const T& c);

	template<typename T>
	friend Matrix operator*(const T& c,const Matrix& alt);

	Matrix operator* (const Matrix& alt) const;
	void operator*=(const Matrix& alt);

	Matrix operator+(const Matrix& alt) const;
	void operator+=(const Matrix& alt);

	Matrix operator-(const Matrix& alt) const;
	void operator-=(const Matrix& alt);
	Matrix operator-() const;

	bool operator==(const Matrix& alt) const;
	bool operator!=(const Matrix& alt) const;

	std::complex<double> operator()(const unsigned int& _m, const unsigned int& _n) const;
	std::complex<double>& operator()(const unsigned int& _m, const unsigned int& _n);


	//Operations
	Matrix Transpose() const;
	static Matrix Transpose(const Matrix& A);

	Matrix Hermitian() const;
	static Matrix Hermitian(const Matrix& A);

	Matrix Anti_Transpose() const;
	static Matrix Anti_Transpose(const Matrix& A);

	void InterchangeRow(unsigned int p, unsigned int q);

	static Matrix UTform(Matrix A, Matrix* B = nullptr);
	void UTform(Matrix* B = nullptr);

	static lu LU(const Matrix& input, unsigned int* GetCount = nullptr);

	static std::complex<double> Determinent(const Matrix& A);
	std::complex<double> Determinent() const;

	std::complex<double> Trace() const;
	static std::complex<double> Trace(const Matrix& alt);

	static Polynomial CharPoly(const Matrix& A);

	static std::vector<std::complex<double>> EigenValues(const Matrix& A, unsigned int* zeros=NULL);

	void EF(Matrix* B = nullptr);
	static Matrix EF(Matrix A, Matrix* B = nullptr);

	void RREF(Matrix* B = nullptr);
	static Matrix RREF(Matrix A,Matrix* B=nullptr);

	static std::vector<Matrix> NullSpace(Matrix A);
	std::vector<Matrix> NullSpace();

	std::complex<double> Cofactor(unsigned int _m, unsigned int _n);

	static std::vector<Matrix> EigenVector(const Matrix& A, const std::complex<double>& lam);
	std::vector<Matrix> EigenVector(const std::complex<double>& lam) const;

	std::vector<Matrix> AllEigenVectors() const;
	static std::vector<Matrix> AllEigenVectors(const Matrix& A);

	static std::complex<double> Dot(const Matrix& A, const Matrix& B);
	std::complex<double> Dot(const Matrix& B) const;

	static std::complex<double> Magnitude(const Matrix& A);
	std::complex<double> Magnitude() const;

	static Matrix Orthonormalize(const Matrix& A);
	void Orthonormalize();

	static void Orthonormalize(std::vector<Matrix>& A);

	qr QR() const;
	static qr QR(const Matrix& A);

	Matrix Cholesky() const;
	static Matrix Cholesky(const Matrix& A);

	static Matrix Invert_U(Matrix U);
	Matrix Invert_U() const;

	static Matrix Invert_L(Matrix L);
	Matrix Invert_L() const;

	static Matrix Invert(const Matrix& A);

	void Normalize();
	static Matrix Normalize(Matrix A);

	static svd SVD(const Matrix& A);

	static Matrix GramSchmidth(Matrix A);
	void GramSchmidth();

	static Matrix Conjugate(const Matrix& A);
	void Conjugate();

};
#include"Matrix_t.h"

struct lu
{
	Matrix P,L,U;
};

struct qr
{
	Matrix Q,R;//Orthonormal,Uppertriangular
};

struct svd
{
	Matrix U,S,V;
};

#endif // !Mat_h