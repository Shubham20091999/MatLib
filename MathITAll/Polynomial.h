#ifndef Poly_h
#define Poly_h

#include"MathIT.h"

class Polynomial
{
private:
	std::vector<std::complex<double>> m_coeffs;
	
	template<typename T>
	void DivideBySolution(const T& c);
	std::complex<double> FindSol_N(const double& StopCriterion) const;
	//std::vector<std::complex<double>> FindSol_A() const;

public:
	//Constructors;
	Polynomial()=delete;

	Polynomial(unsigned int deg);

	template<typename T>
	Polynomial(unsigned int deg, T* ls);

	template<typename T>
	Polynomial(std::initializer_list<T> args);

	template<typename T>
	Polynomial(std::initializer_list<std::initializer_list<T>> args);


	//Destructor
	~Polynomial();

	//Copy Constructor
	Polynomial(const Polynomial& alt);

	//Move Constructor
	Polynomial(Polynomial&& alt) noexcept;

	//Visualization
	friend std::ostream& operator <<(std::ostream& os, const Polynomial& a);

	//Getters
	unsigned int Get_Degree() const;
	std::vector<std::complex<double>> Get_Coeffs() const;

	//Operator Overloaded
	template<typename T>
	std::complex<double> operator()(const T& x) const;

	void operator=(const Polynomial& other);
	
	template<typename T>
	void operator*=(const T& c);
	template<typename T>
	Polynomial operator*(const T& other) const;
	
	template<typename T>
	friend Polynomial operator*(const T& other, const Polynomial& A);

	Polynomial operator*(const Polynomial& alt) const;
	void operator*=(Polynomial& alt);

	template<typename T>
	void operator/=(const T& c);
	template<typename T>
	Polynomial operator/(const T& c) const;


	Polynomial operator+(const Polynomial& other) const;
	void operator += (const Polynomial& other);

	Polynomial operator-(const Polynomial& other) const;
	void operator-=(const Polynomial& other);

	std::complex<double> operator [](unsigned int degree) const;
	std::complex<double>& operator [](unsigned int degree);

	bool operator==(const Polynomial& other) const;
	bool operator!=(const Polynomial& other) const;

	//Operations
	void Differentiate(unsigned int n = 1);
	void Integrate(double x = 0, double y = 0);
	std::vector<std::complex<double>> FindAllSol_N(const double& StopCriterion = 10e-15) const;
	std::vector<std::complex<double>> FindAllSol_A(bool safe = false, const double& StopCriterion = 10e-10) const;

	//Static-Operation
	static Polynomial Differentiate(const Polynomial& alt, unsigned int n);
	static Polynomial Integrate(const Polynomial& alt, double x = 0, double y = 0);
	static std::vector<std::complex<double>> FindAllSol_N(const Polynomial& A, const double& StopCriterion = 10e-15);
	static std::vector<std::complex<double>> FindAllSol_A(const Polynomial& A,bool safe=false, const double& StopCriterion = 10e-10);
};

#include"Polynomial_t.h"


#endif//Poly_h