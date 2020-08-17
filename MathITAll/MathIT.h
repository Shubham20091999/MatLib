#ifndef Math
#define Math


#ifdef UseDebug
#define debug(x) std::cout<<x<<"\n"
#else
#define debug(x)
#endif // debug


#ifndef Precision
#define Precision 10e-5
#endif // !Precision


#include<iostream>
#include<vector>
#include<thread>
#include<initializer_list>
#include<complex>
//#include<algorithm>
//#include<cmath>


#define print(x) std::cout<<x<<"\n"
#define printV(x) for(auto a:x) std::cout<<a<<" "; print(" ")

struct Timer
{
	std::chrono::time_point<std::chrono::steady_clock> start, end;

	std::chrono::duration<double> duration = std::chrono::duration<double>(0);

	const char* name = nullptr;
	Timer(const char* _name) :
		start(std::chrono::high_resolution_clock::now()), name(_name) {}

	~Timer()
	{
		end = std::chrono::high_resolution_clock::now();

		duration = end - start;

		double ms = (double)duration.count() * 1000;
		std::cout << name << " : " << ms << " ms" << "\n";
	}
};

struct Counter
{
	unsigned __int32 count;
	const char* name;
	Counter(const char* _name) :count(0), name(_name) {}

	~Counter() { std::cout << "Count: " << name << " : " << count << "\n"; }

	void Increase() { count += 1; }
};


inline int Permutations(unsigned int m, unsigned int n)
{
	unsigned int ans = 1;
	for (unsigned int i = n + 1; i <= m; i++)
	{
		ans *= i;
	}
	return ans;
}

template<typename T>
inline unsigned int RemoveDuplicates_Real(std::vector<T>& v)
{
	unsigned int count = 0;
	unsigned int i = 0;
	while (i < v.size())
	{
		unsigned int j = i + 1;
		while (j < v.size())
		{
			if (abs(v[i]- v[j]) < Precision)
			{
				count += 1;
				v.erase(v.begin() + i);
				i--;
				break;
			}

			j++;
		}
		i++;
	}
	return count;
}

template<typename T>
inline unsigned int RemoveDuplicates_Complex(std::vector<std::complex<T>>& v)
{
	unsigned int count = 0;
	unsigned int i = 0;
	while (i < v.size())
	{
		unsigned int j = i + 1;
		while (j < v.size())
		{
			if (abs(v[i].real() - v[j].real()) < Precision && abs(v[i].imag() - v[j].imag())<Precision)
			{
				count += 1;
				v.erase(v.begin() + i);
				i--;
				break;
			}

			j++;
		}
		i++;
	}
	return count;
}

template<typename T>
inline bool isComplex(std::complex<T>& v)
{
	if (abs(v.imag()) > Precision)
		return true;
	return false;
}

inline void error(const char* a,const std::exception& c=std::exception("Unk"))
{
	print("Error : "<<a);
	throw(c);
}

inline bool isZero(std::complex<double> x)
{
	if(abs(x.real())>Precision || abs(x.imag())>Precision)
		return false;
	return true;
}

inline unsigned int MoveZeroAtEnd(std::vector<std::complex<double>>& x)
{
	unsigned int changed=0;
	for (unsigned int i = 0; i < x.size()-changed && changed<x.size(); i++)
	{
		if(isZero(x[i]))
		{
			std::rotate(x.begin()+i, x.begin() +i+ 1, x.end());
			i--;
			changed+=1;
		}
	}

	return changed;
}
inline void SameTogether(std::vector<std::complex<double>>& x)
{

	for (unsigned int i = 0; i < x.size()-1; i++)
	{
		unsigned int change = 0;
		for (unsigned int j = i + 1; j < x.size(); j++)
		{
			if (isZero(x[i] - x[j]))
			{
				change+=1;
				auto tmp=x[j];
				x[j]=x[i+change];
				x[i+change]=tmp;
			}
		}
	}
}

#include"Polynomial.h"
#include"Matrix.h"

#endif