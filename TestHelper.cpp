#include "TestHelper.h"

#define NUM_PRECIS 10

void PrintMatrix(const Matrix& A, const std::string& str, std::ostream& out)
{
#ifdef NUM_PRECIS	
	out.precision(NUM_PRECIS);
#endif

	out<<"Matrix:"<<str.c_str()<<" "<<A.nrow()<<" "<<A.ncol()<<std::endl;
	for(int i=0;i<A.nrow();i++)
	{
		for(int j=0;j<A.ncol();j++)
			out<<A(i,j)<<" ";
		out<<std::endl;
	}
}

void PrintVector(const Vector& A, const std::string& str, std::ostream& out)
{
#ifdef NUM_PRECIS	
	out.precision(NUM_PRECIS);
#endif

	out<<"Vector:"<<str.c_str()<<" "<<A.size()<<std::endl;
	for(int i=0;i<A.size();i++)
		out<<A(i)<<" ";
	out<<std::endl;
}

void PrintArray(const Array<double>& A, const std::string& str, std::ostream& out)
{	
#ifdef NUM_PRECIS	
	out.precision(NUM_PRECIS);
#endif

	out<<"Array:"<<str.c_str()<<" "<<A.size()<<std::endl;
	for(int i=0;i<A.size();i++)
		out<<A.get(i)<<" ";
	out<<std::endl;
}


void ReadMatrix(Matrix& A, std::string& str, std::istream& fin)
{
	int m,n;
	std::string str_all;
	fin>>str_all>>m>>n;
	str = str_all.substr(7);
	A.resize(m,n);
	for(int i=0;i<m;i++)
		for(int j=0;j<n;j++)
			fin>>A(i,j);
}

void ReadVector(Vector& A, std::string& str, std::istream& fin)
{
	int m;
	std::string str_all;
	fin>>str_all>>m;
	str=str_all.substr(7);
	A.resize(m);
	for(int i=0;i<m;i++)
		fin>>A(i);
}

void PrintMatrixRowAbsSum(const Matrix& A, const std::string& str, std::ostream& out)
{
#ifdef NUM_PRECIS	
	out.precision(NUM_PRECIS);
#endif

	Vector ARowAbsSum(A.nrow(),0.0);

	for(int i=0;i<A.nrow();i++)
	{
		for(int j=0;j<A.ncol();j++)
			ARowAbsSum(i) += fabs(A.get(i,j));
	}

	out<<"Vector:"<<str.c_str()<<"_RowAbsSum "<<ARowAbsSum.size()<<std::endl;
	for(int i=0;i<ARowAbsSum.size();i++)
		out<<ARowAbsSum(i)<<" ";
	out<<std::endl;
}

Vector GetMatrixRowAbsSum(const Matrix& A)
{
	Vector ARowAbsSum(A.nrow(),0.0);

	for(int i=0;i<A.nrow();i++)
	{
		for(int j=0;j<A.ncol();j++)
			ARowAbsSum(i) += fabs(A.get(i,j));
	}

	return ARowAbsSum;

}

void PrintState(const SimTK::State s,std::ostream& out)
{
#ifdef NUM_PRECIS	
	out.precision(NUM_PRECIS);
#endif

	Vector q = s.getQ();
	Vector u = s.getU();
	Vector z = s.getZ();

	out<<"s.time: "<<s.getTime()<<std::endl;
	PrintVector(q,"s.q",out);
	PrintVector(u,"s.u",out);
	PrintVector(z,"s.z",out);
}

void solveLinearCoefs(const Matrix& xy_pair0, const Matrix& xy_pair1, Vector& A, Vector& B)
{
	A = (xy_pair1.col(1)-xy_pair0.col(1)).elementwiseDivide(xy_pair1.col(0)-xy_pair0.col(0));
	B = xy_pair1.col(1)-A.elementwiseMultiply(xy_pair1.col(0));
}

void solveLinearCoefs(const Vector& x0, const Vector& y0, const Vector& x1, const Vector& y1, Vector& A, Vector& B)
{
	Vector delta_x = x1-x0;
	for(int i=0;i<delta_x.size();i++)
	{
		if(delta_x[i] == 0.0)
			delta_x[i] = 1.0;
	}
	A = (y1-y0).elementwiseDivide(delta_x);
	B = y1-A.elementwiseMultiply(x1);
}

void solveLinearCoefs(int n, const double* x0, const double* y0, const double* x1, const double* y1, Vector& A, Vector& B)
{
	A.resize(n);
	B.resize(n);

	for(int i=0;i<n;i++)
	{
		double delta_x = x1[i]-x0[i];
		if(delta_x == 0.0)
			delta_x = 1.0;
		A[i] = (y1[i]-y0[i])/delta_x;
		B[i] = y1[i]-A[i]*x1[i];
	}
}

double clamp(double lower, double upper, double x)
{
	if(x<lower)
		return lower;
	else if(x>upper)
		return upper;
	else
		return x;
}
