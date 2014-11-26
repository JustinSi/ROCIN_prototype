#ifndef OPENSIM_TESTHELPER_H
#define OPENSIM_TESTHELPER_H
#include "SimTKmath.h"
#include <OpenSim\Common\Array.h>
using SimTK::Matrix;
using SimTK::Vector;
using SimTK::Vec3;
using OpenSim::Array;

void PrintMatrix(const Matrix& A, const std::string& str, std::ostream& out);
void PrintMatrixRowAbsSum(const Matrix& A, const std::string& str, std::ostream& out);
void PrintVector(const Vector& A, const std::string& str, std::ostream& out);
void ReadMatrix(Matrix& A, std::string& str, std::istream& fin);
void ReadVector(Vector& A, std::string& str, std::istream& fin);
Vector GetMatrixRowAbsSum(const Matrix& A);
void PrintState(const SimTK::State s, std::ostream& out);

void PrintArray(const Array<double>& A, const std::string& str, std::ostream& out);
void solveLinearCoefs(const Matrix& xy_pair0, const Matrix& xy_pair1, Vector& A, Vector& B);
void solveLinearCoefs(const Vector& x0, const Vector& y0, const Vector& x1, const Vector& y1, Vector& A, Vector& B);
void solveLinearCoefs(int n, const double* x0, const double* y0, const double* x1, const double* y1, Vector& A, Vector& B);

double clamp(double lower, double upper, double x);
//int round(double d);
#endif