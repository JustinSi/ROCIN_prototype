#ifndef OPENSIM_MPC_H_
#define OPENSIM_MPC_H_

/********************************************************/
/*                  MPC Solver                          */
/********************************************************/

#include <SimTKmath.h>
#include <OpenSim/Common/GCVSplineSet.h>

using SimTK::Array_;
using SimTK::Matrix;
using SimTK::Vector;

namespace OpenSim {

class MPC;

// the class that defines an OptimizerSystem which can be used to solve MPC
class MPCQP : public SimTK::OptimizerSystem {
public:
	MPCQP() { _mpc = NULL; }
	MPCQP(MPC* mpc);
	void setMPC(MPC* mpc) { _mpc = mpc; }

	int objectiveFunc( const Vector& coefficients, bool new_coefficients, SimTK::Real& f ) const;
	int gradientFunc( const Vector& coefficients, bool new_coefficients, Vector &gradient ) const;

	int objectiveFunc_const_dynamics_explicit( const Vector& coefficients, bool new_coefficients, SimTK::Real& f ) const;
	int gradientFunc_const_dynamics_explicit( const Vector& coefficients, bool new_coefficients, Vector &gradient ) const;

	int objectiveFunc_const_dynamics_implicit( const Vector& coefficients, bool new_coefficients, SimTK::Real& f ) const;
	int gradientFunc_const_dynamics_implicit( const Vector& coefficients, bool new_coefficients, Vector &gradient ) const;

	int objectiveFunc_varying_dynamics_explicit( const Vector& coefficients, bool new_coefficients, SimTK::Real& f ) const;
	int gradientFunc_varying_dynamics_explicit( const Vector& coefficients, bool new_coefficients, Vector &gradient ) const;

	int objectiveFunc_varying_dynamics_implicit( const Vector& coefficients, bool new_coefficients, SimTK::Real& f ) const;
	int gradientFunc_varying_dynamics_implicit( const Vector& coefficients, bool new_coefficients, Vector &gradient ) const;


private:

	MPC* _mpc;
};

class MPC {
	friend class MPCQP;
public:
	MPC(int n_controls, int n_y, double ti, double dt, const Array_<Vector>& reference, int windowSize = 1);
	MPC(int n_controls, int n_y, int n_samples, double ti, double dt, FunctionSet* r_spline_input, int windowSize = 1);
	

	~MPC();
	
	void setUseImplicit(bool s) { _USE_IMPLICIT = s; }

	void setABC(const Matrix& A, const Matrix& B, const Vector& C);
	void setABCArray(const Array<Matrix>& A_array, const Array<Matrix>& B_array, const Array<Vector>& C_array);

	void setDandE(const Matrix& D, const Vector&E);
	void setDandEArray(const Array<Matrix>& D_array, const Array<Vector>& E_array);

	void setDiagDandE(const Vector& D, const Vector& E);
	void setDiagDandEAarray(const Array<Vector>& D_array, const Array<Vector>& E_array);

	void setPQR(const Matrix&P, const Matrix& Q, const Matrix& R)  { _P=P;_P_is_diag=false; _Q=Q ;_Q_is_diag=false; _R=R;_R_is_diag=false; }
	void setQd(const Matrix& Qd) { _Qd = Qd;_Qd_is_diag=false; }
	void setSd(const Matrix& Sd) { _Sd = Sd;_Sd_is_diag=false; }

	void setDiagPQR(const Vector&P, const Vector&Q, const Vector&R) { _diag_P=P;_P_is_diag=true; _diag_Q=Q ;_Q_is_diag=true; _diag_R=R;_R_is_diag=true; } 
	void setDiagQd(const Vector& Qd) { _diag_Qd = Qd;_Qd_is_diag=true; }
	void setDiagSd(const Vector& Sd) { _diag_Sd = Sd;_Sd_is_diag=true; }


	Vector DLeftMultiply(const Vector& e) const;
	Vector DTransposeLeftMultiply(const Vector& e) const;
	Vector DiLeftMultiply(int i, const Vector& e) const;
	Vector DiTransposeLeftMultiply(int i, const Vector& e) const;

	Vector PLeftMultiply(const Vector& e) const;
	Vector QLeftMultiply(const Vector& e) const;
	Vector RLeftMultiply(const Vector& e) const;
	Vector QdLeftMultiply(const Vector& e) const;
	Vector SdLeftMultiply(const Vector& e) const;


	void setPenalizeYdot(bool s) { _penalize_ydot = s; }
	void setPenalizeUdot(bool s) { _penalize_udot = s; }
	void setUsingVaryingDynamics(bool s) { _use_varying_dynamics = s; }
	void setWinSize(int h) { _qp_win_size = h; _u_array.resize(_r_array[0].size(),h);}
	void setNaturalFrequency(double s) { _natural_frequency = s; }
	int getCurWinSize() { return _qp_win_size; }
	void updateYdotRef(const Vector& curY);
	
	int getTimeIndex(SimTK::Real t);
	void precomputeU(double t, const Vector& initY);
	Vector getCurrentU() const { return _u; }
	Vector getU(int i) const { return _u_array.col(i); }
	const Matrix& getUArray() const { return _u_array; }
	bool isUpToDate() { return _up_to_date; }
	bool isUpToDate(double t);
	double getDt() { return _dt; }
	double getInitTime() { return _initialTime; }
	double getFinalTime() { return _finalTime; }
	Array<double> getReferenceAtTime(double t);

	void setLowerBounds(const Vector& lowerbounds) { _lowerbounds = lowerbounds; }
	void setUpperBounds(const Vector& upperbounds) { _upperbounds = upperbounds; }

	void testMPC();

private:

	double _initialTime;
	double _finalTime;
	double _dt;
	int _n_samples;


	double _natural_frequency;

	Array_<Vector> _r_array;
	FunctionSet* _r_spline;
	Array_<Vector> _r_dot_array;

	//
	// y_dot = A*y + B*u + C;

	Matrix _A;
	Matrix _B;
	Vector _C;

	Matrix _D;
	Vector _diag_D;
	bool _D_is_diag;

	Vector _E;

	Array<Matrix> _A_array;
	Array<Matrix> _B_array;
	Array<Vector> _C_array;

	Array<Matrix> _D_array;
	Array<Vector> _diag_D_array;

	Array<Vector> _E_array;

	// Z = I-A*dt;
	Matrix _Zinv; 

	Array<Matrix> _Zinv_array;

	Array<Matrix> _Z2inv_array;

	// J = y(t)'*P*y(t)+\integ{(y-y_ref)'*Q*(y-y_ref)+(y_dot-y_dot_ref)'*Qd*(y_dot-y_dot_ref)+u'*R*u+u_dot'*Sd*u_dot}

	Matrix _P;
	Vector _diag_P;
	bool _P_is_diag;

	Matrix _Q;
	Vector _diag_Q;
	bool _Q_is_diag;

	Matrix _Qd;
	Vector _diag_Qd;
	bool _Qd_is_diag;

	Matrix _R;
	Vector _diag_R;
	bool _R_is_diag;

	Matrix _Sd;
	Vector _diag_Sd;
	bool _Sd_is_diag;

	MPCQP _qp;
	SimTK::Optimizer _opt;

	int _n_controls;
	int _n_y;

	int _qp_start_index;
	int _qp_win_size;
	Vector _qp_y0;
	Vector _qp_u0;
	Vector _u;
	Vector _u_next;
	Matrix _u_array;

	double _bound_tracking;
	double _bound_control;

	bool _USE_IMPLICIT;

	bool _up_to_date;
	bool _penalize_ydot;
	bool _penalize_udot;
	bool _use_varying_dynamics;

	Vector _lowerbounds;
	Vector _upperbounds;
};

}


#endif