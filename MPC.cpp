#include "MPC.h"
#include "TestHelper.h"
#include <OpenSim/Common/Constant.h>
#include <fstream>
#include <cstring>

using namespace OpenSim;

// using a FunctionSet to provide the desired trajectory of the observation variables
MPC::MPC(int n_controls, int n_y, int n_samples, double ti, double dt, FunctionSet* r_spline_input, int windowSize)
{
	_n_controls = n_controls;
	_n_y = n_y;

	_initialTime = ti;
	_dt = dt;
	_n_samples = n_samples;
	_finalTime = ti+_dt*(_n_samples-1);

	_r_array.resize(_n_samples);
	_r_dot_array.resize(_n_samples);

	_r_spline = new FunctionSet();
	int n_ref_spline_size = r_spline_input->getSize();
	for(int i=0;i<n_ref_spline_size;i++)
		_r_spline->cloneAndAppend(r_spline_input->get(i));

	Constant zeroFunc(0.0);

	for(int i=n_ref_spline_size;i<_n_y;i++)
		_r_spline->cloneAndAppend(zeroFunc);

	Array<double> r_ref_array(0.0,n_y);
	Array<double> r_dot_ref_array(0.0,n_y);

	double time;

	for(int i=0;i<_n_samples;i++)
	{
		_r_array[i].resize(n_y);
		_r_dot_array[i].resize(n_y);

		time = _initialTime + _dt* double(i);
		_r_spline->evaluate(r_ref_array,0,time);
		_r_spline->evaluate(r_dot_ref_array,1,time);
		for(int j=0;j<n_y;j++)
		{
			_r_array[i][j] = r_ref_array[j];
			_r_dot_array[i][j] = r_dot_ref_array[j];
		}

	}

	_USE_IMPLICIT = true;
	_up_to_date = false;
	_penalize_ydot = false;
	_penalize_udot = false;
	_use_varying_dynamics = false;//true;


	_P_is_diag = false;
	_Q_is_diag = false;
	_Qd_is_diag = false;
	_R_is_diag = false;
	_Sd_is_diag = false;

	_qp_start_index = 0;
	_qp_win_size = windowSize;	//10

	_bound_tracking = 0.03;
	_bound_control = 100;

	_natural_frequency = 10.0;

	_A_array.setSize(_qp_win_size);
	_B_array.setSize(_qp_win_size);
	_C_array.setSize(_qp_win_size);
	_Zinv_array.setSize(_qp_win_size);
	_Z2inv_array.setSize(_qp_win_size);

	_D_array.setSize(_qp_win_size);
	_diag_D_array.setSize(_qp_win_size);
	_E_array.setSize(_qp_win_size);


    _u_array.resize(n_controls,_qp_win_size); 
	_u_array.setToZero(); 
	_u.resize(n_controls);
	_u.setToZero();
	_u_next.resize(n_controls);
	_u_next.setToZero();
	_qp_u0.resize(n_controls);
	_qp_u0.setToZero();

	_lowerbounds.resize(n_controls*_qp_win_size);
	_lowerbounds.setTo(-SimTK::Infinity);
	_upperbounds.resize(n_controls*_qp_win_size);
	_upperbounds.setTo(SimTK::Infinity);

	_qp.setMPC(this);
	_qp.setNumParameters(n_controls*_qp_win_size);

	//_opt.setOptimizerSystem(_qp,SimTK::InteriorPoint);
	_opt.setOptimizerSystem(_qp,SimTK::CFSQP);
	//_opt.setOptimizerSystem(_qp,SimTK::LBFGSB);
	_opt.setConvergenceTolerance(1e-4);	//1e-3
	_opt.setMaxIterations(200);	//1000
	_opt.setAdvancedBoolOption("warm_start",true);
	_opt.setAdvancedRealOption("obj_scaling_factor",1);
	_opt.setAdvancedRealOption("nlp_scaling_max_gradient",100);



}

// using an array to provide the desired trajectory of the observation variables
MPC::MPC(int n_controls, int n_y, double ti, double dt, const Array_<Vector>& reference, int windowSize)
{
	_n_controls = n_controls;
	_n_y = n_y;

	_initialTime = ti;
	_dt = dt;
	
	_n_samples = reference.size();
	_r_array.resize(_n_samples);
	int size_refvec = reference[0].size();

	for(int i=0;i<_n_samples;i++)
	{
		_r_array[i].resize(_n_y);
		_r_array[i].setToZero();
		_r_array[i].updBlock(0,0,size_refvec,1) = reference[i];
		
	}

	
	_finalTime = ti+_dt*(_n_samples-1);	
	//compute _r_dot_array;
	_r_dot_array.resize(_n_samples);

	_r_spline = NULL;

	Storage* rStore = new Storage();
	double time;
	for(int i=0;i<_n_samples;i++)
	{
		time = _initialTime +_dt*double(i);
		rStore->append(time,_r_array[i]);
	}

	_r_spline = new GCVSplineSet(5,rStore);
	delete rStore;


	Array<double> rdot_spline(0.0,_n_y);

	for(int i=0;i<_n_samples;i++)
	{

		time = _initialTime + _dt*double(i);
		_r_spline->evaluate(rdot_spline,1,time);

		_r_dot_array[i].resize(_n_y);
		for(int j=0;j<_n_y;j++)
			_r_dot_array[i].set(j,rdot_spline.get(j));

	}


	//_qp = new MPCQP(this);

	_USE_IMPLICIT = true;
	_up_to_date = false;
	_penalize_ydot = false;
	_penalize_udot = false;

	_P_is_diag = false;
	_Q_is_diag = false;
	_Qd_is_diag = false;
	_R_is_diag = false;
	_Sd_is_diag = false;

	_qp_start_index = 0;
	_qp_win_size = windowSize;	//10

	_bound_tracking = 0.03;
	_bound_control = 100;

	_natural_frequency = 10.0;

	_A_array.setSize(_qp_win_size);
	_B_array.setSize(_qp_win_size);
	_C_array.setSize(_qp_win_size);
	_Zinv_array.setSize(_qp_win_size);
	_Z2inv_array.setSize(_qp_win_size);

	_D_array.setSize(_qp_win_size);
	_diag_D_array.setSize(_qp_win_size);
	_E_array.setSize(_qp_win_size);


    _u_array.resize(n_controls,_qp_win_size); 
	_u_array.setToZero(); 
	_u.resize(n_controls);
	_u.setToZero();
	_u_next.resize(n_controls);
	_u_next.setToZero();
	_qp_u0.resize(n_controls);
	_qp_u0.setToZero();

	_lowerbounds.resize(n_controls*_qp_win_size);
	_lowerbounds.setTo(-SimTK::Infinity);
	_upperbounds.resize(n_controls*_qp_win_size);
	_upperbounds.setTo(SimTK::Infinity);

	_qp.setMPC(this);
	_qp.setNumParameters(n_controls*_qp_win_size);

	//_opt.setOptimizerSystem(_qp,SimTK::InteriorPoint);
	_opt.setOptimizerSystem(_qp,SimTK::CFSQP);
	//_opt.setOptimizerSystem(_qp,SimTK::LBFGSB);
	_opt.setConvergenceTolerance(1e-4);	//1e-3
	_opt.setMaxIterations(200);	//1000
	_opt.setAdvancedBoolOption("warm_start",true);
	_opt.setAdvancedRealOption("obj_scaling_factor",1);
	_opt.setAdvancedRealOption("nlp_scaling_max_gradient",100);

}

MPC::~MPC()
{
	if(_r_spline != NULL)
	{
		delete _r_spline;
	}

}

void MPC::setABCArray(const Array<Matrix>& A_array, const Array<Matrix>& B_array, const Array<Vector>& C_array)
{
	for(int i=0;i<A_array.size();i++)
	{

		if(A_array[i].ncol()!=_n_y || A_array[i].nrow()!=_n_y)
		{
			std::cout<<"A_array matrix size do not fit into MPC setting!"<<std::endl;
			exit(0);
		}

		if(B_array[i].ncol()!=_n_controls || B_array[i].nrow()!=_n_y)
		{		
			std::cout<<"B_array matrix size do not fit into MPC setting!"<<std::endl;
			exit(0);
		}

		if(C_array[i].size() != _n_y)
		{
			std::cout<<"C_array vector size do not fit into MPC setting!" <<std::endl;
			exit(0);
		}

		_A_array[i] = A_array[i];
		_B_array[i] = B_array[i];
		_C_array[i] = C_array[i];

		if(_USE_IMPLICIT)
		{

			Matrix Zi(_A_array[i].nrow(),_A_array[i].ncol());
			Zi.setToZero();
			Zi.diag().setTo(1.0);

			Zi -= _A_array[i]*_dt;

			_Zinv_array[i] = Zi.invert();

			Matrix Z2i(_A_array[i].nrow(),_A_array[i].ncol());
			Z2i.setToZero();
			Z2i.diag().setTo(1.0);

			Z2i -= _A_array[i]*(_dt*0.5);
			_Z2inv_array[i] = Z2i.invert();
		}
	}

}

void MPC::setABC(const Matrix& A, const Matrix& B, const Vector& C) 
{ 
	if(A.ncol()!=_n_y || A.nrow()!=_n_y)
	{
		std::cout<<"A size do not fit into MPC setting!"<<std::endl;
		exit(0);
	}

	if(B.ncol()!=_n_controls || B.nrow()!=_n_y)
	{		
		std::cout<<"B size do not fit into MPC setting!"<<std::endl;
		exit(0);
	}

	if(C.size() != _n_y)
	{
		std::cout<<"C size do not fit into MPC setting!" <<std::endl;
		exit(0);
	}


	_A=A;
	_B=B;
	_C=C; 

	if(_USE_IMPLICIT)
	{
		Matrix Z(_A.nrow(),_A.ncol());
		Z.setToZero();
		Z.diag().setTo(1.0);

		Z -= _A*_dt;

		_Zinv = Z.invert();
	}
}

void MPC::setDandE(const Matrix& D, const Vector&E)
{
	if(D.ncol() != _n_controls || D.nrow() != E.size())
	{
		std::cout<<"D and Esize do not fit into MPC setting!"<<std::endl;
		exit(0);
	}

	_D = D;
	_E = E;

	_D_is_diag = false;
}

void MPC::setDiagDandE(const Vector& D, const Vector& E)
{
	if(D.size() != _n_controls || D.size() != E.size())
	{
		std::cout<<"D and Esize do not fit into MPC setting!"<<std::endl;
		exit(0);		
	}

	_diag_D = D;
	_E = E;

	_D_is_diag = true;
}

void MPC::setDandEArray(const Array<Matrix>& D_array, const Array<Vector>& E_array)
{
	for(int i=0;i<D_array.size();i++)
	{
		if(D_array[i].ncol() != _n_controls || D_array[i].nrow() != E_array[i].size())
		{
			std::cout<<"D and Esize do not fit into MPC setting!"<<std::endl;
			exit(0);
		}

		_D_array[i] = D_array[i];
		_E_array[i] = E_array[i];
	}

	_D_is_diag = false;
}

void MPC::setDiagDandEAarray(const Array<Vector>& D_array, const Array<Vector>& E_array)
{
	for(int i=0;i<D_array.size();i++)
	{
		if(D_array[i].size() != _n_controls || D_array[i].size() != E_array[i].size())
		{
			std::cout<<"D and Esize do not fit into MPC setting!"<<std::endl;
			exit(0);
		}

		_diag_D_array[i] = D_array[i];
		_E_array[i] = E_array[i];
	}

	_D_is_diag = true;

}

Vector MPC::DLeftMultiply(const Vector& e) const
{
	if(_D_is_diag)
		return _diag_D.elementwiseMultiply(e);
	else
		return _D*e;
}

Vector MPC::DiLeftMultiply(int i, const Vector& e) const
{
	if(_D_is_diag)
		return _diag_D_array[i].elementwiseMultiply(e);
	else
		return _D_array[i]*e;
}

Vector MPC::DTransposeLeftMultiply(const Vector& e) const
{
	if(_D_is_diag)
		return _diag_D.elementwiseMultiply(e);
	else
		return _D.transpose()*e;
}

Vector MPC::DiTransposeLeftMultiply(int i, const Vector& e) const
{
	if(_D_is_diag)
		return _diag_D_array[i].elementwiseMultiply(e);
	else
		return _D_array[i].transpose()*e;
}

Vector MPC::PLeftMultiply(const Vector& e) const
{
	if(_P_is_diag)
		return _diag_P.elementwiseMultiply(e);
	else
		return _P*e;
}

Vector MPC::QLeftMultiply(const Vector& e) const
{
	if(_Q_is_diag)
		return _diag_Q.elementwiseMultiply(e);
	else
		return _Q*e;
}

Vector MPC::RLeftMultiply(const Vector& e) const
{
	if(_R_is_diag)
		return _diag_R.elementwiseMultiply(e);
	else
		return _R*e;
}

Vector MPC::QdLeftMultiply(const Vector& e) const
{
	if(_Qd_is_diag)
		return _diag_Qd.elementwiseMultiply(e);
	else
		return _Qd*e;
}

Vector MPC::SdLeftMultiply(const Vector& e) const
{
	if(_Sd_is_diag)
		return _diag_Sd.elementwiseMultiply(e);
	else
		return _Sd*e;
}

// test function to see whether the provided gradient functions are correct (compare them with numerical gradients)
void MPC::testMPC()
{


	int n_u = _n_controls;
	int n_y = _n_y;

	int n_vars = n_u*_qp_win_size;


	Vector u(n_vars);
	//u.setToZero();
	//u.setTo(10.0);

	for(int i=0;i<n_vars;i++)
		u.setTo(double(i)*0.1);
		//u.setTo(sin(double(i)*0.1));
	//u.setTo(0.0);
	double delta = 0.000001;
	Vector du(n_vars);


	//	test object and gradient
	Vector gradient_analytic(n_vars);
	Vector gradient_numeric(n_vars);
	//numeric gradient
	double f=0.0,f_new =0.0;
	_qp.objectiveFunc(u,true,f);
	_qp.gradientFunc(u,true,gradient_analytic);
	
	for(int i=0;i<n_vars;i++)
	{
		du.setToZero();
		du(i) = delta;
		_qp.objectiveFunc(u+du,true,f_new);
		double df = f_new-f;
		gradient_numeric(i) = df/delta; 
	}


	Matrix gradient_comparison(n_vars,3);
	gradient_comparison.updCol(0) = gradient_analytic;
	gradient_comparison.updCol(1) = gradient_numeric;
	gradient_comparison.updCol(2) = (gradient_analytic-gradient_numeric).elementwiseDivide(gradient_analytic);


	PrintMatrix(gradient_comparison,"gradient_comparison",std::cout);
	//std::cout<<"gradient_ankle_r: "<<gradient_analytic.get(13)<<" gradient_ankle_l: "<<gradient_analytic.get(22)<<std::endl;
	//PrintVector(gradient_analytic-gradient_numeric,"gradient_difference",std::cout);


	//test constraint and jacobian

/*	int n_c = 0;
	n_c = n_y;//6*_qp_win_size;

	Matrix Jacob_analytic(n_c,n_vars);
	Matrix Jacob_numeric(n_c,n_vars);

	Vector c(n_c),c_new(n_c);
	_qp.constraintFunc(u,true,c);
	_qp.constraintJacobian(u,true,Jacob_analytic);

	for(int i=0;i<n_vars;i++)
	{
		du.setToZero();
		du(i) = delta;
		_qp.constraintFunc(u+du,true,c_new);
		Jacob_numeric.updCol(i) = (c_new-c)/delta;
	}

	Matrix Jacob_diff = Jacob_analytic - Jacob_numeric;
	Matrix Jacob_diff_normalize = Jacob_diff.elementwiseDivide(Jacob_analytic);

	PrintMatrix(Jacob_diff,"Jacob_diff",std::cout);
	PrintMatrix(Jacob_diff_normalize,"Jacob_diff_normalize",std::cout);
	*/

}

int MPC::getTimeIndex(SimTK::Real t)
{
	int k= round((t-_initialTime)/_dt);
	if(k<0)
		k=0;
	if(k>=_n_samples)
		k=_n_samples-1;

	return k;
}

bool MPC::isUpToDate(double t)
{
	int cur_k = getTimeIndex(t);
	if(cur_k > _qp_start_index)
		return false;

	return true;
}

// this is the core function that solves MPC problem
void MPC::precomputeU(double t, const Vector& initY)
{
	int n_u = 0;
	if(_B_array.size()>0)
		n_u = _B_array[0].ncol();
	else
		n_u = _B.ncol();

	int n_y = initY.size();


	int cur_k = getTimeIndex(t);
	if(cur_k > _qp_start_index)
		_up_to_date = false;

	

	_qp_y0 = initY;
	_qp_start_index = cur_k;
	if(_qp_win_size>_n_samples-1-cur_k)
		_qp_win_size = _n_samples-1-cur_k;

	if(_qp_win_size<1)
	{
		_up_to_date = true;
		return;
	}

	//using PD law to update Ydot reference
	if(_penalize_ydot)
		updateYdotRef(initY);

	int n_para = n_u*_qp_win_size;		

	_qp.setNumParameters(n_para);
	_qp.setNumEqualityConstraints(0);
	_qp.setNumInequalityConstraints(0);

	_qp.setParameterLimits(_lowerbounds.block(0,0,n_para,1).getAsVector(),_upperbounds.block(0,0,n_para,1).getAsVector());

	Vector result_u(n_para);


	result_u.setTo(0.0);

	for(int i=0;i<_qp_win_size;i++)
	{
		result_u.updBlock(n_u*i,0,n_u,1) = _u_array.col(i);
	}
	
	SimTK::Real f = 0.0;
	

	try{
		f = _opt.optimize(result_u);
	}
	catch(const std::exception& ex)
	{
		std::cout<<ex.what()<<std::endl;
	}

	std::cout<<"t = "<<t<<std::endl;
	std::cout<<"optimization error: "<<f<<std::endl;




	for(int i=0;i<_qp_win_size;i++)
	{
		_u_array.updCol(i) = result_u.block(n_u*i,0,n_u,1).getAsVector();
	}

	_u = _u_array.col(0);
	if(_qp_win_size>1)
		_u_next = _u_array.col(1);
	else
		_u_next = _u;

	_qp_u0 = _u;

	_up_to_date = true;

}

// update the YdotRef (in our case, it is mainly used to update the desired acceleration by PD rule)
void MPC::updateYdotRef(const Vector& curY)
{

	//critical damping, natural frequency
	double omega = _natural_frequency;

	//std::cout<<"omega = "<<omega<<std::endl;

	double Ks = omega*omega, Kd = 2.0*omega;

	int n_y = curY.size();

	int n_coords = n_y/2;
	Matrix acc_exp(n_coords,1);
	Matrix vel_exp(n_coords,1);

	Array<double> rdot_spline;
	Array<double> rddot_spline;


	int k_ydot = 0;
	int k_cur = 0;
	double time = 0.0;

	Vector pErr_init, vErr_init;
	Vector A, B;

	for(int i=0;i<_qp_win_size;i++)	
	{
		k_cur = _qp_start_index+i;
		if(_USE_IMPLICIT)
			k_ydot = k_cur+1;
		else
			k_ydot = k_cur;

		time = _initialTime + _dt*double(k_cur+1);

		_r_spline->evaluate(rdot_spline,1,time);

		for(int j=0;j<n_coords;j++)
		{
			vel_exp.set(j,0,rdot_spline.get(j));
			acc_exp.set(j,0,rdot_spline.get(j+n_coords));
		}


		if(i==0)
		{
			pErr_init = (_r_array[k_cur].block(0,0,n_coords,1)-curY.block(0,0,n_coords,1)).getAsVector();
			vErr_init = (_r_array[k_cur].block(n_coords,0,n_coords,1)-curY.block(n_coords,0,n_coords,1)).getAsVector();
			A = pErr_init;
			B = vErr_init+pErr_init*omega;
		}

		double t_damp = _dt*i;
		double exp_damp = exp(-omega*t_damp);

		Vector pErr = (A+B*t_damp)*exp_damp;
		Vector vErr = (B*(1.0-omega*t_damp)-A*omega)*exp_damp;

		_r_dot_array[k_ydot].updBlock(0,0,n_coords,1) = vel_exp+Ks*pErr;
		_r_dot_array[k_ydot].updBlock(n_coords,0,n_coords,1) = acc_exp+Ks*pErr+Kd*vErr;
		
	}

}

MPCQP::MPCQP(MPC* mpc)
{
	_mpc = mpc;
}

// override the objective function
int MPCQP::objectiveFunc( const Vector& coefficients, bool new_coefficients, SimTK::Real& f ) const
{
	if(_mpc->_USE_IMPLICIT)
	{
		if(_mpc->_use_varying_dynamics)
			return objectiveFunc_varying_dynamics_implicit(coefficients,new_coefficients,f);
		else
			return objectiveFunc_const_dynamics_implicit(coefficients,new_coefficients,f);
	}
	else
	{
		if(_mpc->_use_varying_dynamics)
			return objectiveFunc_varying_dynamics_explicit(coefficients,new_coefficients,f);
		else
			return objectiveFunc_const_dynamics_explicit(coefficients,new_coefficients,f);

	}
}

// override the gradient of the objective function
int MPCQP::gradientFunc( const Vector& coefficients, bool new_coefficients, Vector &gradient ) const
{
	if(_mpc->_USE_IMPLICIT)
	{
		if(_mpc->_use_varying_dynamics)
			return gradientFunc_varying_dynamics_implicit(coefficients,new_coefficients,gradient);
		else
			return gradientFunc_const_dynamics_implicit(coefficients,new_coefficients,gradient);
	}
	else
	{
		if(_mpc->_use_varying_dynamics)
			return gradientFunc_varying_dynamics_explicit(coefficients,new_coefficients,gradient);
		else
			return gradientFunc_const_dynamics_explicit(coefficients,new_coefficients,gradient);
	}
}

// using varying dynamics (instead of constant dynamics) and explicit formulation
int MPCQP::objectiveFunc_varying_dynamics_explicit( const Vector& coefficients, bool new_coefficients, SimTK::Real& f ) const
{
	f = 0.0;
	Vector y = _mpc->_qp_y0;
	Vector u_pre = _mpc->_qp_u0;
	int n_y = y.size();
	int n_u = _mpc->_B_array[0].ncol();

	int k_y = 0;

	for(int i=0;i<_mpc->_qp_win_size;i++)
	{
		int k_u = _mpc->_qp_start_index+i;
		k_y = k_u+1;
		Vector u(n_u);

		u = coefficients.block(i*n_u,0,n_u,1).getAsVector();

		Vector ydot = _mpc->_A_array[i]*y+_mpc->_B_array[i]*u+_mpc->_C_array[i];

		y += ydot*_mpc->_dt;
		Vector e = y-_mpc->_r_array[k_y];
		if(_mpc->_penalize_ydot)
		{
			Vector e_ydot = ydot - _mpc->_r_dot_array[k_u];
			f += e_ydot.elementwiseMultiply(_mpc->QdLeftMultiply(e_ydot)).sum();
		}

		f += e.elementwiseMultiply(_mpc->QLeftMultiply(e)).sum();

		Vector z = _mpc->DiLeftMultiply(i,u)+_mpc->_E_array[i];
		f += z.elementwiseMultiply(_mpc->RLeftMultiply(z)).sum();

		if(_mpc->_penalize_udot)
		{
			Vector udot = (u-u_pre)/_mpc->_dt;
			f += udot.elementwiseMultiply(_mpc->SdLeftMultiply(udot)).sum();
			u_pre = u;
		}
	}

	f *= _mpc->_dt;

	Vector e = y-_mpc->_r_array[k_y];

	f += e.elementwiseMultiply(_mpc->PLeftMultiply(e)).sum();

	return 0;

}

int MPCQP::gradientFunc_varying_dynamics_explicit(const Vector& coefficients, bool new_coefficients, Vector &gradient ) const
{
	gradient.setToZero();
	Vector y = _mpc->_qp_y0;
	Vector u_pre = _mpc->_qp_u0;
	int n_y = y.size();
	int n_u = _mpc->_B_array[0].ncol();

	int k_y = 0;

	Matrix dydu(n_y,n_u);
	dydu.setToZero();
	Matrix dydotdu(n_y,n_u);
	dydotdu.setToZero();


	Matrix eye(n_y,n_y);
	eye.setToZero();
	eye.diag().setTo(1.0);


	for(int i=0;i<_mpc->_qp_win_size;i++)
	{
		int k_u = _mpc->_qp_start_index+i;
		k_y = k_u+1;

		Vector u(n_u);
		u = coefficients.block(i*n_u,0,n_u,1).getAsVector();
		Vector ydot = _mpc->_A_array[i]*y+_mpc->_B_array[i]*u+_mpc->_C_array[i];
		y += ydot*_mpc->_dt;

		Vector e_ydot = ydot - _mpc->_r_dot_array[k_u];
		Vector e_y = y - _mpc->_r_array[k_y];

		for(int j=i;j>=0;j--)
		{
			if(j==i)
			{
				dydotdu = _mpc->_B_array[i];
				dydu = dydotdu*_mpc->_dt;
			}
			else
			{
				//dydotdu = _mpc->_B_array[j]*
				dydotdu = _mpc->_B_array[j]*_mpc->_dt;
				for(int k=j+1;k<i;k++)
					dydotdu = (eye+_mpc->_A_array[k]*_mpc->_dt)*dydotdu;
				
				dydu = (eye+_mpc->_A_array[i]*_mpc->_dt)*dydotdu;
				dydotdu = _mpc->_A_array[i]*dydotdu;
			}


			gradient.updBlock(j*n_u,0,n_u,1) += 2.0*dydu.transpose()*_mpc->QLeftMultiply(e_y);

			if(_mpc->_penalize_ydot)
				gradient.updBlock(j*n_u,0,n_u,1) += 2.0*dydotdu.transpose()*_mpc->QdLeftMultiply(e_ydot);			

			if(i==_mpc->_qp_win_size-1)
				gradient.updBlock(j*n_u,0,n_u,1) += 2.0*dydu.transpose()*_mpc->PLeftMultiply(e_y)/_mpc->_dt;

		}

		Vector z = _mpc->DiLeftMultiply(i,u)+_mpc->_E_array[i];
		gradient.updBlock(i*n_u,0,n_u,1) += _mpc->DiTransposeLeftMultiply(i,2.0*_mpc->RLeftMultiply(z));

		if(_mpc->_penalize_udot)
		{
			Vector udot = (u-u_pre)/_mpc->_dt;
			Vector temp = 2.0*_mpc->SdLeftMultiply(udot)/_mpc->_dt;

			if(i==0)
				gradient.updBlock(i*n_u,0,n_u,1) += temp;
			else
			{
				gradient.updBlock(i*n_u,0,n_u,1) += temp;
				gradient.updBlock((i-1)*n_u,0,n_u,1) -= temp;
			}

			u_pre = u;
		}


	}

	gradient *= _mpc->_dt;


	return 0;

}

// using varying dynamics and implicit formulation
int MPCQP::objectiveFunc_varying_dynamics_implicit( const Vector& coefficients, bool new_coefficients, SimTK::Real& f ) const
{
	f = 0.0;
	Vector y = _mpc->_qp_y0;
	Vector u_pre = _mpc->_qp_u0;
	int n_y = y.size();
	int n_u = _mpc->_B.ncol();

	int k_y = 0;

	for(int i=0;i<_mpc->_qp_win_size;i++)
	{
		int k_u = _mpc->_qp_start_index+i;
		k_y = k_u+1;
		Vector u(n_u);

		u = coefficients.block(i*n_u,0,n_u,1).getAsVector();

		Vector tmp_BuPlusC = _mpc->_B_array[i]*u+_mpc->_C_array[i];

		//y = _mpc->_Zinv_array[i]*(y+_mpc->_B_array[i]*_mpc->_dt*u+_mpc->_C_array[i]*_mpc->_dt);
		y = _mpc->_Zinv_array[i]*(y+tmp_BuPlusC*_mpc->_dt);

		//Vector ydot = _mpc->_A_array[i]*y+_mpc->_B_array[i]*u+_mpc->_C_array[i];
		Vector ydot = _mpc->_A_array[i]*y+tmp_BuPlusC;

		Vector e = y-_mpc->_r_array[k_y];
		if(_mpc->_penalize_ydot)
		{
			Vector e_ydot = ydot - _mpc->_r_dot_array[k_y];
			f += e_ydot.elementwiseMultiply(_mpc->QdLeftMultiply(e_ydot)).sum();
		}

		f += e.elementwiseMultiply(_mpc->QLeftMultiply(e)).sum();

		Vector z = _mpc->DiLeftMultiply(i,u)+_mpc->_E_array[i];
		f += z.elementwiseMultiply(_mpc->RLeftMultiply(z)).sum();

		if(_mpc->_penalize_udot)
		{
			Vector udot = (u-u_pre)/_mpc->_dt;
			f += udot.elementwiseMultiply(_mpc->SdLeftMultiply(udot)).sum();
			u_pre = u;
		}

	}

	f *= _mpc->_dt;

	Vector e = y-_mpc->_r_array[k_y];

	f += e.elementwiseMultiply(_mpc->PLeftMultiply(e)).sum();

	return 0;

}

int MPCQP::gradientFunc_varying_dynamics_implicit( const Vector& coefficients, bool new_coefficients, Vector &gradient ) const
{
	gradient.setToZero();
	Vector y = _mpc->_qp_y0;
	Vector u_pre = _mpc->_qp_u0;
	int n_y = y.size();
	int n_u = _mpc->_B_array[0].ncol();

	int k_y = 0;
	Matrix eye(n_y,n_y);
	eye.setToZero();
	eye.diag().setTo(1.0);

	Matrix dydu(n_y,n_u);
	dydu.setToZero();
	Matrix dydotdu(n_y,n_u);
	dydotdu.setToZero();

	for(int i=0;i<_mpc->_qp_win_size;i++)
	{
		int k_u = _mpc->_qp_start_index+i;
		k_y = k_u+1;

		Vector u(n_u);
		u = coefficients.block(i*n_u,0,n_u,1).getAsVector();

		Vector tmp_BuPlusC = _mpc->_B_array[i]*u+_mpc->_C_array[i];

		//y = _mpc->_Zinv_array[i]*(y+_mpc->_B_array[i]*_mpc->_dt*u+_mpc->_C_array[i]*_mpc->_dt);
		y = _mpc->_Zinv_array[i]*(y+tmp_BuPlusC*_mpc->_dt);

		//Vector ydot = _mpc->_A_array[i]*y+_mpc->_B_array[i]*u+_mpc->_C_array[i];
		Vector ydot = _mpc->_A_array[i]*y+tmp_BuPlusC;

		Vector e_y = y-_mpc->_r_array[k_y];
		Vector e_ydot = ydot-_mpc->_r_dot_array[k_y];

		for(int j=i;j>=0;j--)
		{

			dydu = _mpc->_Zinv_array[j]*_mpc->_B_array[j]*_mpc->_dt;
			for(int k=j+1;k<=i;k++)
				dydu = _mpc->_Zinv_array[k]*dydu;

			gradient.updBlock(j*n_u,0,n_u,1) += 2.0*dydu.transpose()*_mpc->QLeftMultiply(e_y);

			if(_mpc->_penalize_ydot)
			{
				if(j==i)
					dydotdu = _mpc->_A_array[i]*dydu + _mpc->_B_array[i];
				else
					dydotdu = _mpc->_A_array[i]*dydu;

				gradient.updBlock(j*n_u,0,n_u,1) += 2.0*dydotdu.transpose()*_mpc->QdLeftMultiply(e_ydot);
			}

			if(i==_mpc->_qp_win_size-1)
				gradient.updBlock(j*n_u,0,n_u,1) += 2.0*dydu.transpose()*_mpc->PLeftMultiply(e_y)/_mpc->_dt;
		}

		Vector z = _mpc->DiLeftMultiply(i,u)+_mpc->_E_array[i];
		gradient.updBlock(i*n_u,0,n_u,1) += _mpc->DiTransposeLeftMultiply(i,2.0*_mpc->RLeftMultiply(z));


		if(_mpc->_penalize_udot)
		{
			Vector udot = (u-u_pre)/_mpc->_dt;
			Vector temp = 2.0*_mpc->SdLeftMultiply(udot)/_mpc->_dt;

			if(i==0)
				gradient.updBlock(i*n_u,0,n_u,1) += temp;
			else
			{
				gradient.updBlock(i*n_u,0,n_u,1) += temp;
				gradient.updBlock((i-1)*n_u,0,n_u,1) -= temp;
			}

			u_pre = u;
		}

	}

	gradient *= _mpc->_dt;

	return 0;
}

// using constant dynamics and explicit formulation
int MPCQP::objectiveFunc_const_dynamics_explicit( const Vector& coefficients, bool new_coefficients, SimTK::Real& f ) const
{
//	PrintVector(coefficients,"coefficients",std::cout);
//	PrintVector(coefficients.block(92,0,39,1).getAsVector(),"coefficients_subblock",std::cout);

	f=0.0;
	Vector y = _mpc->_qp_y0;
	Vector u_pre = _mpc->_qp_u0;
	int n_y = y.size();
	int n_u = _mpc->_B.ncol();

	int k_y = 0;


	for(int i=0;i<_mpc->_qp_win_size;i++)
	{		
		int k_u = _mpc->_qp_start_index+i;
		k_y = k_u+1;
		Vector u(n_u);
		
        u = coefficients.block(i*n_u,0,n_u,1).getAsVector();

		Vector ydot = _mpc->_A*y+_mpc->_B*u+_mpc->_C;
		
		y += ydot*_mpc->_dt;
		Vector e = y-_mpc->_r_array[k_y];
		if(_mpc->_penalize_ydot)
		{
			Vector e_ydot  = ydot - _mpc->_r_dot_array[k_u];
			f += e_ydot.elementwiseMultiply(_mpc->QdLeftMultiply(e_ydot)).sum();
		}
		f += e.elementwiseMultiply(_mpc->QLeftMultiply(e)).sum();

		Vector z = _mpc->DLeftMultiply(u)+_mpc->_E;
		f += z.elementwiseMultiply(_mpc->RLeftMultiply(z)).sum();
		if(_mpc->_penalize_udot)
		{
			Vector udot = (u-u_pre)/_mpc->_dt;
			f += udot.elementwiseMultiply(_mpc->SdLeftMultiply(udot)).sum();
			u_pre = u;
		}

		
	}

	f *= _mpc->_dt;
	
	Vector e = y-_mpc->_r_array[k_y];

	f += e.elementwiseMultiply(_mpc->PLeftMultiply(e)).sum();

//	std::cout<<"*********************************************************"<<std::endl;
//	std::cout<<"x: "<<coefficients<<std::endl;
//	std::cout<<"objVal: "<<f<<std::endl;
//	std::cout<<"*********************************************************"<<std::endl;


	return 0;
}


int MPCQP::gradientFunc_const_dynamics_explicit( const Vector& coefficients, bool new_coefficients, Vector &gradient ) const
{

	gradient.setToZero();
	Vector y = _mpc->_qp_y0;
	Vector u_pre = _mpc->_qp_u0;
	int n_y = y.size();
	int n_u = _mpc->_B.ncol();

	

		int k_y = 0;
		Matrix eye(n_y,n_y);
		eye.setToZero();
		eye.diag().setTo(1.0);

		Matrix AdtPlusI = _mpc->_A*_mpc->_dt+eye;

		Matrix dydu(_mpc->_B.nrow(),_mpc->_B.ncol());
		dydu.setToZero();
		Matrix dydotdu(_mpc->_B.nrow(),_mpc->_B.ncol());
		dydotdu.setToZero();

//		PrintMatrixRowAbsSum(_mpc->_Q,"_Q",std::cout);
//		PrintMatrixRowAbsSum(_mpc->_B.transpose(),"_B.transpose",std::cout);
//		PrintMatrixRowAbsSum(_mpc->_R,"_R",std::cout);

		for(int i=0;i<_mpc->_qp_win_size;i++)
		{
			int k_u = _mpc->_qp_start_index+i;
			k_y = k_u+1;
			Vector u(n_u);
			u = coefficients.block(i*n_u,0,n_u,1).getAsVector();
			Vector ydot = _mpc->_A*y+_mpc->_B*u+_mpc->_C;			
			y += ydot*_mpc->_dt;

			if(_mpc->_penalize_ydot)
			{
				Vector e_ydot  = ydot - _mpc->_r_dot_array[k_u];
				Matrix tmp(_mpc->_A.nrow(),_mpc->_A.ncol(),0.0);
				tmp.diag().setTo(1.0);
				for(int j=i;j>=0;j--)
				{
					if(j==i)
						dydotdu = _mpc->_B;
					else if(j==i-1)
					{
						tmp = _mpc->_A*_mpc->_dt;
						dydotdu = tmp*_mpc->_B;
					}
					else
					{
						tmp = tmp*AdtPlusI;
						dydotdu = tmp*_mpc->_B;
					}

					gradient.updBlock(j*n_u,0,n_u,1) += 2.0*dydotdu.transpose()*_mpc->QdLeftMultiply(e_ydot);

				}
			}

			for(int j=i;j>=0;j--)
			{
				if(j==i)
					dydu = _mpc->_B*_mpc->_dt;
				else
					dydu = AdtPlusI*dydu;

				
//				PrintMatrixRowAbsSum(dydu.transpose(),"dydu.transpose",std::cout);

				gradient.updBlock(j*n_u,0,n_u,1) += 2.0*dydu.transpose()*_mpc->QLeftMultiply(y-_mpc->_r_array[k_y]);

//				PrintVector(gradient.block(92,0,39,1).getAsVector(),"gradient_subvector",std::cout);
			}

			Vector z = _mpc->DLeftMultiply(u)+_mpc->_E;
			gradient.updBlock(i*n_u,0,n_u,1) += _mpc->DTransposeLeftMultiply(2.0*_mpc->RLeftMultiply(z));


			if(_mpc->_penalize_udot)
			{
				Vector udot = (u-u_pre)/_mpc->_dt;
				Vector temp = 2.0*_mpc->SdLeftMultiply(udot)/_mpc->_dt;

				if(i==0)
					gradient.updBlock(i*n_u,0,n_u,1) += temp;
				else
				{
					gradient.updBlock(i*n_u,0,n_u,1) += temp;
					gradient.updBlock((i-1)*n_u,0,n_u,1) -= temp;
				}

//				PrintVector(gradient.block(92,0,39,1).getAsVector(),"gradient_subvector",std::cout);

				u_pre = u;
			}
		}

		gradient *= _mpc->_dt;
		//gradient *= _mpc->_dt/double(_mpc->_qp_win_size); //decouple the window size effect;

		Vector e = y-_mpc->_r_array[k_y];
		for(int i=_mpc->_qp_win_size-1;i>=0;i--)
		{
			if(i==_mpc->_qp_win_size-1)
				dydu = _mpc->_B*_mpc->_dt;
			else
				dydu = AdtPlusI*dydu;


			gradient.updBlock(i*n_u,0,n_u,1) += 2.0*dydu.transpose()*_mpc->PLeftMultiply(e);

		}


	return 0;
}

// using constant dynamics and implicit formulation
int MPCQP::objectiveFunc_const_dynamics_implicit( const Vector& coefficients, bool new_coefficients, SimTK::Real& f ) const
{
//	PrintVector(coefficients,"coefficients",std::cout);
//	PrintVector(coefficients.block(92,0,39,1).getAsVector(),"coefficients_subblock",std::cout);

	f=0.0;
	Vector y = _mpc->_qp_y0;
	Vector u_pre = _mpc->_qp_u0;
	int n_y = y.size();
	int n_u = _mpc->_B.ncol();

	int k_y = 0;


	for(int i=0;i<_mpc->_qp_win_size;i++)
	{		
		int k_u = _mpc->_qp_start_index+i;
		k_y = k_u+1;
		Vector u(n_u);

		u = coefficients.block(i*n_u,0,n_u,1).getAsVector();

		Vector tmp_BuPlusC = _mpc->_B*u+_mpc->_C;

		//y = _mpc->_Zinv*(y+_mpc->_B*_mpc->_dt*u+_mpc->_C*_mpc->_dt);
		y = _mpc->_Zinv*(y+tmp_BuPlusC*_mpc->_dt);

		//Vector ydot = _mpc->_A*y+_mpc->_B*u+_mpc->_C;
		Vector ydot = _mpc->_A*y+tmp_BuPlusC;
		
		//y += ydot*_mpc->_dt;
		Vector e = y-_mpc->_r_array[k_y];
		if(_mpc->_penalize_ydot)
		{
			Vector e_ydot  = ydot - _mpc->_r_dot_array[k_y];
			f += e_ydot.elementwiseMultiply(_mpc->QdLeftMultiply(e_ydot)).sum();
		}
		f += e.elementwiseMultiply(_mpc->QLeftMultiply(e)).sum();

		Vector z = _mpc->DLeftMultiply(u)+_mpc->_E;
		f += z.elementwiseMultiply(_mpc->RLeftMultiply(z)).sum();

		if(_mpc->_penalize_udot)
		{
			Vector udot = (u-u_pre)/_mpc->_dt;
			f += udot.elementwiseMultiply(_mpc->SdLeftMultiply(udot)).sum();
			u_pre = u;
		}

		
	}

	f *= _mpc->_dt;
	//f *= _mpc->_dt/double(_mpc->_qp_win_size);	//decouple the window size effect;
	
	Vector e = y-_mpc->_r_array[k_y];

	f += e.elementwiseMultiply(_mpc->PLeftMultiply(e)).sum();

//	std::cout<<"objective: "<<f<<std::endl;

	return 0;
}

int MPCQP::gradientFunc_const_dynamics_implicit( const Vector& coefficients, bool new_coefficients, Vector &gradient ) const
{

	gradient.setToZero();
	Vector y = _mpc->_qp_y0;
	Vector u_pre = _mpc->_qp_u0;
	int n_y = y.size();
	int n_u = _mpc->_B.ncol();

	

	int k_y = 0;
	Matrix eye(n_y,n_y);
	eye.setToZero();
	eye.diag().setTo(1.0);

	Matrix AdtPlusI = _mpc->_A*_mpc->_dt+eye;

	Matrix dydu(_mpc->_B.nrow(),_mpc->_B.ncol());
	dydu.setToZero();
	Matrix dydotdu(_mpc->_B.nrow(),_mpc->_B.ncol());
	dydotdu.setToZero();


	for(int i=0;i<_mpc->_qp_win_size;i++)
	{
		int k_u = _mpc->_qp_start_index+i;
		k_y = k_u+1;
		Vector u(n_u);
		u = coefficients.block(i*n_u,0,n_u,1).getAsVector();

		Vector tmp_BuPlusC = _mpc->_B*u+_mpc->_C;

		//y = _mpc->_Zinv*(y+_mpc->_B*_mpc->_dt*u+_mpc->_C*_mpc->_dt);
		y = _mpc->_Zinv*(y+tmp_BuPlusC*_mpc->_dt);
		
		//Vector ydot = _mpc->_A*y+_mpc->_B*u+_mpc->_C;
		Vector ydot = _mpc->_A*y+tmp_BuPlusC;

		Vector e_y = y-_mpc->_r_array[k_y];
		Vector e_ydot  = ydot - _mpc->_r_dot_array[k_y];

		for(int j=i;j>=0;j--)
		{
			if(j==i)
				dydu = _mpc->_Zinv*_mpc->_B*_mpc->_dt;
			else
				dydu = _mpc->_Zinv*dydu;

			gradient.updBlock(j*n_u,0,n_u,1) += 2.0*dydu.transpose()*_mpc->QLeftMultiply(e_y);

			if(_mpc->_penalize_ydot)
			{
				if(j==i)
					dydotdu = _mpc->_A*dydu+_mpc->_B;
				else
					dydotdu = _mpc->_A*dydu;

				gradient.updBlock(j*n_u,0,n_u,1) += 2.0*dydotdu.transpose()*_mpc->QdLeftMultiply(e_ydot);
			}

		}

		Vector z = _mpc->DLeftMultiply(u)+_mpc->_E;
		gradient.updBlock(i*n_u,0,n_u,1) += _mpc->DTransposeLeftMultiply(2.0*_mpc->RLeftMultiply(z));


		if(_mpc->_penalize_udot)
		{
			Vector udot = (u-u_pre)/_mpc->_dt;
			Vector temp = 2.0*_mpc->SdLeftMultiply(udot)/_mpc->_dt;

			if(i==0)
				gradient.updBlock(i*n_u,0,n_u,1) += temp;
			else
			{
				gradient.updBlock(i*n_u,0,n_u,1) += temp;
				gradient.updBlock((i-1)*n_u,0,n_u,1) -= temp;
			}

//				PrintVector(gradient.block(92,0,39,1).getAsVector(),"gradient_subvector",std::cout);

			u_pre = u;
		}
	}

	gradient *= _mpc->_dt;
	//gradient *= _mpc->_dt/double(_mpc->_qp_win_size); //decouple the window size effect;

	Vector e = y-_mpc->_r_array[k_y];
	for(int i=_mpc->_qp_win_size-1;i>=0;i--)
	{
		if(i==_mpc->_qp_win_size-1)
			dydu =  _mpc->_Zinv*_mpc->_B*_mpc->_dt;
		else
			dydu = _mpc->_Zinv*dydu;

		gradient.updBlock(i*n_u,0,n_u,1) += 2.0*dydu.transpose()*_mpc->PLeftMultiply(e);
	}

//	PrintVector(gradient,"gradient",std::cout);
//	PrintVector(gradient.block(92,0,39,1).getAsVector(),"gradient_subvector",std::cout);

/*	Matrix dfdF(n_u,n_y);
	Vector dfdv(n_u);

	memcpy(dfdF.updContiguousScalarData(),gradient.getContiguousScalarData(),sizeof(SimTK::Real)*n_u*n_y);
	memcpy(dfdv.updContiguousScalarData(),gradient.getContiguousScalarData()+n_u*n_y,sizeof(SimTK::Real)*n_u);

	PrintMatrix(dfdF.block(0,9,n_u,3),"dfdF_part",std::cout);
	PrintVector(dfdv,"dfdv",std::cout);*/

	return 0;
}

Array<double> MPC::getReferenceAtTime(double t)
{
	int size_refvec = _r_array[0].size();
	Array<double> r(0.0,size_refvec);

	int cur_k = getTimeIndex(t);

	for(int i=0;i<size_refvec;i++)
		r.set(i,_r_array[cur_k].get(i));

	return r;
}
