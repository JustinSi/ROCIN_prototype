#include "StateInitializer.h"
#include <OpenSim/Simulation/Model/Muscle.h>
#include "ROCINForceController.h"
#include "TestHelper.h"

StateInitializer::StateInitializer(Model* model, ConstantController* controller, const SimTK::State& initState)
{
	_model = model;
	_controller = controller;
	_initState = initState;

	_n_controls = _model->getActuators().getSize();

	_n_muscles = 0;

	for(int i=0;i<_n_controls;i++)
	{
		Actuator& act = _model->getActuators().get(i);
		Muscle* m = dynamic_cast<Muscle*>(&act);
		if(m!= NULL)
			_n_muscles++;
	}


	_n_Z = _initState.getZ().size();
	_w_controls.resize(_n_controls);
	_w_actuator_forces.resize(_n_controls);
	

	_w_controls.setTo(1.0);
	_w_actuator_forces.setTo(1.0);
	_w_acts = 1.0;
	_w_acc  = 1.0e1;//1.0e3;

	_lower_bd_actuatorforces.resize(_n_controls);
	_upper_bd_actuatorforces.resize(_n_controls);

	_lower_bd_actuatorforces.setTo(-SimTK::Infinity);
	_upper_bd_actuatorforces.setTo(SimTK::Infinity);
}


SimTK::State StateInitializer::getOptimizedInitState()
{
	if(_n_Z == 0)
		return _initState;


	StateInitializerQP qp(this);

	qp.setNumParameters(_n_controls);

	PrintVector(_lower_bd_actuatorforces,"_lower_bd_actuatorforces",std::cout);
	PrintVector(_upper_bd_actuatorforces,"_upper_bd_actuatorforces",std::cout);

	qp.setParameterLimits(_lower_bd_actuatorforces,_upper_bd_actuatorforces);
	qp.setNumEqualityConstraints(_UDotRef.size());
	qp.setNumInequalityConstraints(0);
	

	//SimTK::Optimizer opt(qp,SimTK::CFSQP);
	SimTK::Optimizer opt(qp, SimTK::InteriorPoint);
	opt.setConvergenceTolerance(1e-4);	//tol
	opt.setMaxIterations(200);	//200

	Vector result(_n_controls);
	result.setToZero();

	SimTK::Real f = 0.0;

//	qp.testQP();

	try{
		f = opt.optimize(result);
	}
	catch(const std::exception& ex)
	{
		std::cout<<ex.what()<<std::endl;
	}

	std::cout<<"Initial State Optimization Error: "<<f<<std::endl;

	SimTK::State stateOpt = _initState;
	
	Vector muscleForces(_n_muscles);
	Vector result_forces = result.elementwiseMultiply(_opt_actuator_forces);

	muscleForces = result_forces.block(0,0,_n_muscles,1).getAsVector();
	Vector tol_lm(_n_muscles), tol_acts(_n_muscles);
	Vector muscleFiberLens(_n_muscles), activations(_n_muscles);
	
	tol_lm.setTo(1.0e-6);
	tol_acts.setTo(1.0e-3);

	PrintVector(result,"result",std::cout);

	rootSolveMuscleFiberLength(stateOpt,muscleForces,_lm_min,_lm_max,tol_lm,muscleFiberLens);

	ROCINForceController::setStateFiberLength(*_model,muscleFiberLens,stateOpt);

	Vector lm_dot(_n_muscles);
	lm_dot.setToZero();

	_model->getMultibodySystem().realize(_initState,SimTK::Stage::Dynamics);
	

	int idx_musc = 0;
	for(int i=0;i<_n_controls;i++)
	{
		Actuator& act = _model->getActuators().get(i);
		Muscle* m = dynamic_cast<Muscle*>(&act);
		if(m!= NULL)
		{
			lm_dot[idx_musc] = m->getSpeed(_initState);//m->getSpeed(_initState)/m->getCosPennationAngle(_initState);//0;//m->getSpeed(_initState);
			idx_musc++;
		}
		
	}



	rootSolveActivations(stateOpt,muscleForces,muscleFiberLens,lm_dot,_acts_min,_acts_max,tol_acts,activations);

	ROCINForceController::setStateActivation(*_model,activations,stateOpt);	

	return stateOpt;
}

int StateInitializerQP::objectiveFunc(const Vector& coefficients, bool new_coefficients, SimTK::Real& f) const
{	
	f = 0.0;
	//penalize actuator foces
	f += coefficients.elementwiseMultiply(coefficients).elementwiseMultiply(_stateInitializer->_w_actuator_forces).sum();
	std::cout<<"objective function f: "<<f<<std::endl;
	return 0;
}

int StateInitializerQP::gradientFunc( const Vector& coefficients, bool new_coefficients, Vector &gradient ) const
{
	gradient.setToZero();

	//gradient for acutator forces
	gradient += coefficients.elementwiseMultiply(_stateInitializer->_w_actuator_forces)*2.0;
	

	return 0;
}

int StateInitializerQP::constraintFunc( const Vector& coefficients, bool new_coefficients, Vector &constraints) const
{
	constraints = _stateInitializer->_A_f*coefficients+_stateInitializer->_B_f-_stateInitializer->_UDotRef;

	return 0;
}

int StateInitializerQP::constraintJacobian( const Vector& coefficients, bool new_coefficients, Matrix& jac)  const
{
	jac = _stateInitializer->_A_f;
	return 0;
}

void StateInitializer::rootSolveActivations(const SimTK::State& s, const Vector& muscleForces, const Vector& lm, const Vector& lm_dot, const Vector& acts_min, const Vector& acts_max, const Vector& tol_acts, Vector& activations) const
{
	int numMuscles = lm.size();
	Array<double> para(0.0,numMuscles*2);
	for(int i=0;i<numMuscles;i++)
	{
		para[i] = lm[i];
		para[numMuscles+i] = lm_dot[i];
	}

	ROCINForceController::rootSolve(&(ROCINForceController::evaluateMuscleForceBasedOnActivation),*_model,s,muscleForces,acts_min,acts_max,tol_acts,activations,&para[0]);
}
void StateInitializer::rootSolveMuscleFiberLength(const SimTK::State& s, const Vector& muscleForces, const Vector& lm_min, const Vector& lm_max, const Vector& tol, Vector& muscleFiberLengthVec) const
{
	ROCINForceController::rootSolve(&(ROCINForceController::evaluateMuscleForceBasedOnFiberLength),*_model,s,muscleForces,lm_min,lm_max,tol,muscleFiberLengthVec);
}

void StateInitializerQP::testQP()
{
	int n_vars = _stateInitializer->_n_controls;


	Vector u(n_vars);

	for(int i=0;i<n_vars;i++)
		u.setTo(double(i)*1.47);


	double delta = 0.000001;
	Vector du(n_vars);

	Vector gradient_analytic(n_vars);
	Vector gradient_numeric(n_vars);

	//numeric gradient
	double f=0.0,f_new =0.0;
	objectiveFunc(u,true,f);
	gradientFunc(u,true,gradient_analytic);
	
	for(int i=0;i<n_vars;i++)
	{
		du.setToZero();
		du(i) = delta;
		objectiveFunc(u+du,true,f_new);
		double df = f_new-f;
		gradient_numeric(i) = df/delta; 
	}


	Matrix gradient_comparison(n_vars,3);
	gradient_comparison.updCol(0) = gradient_analytic;
	gradient_comparison.updCol(1) = gradient_numeric;
	gradient_comparison.updCol(2) = (gradient_analytic-gradient_numeric).elementwiseDivide(gradient_analytic);


	PrintMatrix(gradient_comparison,"gradient_comparison",std::cout);

	int n_c = _stateInitializer->_UDotRef.size();

	Matrix Jacob_analytic(n_c,n_vars);
	Matrix Jacob_numeric(n_c,n_vars);

	Vector c(n_c),c_new(n_c);
	constraintFunc(u,true,c);
	constraintJacobian(u,true,Jacob_analytic);

	for(int i=0;i<n_vars;i++)
	{
		du.setToZero();
		du(i) = delta;
		constraintFunc(u+du,true,c_new);
		Jacob_numeric.updCol(i) = (c_new-c)/delta;
	}

	Matrix Jacob_diff = Jacob_analytic - Jacob_numeric;
	Matrix Jacob_diff_normalize = Jacob_diff.elementwiseDivide(Jacob_analytic);

	PrintMatrix(Jacob_diff,"Jacob_diff",std::cout);
	PrintMatrix(Jacob_diff_normalize,"Jacob_diff_normalize",std::cout);

}

void StateInitializer::setMBSDynamicsAndOptForces(const Matrix& A_mbs, const Vector& B_mbs, const Vector& Opt_f)
{
	_opt_actuator_forces = Opt_f;
	Matrix Opt_f_diag(Opt_f.size(),Opt_f.size());
	Opt_f_diag.setToZero();
	Opt_f_diag.updDiag() = Opt_f;

	_A_f = A_mbs*Opt_f_diag;
	_B_f = B_mbs;
}