#include "ControllerHelper.h"
#include <OpenSim/Simulation/Model/Muscle.h>
#include "TestHelper.h"

using SimTK::Matrix;
using SimTK::Vector;
using SimTK::Vec3;
using namespace OpenSim;

void ROCINActuatorSubsystemRep::setCompleteState(const SimTK::State& state) {
	_completeState = state;
}
const SimTK::State& ROCINActuatorSubsystemRep::getCompleteState() const{
	return(_completeState);
}


FunctionSet* ROCINActuatorSubsystemRep::getCoordinateTrajectories() const {
   return( _qSet);
}

FunctionSet* ROCINActuatorSubsystemRep::getSpeedTrajectories() const {
	return ( _uSet);
}

void ROCINActuatorSubsystemRep::holdCoordinatesConstant(double t) {
	_holdCoordinatesConstant = true;
	_holdTime = t;
}

void ROCINActuatorSubsystemRep::releaseCoordinates() {
	_holdCoordinatesConstant = false;
}


//void ROCINActuatorSubsystemRep::setDesiredCoords(const double* desiredCoords)
//{	
//	int size = _qDesired.getSize();
//	for(int i=0;i<size;i++)
//		_qDesired[i] = desiredCoords[i];
//}

//void ROCINActuatorSubsystemRep::setDesiredSpeeds(const double* desiredSpeeds)
//{
//	int size = _uDesired.getSize();
//	for(int i=0;i<size;i++)
//		_uDesired[i] = desiredSpeeds[i];
//}

void ROCINActuatorSubsystemRep::setCoordinateCorrections(const double* aCorrections) {
    int i ;
    int size = _qCorrections.getSize();
    for(i=0;i<size;i++) {
         _qCorrections[i] = aCorrections[i];
    }
  
}

void ROCINActuatorSubsystemRep::setSpeedCorrections(const double* aCorrections) {
    int i ;
    int size = _uCorrections.getSize();
    for(i=0;i<size;i++) {
         _uCorrections[i] = aCorrections[i];
    }
}
  
void ROCINActuatorSubsystemRep::setCoordinateTrajectories(FunctionSet *aSet) {
	// ERROR CHECKING
	if(aSet == NULL) {
		std::string msg = "ROCINActuatorSubsystemRep.setCoordinateTrajectories:";
		msg += " ERR- NULL function set.\n";
		throw( Exception(msg,__FILE__,__LINE__) );
	}
	if(aSet->getSize() != _model->getNumCoordinates()) {
		std::string msg = "ROCINActuatorSubsystemRep.setCoordinateTrajectories:";
		msg += " ERR- incorrect number of trajectories.\n";
		throw( Exception(msg,__FILE__,__LINE__) );
	}

	_qSet = aSet;
}

void ROCINActuatorSubsystemRep::setSpeedTrajectories(FunctionSet *aSet) {
	// ERROR CHECKING
	if(aSet == NULL) {
		std::string msg = "ROCINActuatorSubsystemRep.setSpeedTrajectories:";
		msg += " ERR- NULL function set.\n";
		throw( Exception(msg,__FILE__,__LINE__) );
	}
	if(aSet->getSize() != _model->getNumSpeeds()) {
		std::string msg = "ROCINActuatorSubsystemRep.setSpeedTrajectories:";
		msg += " ERR- incorrect number of trajectories.\n";
		throw( Exception(msg,__FILE__,__LINE__) );
	}

	_uSet = aSet;
}


ROCINActuatorSubsystemRep::ROCINActuatorSubsystemRep(Model* model) 
       : SimTK::Subsystem::Guts( "ROCINActuatorSubsystem", "2.0"),
	   _holdCoordinatesConstant(false),
	   _holdTime(0.0),
	   _qSet(NULL),
       _model(model) {

       _qCorrections.setSize(_model->getNumCoordinates() );
       _uCorrections.setSize(_model->getNumSpeeds() );
       _qWork.setSize(_model->getNumCoordinates());
       _uWork.setSize(_model->getNumSpeeds());

	   for(int i=0;i<_model->getNumCoordinates();i++)
		   _qCorrections[i] = 0.0;

	   for(int i=0;i<_model->getNumSpeeds();i++)
		   _uCorrections[i] = 0.0;

//	   _qDesired.setSize(_model->getNumCoordinates());
//	   _uDesired.setSize(_model->getNumSpeeds());

	   _targetTime_transit = -SimTK::Infinity;
	   _interval_transit = 1.0;
	   
   }

ROCINActuatorSubsystemRep* ROCINActuatorSubsystemRep::cloneImpl() const { return new ROCINActuatorSubsystemRep(*this); }

ROCINActuatorSubsystemRep::~ROCINActuatorSubsystemRep() { }

int ROCINActuatorSubsystemRep::realizeSubsystemTopologyImpl(SimTK::State& s) const {
     // Make sure the CMC Actuator subsystem has the same number of Z's as
     // the model as a whole so we don't miss any of them. This will include
     // some that aren't muscle states, but who cares?

     Vector z(_model->getMultibodySystem().getDefaultState().getNZ());
     s.allocateZ( getMySubsystemIndex(), z );






     return 0;
  }

Model*  ROCINActuatorSubsystemRep::getModel() const {
       return( _model);
  }

void ROCINActuatorSubsystemRep::getEstimatedUDot(double t, SimTK::Vector& uDot) const
{
	int nu = _model->getNumSpeeds();
	uDot.resize(nu);	

    double alpha = 0.0;
	if(_interval_transit == 0.0)
		alpha = 1.0;
	else
	{
		alpha = (_targetTime_transit-t)/_interval_transit;

		if(alpha<0.0)
			alpha = 0.0;
		else if(alpha>1.0)
			alpha = 1.0;
	}

	if(alpha>1e-6)
	{
		Array<double> uDotTransit;

		if(_uSet != NULL)
		{
			_uSet->evaluate(_uWork,0,t);
			_uSet->evaluate(uDotTransit,0,_targetTime_transit);
		}
		else
		{
			_qSet->evaluate(_uWork,1,t);
			_qSet->evaluate(uDotTransit,1,_targetTime_transit);
		}

		double t_transit = _targetTime_transit - t;

		for(int i=0;i<nu;i++)
			uDot[i] = (uDotTransit[i]-(_uWork[i] + _uCorrections[i]*alpha))/t_transit;
	}
	else
	{
		Array<double> uDotWork;
		if(_uSet != NULL)
			_uSet->evaluate(uDotWork,1,t);
		else
			_qSet->evaluate(uDotWork,2,t);

		for(int i=0;i<nu;i++)
			uDot[i] = uDotWork[i];
	}

}

void ROCINActuatorSubsystemRep::assignQU(double t, Vector& q, Vector& u) const
{
    int nq = _model->getNumCoordinates();
    int nu = _model->getNumSpeeds();

    _qSet->evaluate(_qWork,0,t);

	if(_uSet != NULL)
		_uSet->evaluate(_uWork,0,t);
	else
		_qSet->evaluate(_uWork,1,t);


    double alpha = 0;
	if(_interval_transit == 0.0)
		alpha = 1.0;
	else
	{
		alpha = (_targetTime_transit-t)/_interval_transit;

		if(alpha<0.0)
			alpha = 0.0;
		else if(alpha>1.0)
			alpha = 1.0;
	}


    for(int i=0;i<nq;i++) q[i] = _qWork[i] + _qCorrections[i]*alpha;
    for(int i=0;i<nu;i++) u[i] = _uWork[i] + _uCorrections[i]*alpha;

    //for(i=0;i<nq;i++) q[i] = _qWork[i];
    //for(i=0;i<nu;i++) u[i] = _uWork[i]+_uCorrections[i];

    //for(int i=0;i<nq;i++) q[i] = _qWork[i];
    //for(int i=0;i<nu;i++) u[i] = _uWork[i];

/*
	for(int i=0;i<nq;i++) q[i] = _qCorrections[i]*alpha+_qDesired[i];
	for(int i=0;i<nu;i++) u[i] = _uCorrections[i]*alpha+_uDesired[i];*/
	

	//for(int i=0;i<nq;i++) q[i] = _qDesired[i];
	//for(int i=0;i<nu;i++) u[i] = _uDesired[i];


}

void ROCINActuatorSubsystemRep::assignTimeAndQU(double t, SimTK::State& s) const
{
	s.updTime() = t;
	assignQU(t,s.updQ(),s.updU());

}

int ROCINActuatorSubsystemRep::realizeSubsystemDynamicsImpl(const SimTK::State& s) const  {

	Vector& q = _model->getMultibodySystem().getMatterSubsystem().updQ(const_cast<SimTK::State&>(_completeState));
	Vector& u = _model->getMultibodySystem().getMatterSubsystem().updU(const_cast<SimTK::State&>(_completeState));

	double t;
	if(_holdCoordinatesConstant)
		t = _holdTime;
	else	
		t = s.getTime();

	assignQU(t,q,u);
	

    const_cast<SimTK::State&>(_completeState).updZ() = s.getZ();
	const_cast<SimTK::State&>(_completeState).updTime() = t;

//	std::cout<<"completeState.q:"<<_completeState.getQ()<<std::endl;
//	std::cout<<"completeState.u:"<<_completeState.getU()<<std::endl;
//	std::cout<<"completeState.z:"<<_completeState.getZ()<<std::endl;


//    std::cout << "_completeState before realization" << _completeState << std::endl;
     _model->getMultibodySystem().realize(_completeState, SimTK::Stage::Acceleration);
//     SimTK::Vector myControls = _model->getControls(_completeState);
//     PrintVector(myControls, "myControls",std::cout);
     //std::cout << "_completeState after realization" << _completeState << std::endl;
//     exit(0);

//	Vector controls_instant = _model->getControls(_completeState);
//	std::cout<<"controls_instant at t="<<t<<controls_instant<<std::endl;

//	 std::cout<<"completeState.zDot:"<<_completeState.getZDot()<<std::endl;
     
	s.updZDot() = _completeState.getZDot();
	 

 //   std::cout << "s.t = " << s.getTime() << std::endl;
 //   PrintVector(s.getZDot(), "s.ZDot", std::cout);
 //   exit(0);


	//std::cout << "_qWork=" << _qWork <<std::endl;
    //std::cout << "_uWork=" << _uWork <<std::endl;
	//std::cout << "q=" << q <<std::endl;
    //std::cout << "u=" << u <<std::endl;
	//std::cout << "actuatorStates=" << s.getZ() <<std::endl;
	//std::cout << " CMCrealize:Dynamics  time=" << s.getTime(); 
    //std::cout << " Actuator dydt=" << _completeState.getZDot() << std::endl;


     return 0;
  }      

ROCINActuatorSubsystem::ROCINActuatorSubsystem(ROCINActuatorSystem& system, Model* model) {
    adoptSubsystemGuts( rep = new ROCINActuatorSubsystemRep(model) );
    system.adoptSubsystem( *this );
}




ConstantController::ConstantController(int N)
{
	//_controls.resize(N);
	//_controls.setToZero();
	for(int i=0;i<N;i++)
	{
		//ControlLinear* control = new ControlLinear();
		//control->setUseSteps(true);
		ControlConstant* control = new ControlConstant(0.0);
		control->setIsModelControl(true);
		_controlSet.adoptAndAppend(control);
	}

}

void ConstantController::setControls(const double* controls, double t)
{
	//_controls = controls;
	//std::cout<<"setControls: size("<<_controlSet.getSize()<<"): ";

	//for(int i=0;i<_controlSet.getSize();i++)
	//	std::cout<<controls[i]<<" ";
	//std::cout<<std::endl;

	_controlSet.setControlValues(t,controls);

}

void ConstantController::computeControls(const SimTK::State& s, SimTK::Vector& controls) const
{

	for(int iActuator = 0; iActuator<getActuatorSet().getSize(); iActuator++)
	{
			
		const Actuator& thisActuator = getActuatorSet().get(iActuator);
		//Vector thisActuatorsControls(1, _controls[iActuator]);
		Vector thisActuatorsControls(1,_controlSet[iActuator].getControlValue(s.getTime()));
		thisActuator.addInControls(thisActuatorsControls,controls);		
	}

//	PrintVector(controls,"internal_controls",std::cout);
//    exit(0);

//	controls = _controls;

	return;
}

