#include "ROCINForceController.h"
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/Model/ForceSet.h>
#include <OpenSim/Simulation/SimbodyEngine/Body.h>
#include <OpenSim/Simulation/Model/Marker.h>
#include <OpenSim/Simulation/Model/MarkerSet.h>
#include <OpenSim/Simulation/Model/BodySet.h>
#include "SimTKmath.h"
#include <SimTKcommon/internal/Serialize.h>
#include <OpenSim/Common/GCVSplineSet.h>
#include <OpenSim/Common/Signal.h>
#include <OpenSim/Simulation/SimbodyEngine/CoordinateCouplerConstraint.h>
#include <OpenSim/Simulation/Model/HuntCrossleyForce.h>
#include <OpenSim/Simulation/Model/ContactGeometrySet.h>
#include <OpenSim/Simulation/Model/ContactSphere.h>
#include <OpenSim/Simulation/Model/AbstractTool.h>
#include <OpenSim/Common/IO.h>
#include <OpenSim/Simulation/Model/Force.h>
#include <OpenSim/Actuators/CoordinateActuator.h>
#include <OpenSim/Simulation/Model/CoordinateSet.h>
#include <OpenSim/Common/RootSolver.h>
#include <OpenSim/Simulation/Model/ControllerSet.h>
#include <OpenSim/Simulation/Control/Controller.h>

#include <OpenSim/Simulation/Control/ControlConstant.h>
#include <OpenSim/Simulation/Control/ControlLinear.h>

#include <OpenSim/Common/FunctionSet.h>
#include <OpenSim/Simulation/osimSimulationDLL.h>
#include "SimTKsimbody.h"
#include "SimTKcommon/internal/SubsystemGuts.h"
#include "SimTKcommon/internal/SystemGuts.h"

#include <OpenSim/Simulation/Model/ActivationFiberLengthMuscle.h>
#include <OpenSim/Actuators/Thelen2003Muscle.h>

#include "TestHelper.h"
#include "ControllerHelper.h"
#include "VirtualActuatorSet.h"

#include "StateInitializer.h"
#include <time.h>

using SimTK::Matrix;
using SimTK::Vector;
using SimTK::Vec3;

using namespace OpenSim;

#define VERBOSE_PRINT 0

class ROCINControlEventHandler : public SimTK::PeriodicEventHandler {
public:
	ROCINControlEventHandler(ROCINForceController* controller) :
	  PeriodicEventHandler(SimTK::Stage::Time), _controller(controller) {}

    // the event is calling doControlPrecomputation periodically
	void handleEvent (SimTK::State& s, SimTK::Real accuracy, bool& terminate) const {
		terminate = false;
		_controller->doControlPrecomputation(s);
		_controller->setTargetTime(s.getTime()+_controller->getTargetDT());
	}

	SimTK::Real getNextEventTime(const SimTK::State& s, bool includeCurrent) const
	{
		if(_controller->getCheckTargetTime())
			return _controller->getTargetTime();
		else
			return std::numeric_limits<SimTK::Real>::infinity();
	}

	ROCINForceController* _controller;

};


void ROCINForceController::initController(const Model& aModel, const std::string& dataFileName, bool track_coords, double lookahead_window, int lookahead_number)
{
	_track_coords = track_coords;
	_lookahead_window = lookahead_window;
	_lookahead_number = lookahead_number;
	setNull();


	_numNonMuscleActuators = 0;
	int n_u = aModel.getNumSpeeds();
	_mapping_nonMuscleActuators.resize(n_u,_numNonMuscleActuators);
	_mapping_nonMuscleActuators.setToZero();

	_checkTargetTime = false;

	if(_track_coords)
	{
		_coordData = new Storage(dataFileName);

		//TODO: make numCoord auto select only level q;

		if(aModel.getNumSpeeds() != aModel.getNumCoordinates())
		{
			std::cout<<"The number of speeds and the number of coordinates are not equal, does not support this case at this point!!"<<std::endl;
			exit(1);
		}

		_numCoords = aModel.getNumSpeeds();

	
	}
	else
	{
        std::cout << "Current implementation does not support tracking marker positions!" << std::endl;
        exit(1);
	}

}




ROCINForceController::ROCINForceController(const Model& aModel, const std::string& dataFileName, bool track_coords,  double lookahead_window, int lookahead_number)
{
	initController(aModel,dataFileName,track_coords,lookahead_window,lookahead_number);
}

void ROCINForceController::doControlPrecomputation(const SimTK::State& s)
{
    doPrecomputation_Varying_Dynamics(s);
}

ROCINForceController::~ROCINForceController()
{
	if(_mpcSolver != NULL)
		delete _mpcSolver;

	if(_markerData != NULL)
		delete _markerData;

	if(_coordData != NULL)
		delete _coordData;

	if(_externalLoads != NULL)
		delete _externalLoads;

	if(_momentArmSolver != NULL)
		delete _momentArmSolver;

	if(_qSet != NULL)
		delete _qSet;

	if(_uSet != NULL)
		delete _uSet;

	if(_uDotSet != NULL)
		delete _uDotSet;

	if(_internalIntegrator != NULL)
		delete _internalIntegrator;

	if(_internalManager != NULL)
		delete _internalManager;

//	if(_internalController != NULL)
//		delete _internalController;

	if(_internalModel != NULL)
		delete _internalModel;

}



void ROCINForceController::setNull()
{
	_n_q=0;
	_n_u=0;
	_numMuscles=0;
	_numMarkers=0;

	_coord_tracking_penalty = 50.0;
	_speed_tracking_penalty = 50.0;
	_activation_penalty = 1.0;
	_PD_penalty = 0.01;
	_lowpass_freq = -1.0;

	_use_taylor_expansion = false;

	_markerData = NULL;
	_mpcSolver = NULL;

	_coordData = NULL;
	_markerData = NULL;
	_qSet = NULL;
	_uSet = NULL;
	_uDotSet = NULL;

	_externalLoads = NULL;
	_momentArmSolver = NULL;
	_internalModel = NULL;
	_internalIntegrator = NULL;
	_internalManager = NULL;
	_internalController = NULL;
	
	_internalSys = NULL;
	_internalSubsys = NULL;

}


void ROCINForceController::evalModelCoordTrackingErrs(const Model& aModel, const std::string& coordFileName1, const std::string& coordFileName2, const std::string& outFileName)
{
	Storage* refCoordData = new Storage(coordFileName1);
	Storage* cmpCoordData = new Storage(coordFileName2);

	if(refCoordData->isInDegrees())
		aModel.getSimbodyEngine().convertDegreesToRadians(*refCoordData);

	if(cmpCoordData->isInDegrees())
		aModel.getSimbodyEngine().convertDegreesToRadians(*cmpCoordData);

	std::ofstream fout(outFileName);
	if(!fout.is_open())
	{
		std::cout<<"Fail to open file "<<outFileName.c_str()<<" to write!!!"<<std::endl;
		exit(1);
	}

	int n_samples_cmp = cmpCoordData->getSize();
	int numCoords_ref = refCoordData->getColumnLabels().size()-1;

    //write the header
	fout<<"CoordErrors"<<std::endl;
	fout<<"version=1"<<std::endl;
	fout<<"nRows="<<n_samples_cmp<<std::endl;
	fout<<"nColumns="<<1+numCoords_ref<<std::endl;
	fout<<"inDegrees=no"<<std::endl;
	fout<<"endheader"<<std::endl;

	const Array<std::string>& strs_labels = refCoordData->getColumnLabels();
	const Array<std::string>& strs_labels_cmp = cmpCoordData->getColumnLabels();

	int numCoords_cmp = strs_labels_cmp.size()-1;

	Array<int> corresponding_indices_in_cmp(-1,numCoords_ref);

	for(int i=1;i<=numCoords_ref;i++)
	{
		for(int j=1;j<=numCoords_cmp;j++)
		{
			if(strs_labels[i] == strs_labels_cmp[j])
			{
				corresponding_indices_in_cmp[i-1] = j-1;
				break;
			}
		}
	}
	
	fout<<"time";
	for(int i=1;i<=numCoords_ref;i++)
		fout<<" "<<strs_labels.get(i).c_str();
	fout<<std::endl;

	SimTK::Real t;
	SimTK::Vector q(numCoords_cmp);
	SimTK::Vector q_ref(numCoords_ref);
	Vector err_q(numCoords_ref);

	for(int i=0;i<n_samples_cmp;i++)
	{
		
		cmpCoordData->getTime(i,t);	
		cmpCoordData->getData(i,numCoords_cmp,q);
		refCoordData->getDataAtTime(t,numCoords_ref,q_ref);

		for(int j=0;j<numCoords_ref;j++)
		{
			if(corresponding_indices_in_cmp[j]<0)
				err_q[j] = 0;
			else
				err_q[j] = q[corresponding_indices_in_cmp[j]]-q_ref[j];
		}
		
		fout<<t;
		for(int j=0;j<numCoords_ref;j++)
			fout<<" "<<err_q[j];

		fout<<std::endl;
		
	}

	fout.close();

	delete refCoordData;	
	delete cmpCoordData;

}


void ROCINForceController::evalCoordTrackingErrs(const std::string& coordFileName1, const std::string& coordFileName2, const std::string& outFileName)
{
	evalModelCoordTrackingErrs(getModel(),coordFileName1,coordFileName2,outFileName);
}

void ROCINForceController::initMPCSolver(int numExtraObservations)
{
    // number of observations
	int n_y = _numCoords*2+numExtraObservations;

    // number of control variables
	int n_controls = _numMuscles+_numNonMuscleActuators;

	double initTime, finalTime;
	Array_<Vector> y_reference;

	const SimTK::State& initState = getModel().getWorkingState();

	int n_samples = 0;
	
	{

        // get the trajectories of the joint angles and joint velocities
		if(_coordData->isInDegrees())
			getModel().getSimbodyEngine().convertDegreesToRadians(*_coordData);

		initTime = _coordData->getFirstTime();
		finalTime = _coordData->getLastTime();

		n_samples = floor((finalTime-initTime)/_lookahead_window)+1;


		Storage _coordDataPad = *_coordData;
		_coordDataPad.pad(60);

		if(_lowpass_freq>0.0)
			_coordDataPad.lowpassFIR(50,_lowpass_freq);


        // form the complete storage
		Storage* qStore = NULL;
		Storage* uStore = NULL;

		getModel().getSimbodyEngine().formCompleteStorages(initState,_coordDataPad,qStore,uStore);

		if(qStore->isInDegrees())
		{
			getModel().getSimbodyEngine().convertDegreesToRadians(*qStore);
			getModel().getSimbodyEngine().convertDegreesToRadians(*uStore);
		}
		
		_qSet = new GCVSplineSet(5,qStore);		
		_uSet = new GCVSplineSet(5,uStore);
		Storage* dudtStore = _uSet->constructStorage(1);
		_uDotSet = new GCVSplineSet(5,dudtStore);


		//_qSet->print("test_qSet_ROCIN.xml");

		delete qStore;
		delete uStore;
		delete dudtStore;

		FunctionSet * yRefSet = new FunctionSet();
		for(int i=0;i<_qSet->getSize();i++)
			yRefSet->cloneAndAppend(_qSet->get(i));
		
		for(int i=0;i<_uSet->getSize();i++)
			yRefSet->cloneAndAppend(_uSet->get(i));

        // construct the MPC solver
		_mpcSolver = new MPC(n_controls,n_y,n_samples,initTime,_lookahead_window,yRefSet,_lookahead_number);
	}

    // initialize the variables
	_actuator_forces.setSize(n_controls);
	_actuator_controls.setSize(n_controls);
	_actuator_controls_new.setSize(n_controls);
	for(int i=0;i<n_controls;i++)
	{
		_actuator_controls.set(i,0.0);	//set to be 0.01 
		_actuator_controls_new.set(i,0.0);
	}

	_gen_forces_array.setSize(_lookahead_number);
	_actuator_forces_array.setSize(_lookahead_number);

	for(int i=0;i<_lookahead_number;i++)
	{
		_actuator_forces_array[i].resize(n_controls);
		_actuator_forces_array[i].setToZero();
		_gen_forces_array[i].resize(_n_u);
		_gen_forces_array[i].setToZero();
	}

	_muscle_min_control_array.resize(_numMuscles);
	_muscle_max_control_array.resize(_numMuscles);

	int idx_musc=0;
	for(int i=0;i<_numMuscles+_numNonMuscleActuators;i++)
	{
		Actuator& act = getModel().getActuators().get(i);
		PathActuator* m = dynamic_cast<PathActuator*>(&act);
		if(m!=NULL)
		{

			_muscle_min_control_array[idx_musc] = m->getMinControl();
			_muscle_max_control_array[idx_musc] = m->getMaxControl();

			idx_musc++;	
		}

	}

	setTargetDT(_lookahead_window);
	

    // set the lower and upper bounds for the control variables
	_control_lowerbounds.resize(n_controls*_lookahead_number);
	_control_upperbounds.resize(n_controls*_lookahead_number);

	if(_numMuscles>0)
	{
		for(int i=0;i<_lookahead_number;i++)
		{
			_control_lowerbounds.updBlock(i*n_controls,0,_numMuscles,1) = _muscle_min_control_array;
			_control_upperbounds.updBlock(i*n_controls,0,_numMuscles,1) = _muscle_max_control_array;
		}
	}
	
	if(_numNonMuscleActuators>0)
	{
		for(int i=0;i<_lookahead_number;i++)
		{
			_control_lowerbounds.updBlock(i*n_controls+_numMuscles,0,_numNonMuscleActuators,1) = _virtual_actuator_min_control_array;		//200 for running	//-40 for standing
			_control_upperbounds.updBlock(i*n_controls+_numMuscles,0,_numNonMuscleActuators,1) = _virtual_actuator_max_control_array;			//30
		}
	}

	_mpcSolver->setLowerBounds(_control_lowerbounds);
	_mpcSolver->setUpperBounds(_control_upperbounds);


    // initialize the matrix coefficients for the dynamics equation: y_dot = A*y + B*u+ C (y = x in our case)
    // these matrix values are not really used since we will assign them again in the doControlPrecomputation
	Matrix A(n_y,n_y);
	Matrix B(n_y,n_controls);
	Vector C(n_y);

	A.setToZero();
	A.diag().setTo(1.0);

	B.setToZero();
	B.diag().setTo(1.0);

	C.setToZero();
	//C.setTo(0.3);

	Array<Matrix> A_array;
	Array<Matrix> B_array;
	Array<Vector> C_array;

	A_array.setSize(_lookahead_number);
	B_array.setSize(_lookahead_number);
	C_array.setSize(_lookahead_number);

	for(int i=0;i<_lookahead_number;i++)
	{
		A_array[i] = A;
		B_array[i] = B;
		C_array[i] = C;
	}

	_mpcSolver->setABCArray(A_array,B_array,C_array);
	_mpcSolver->setABC(A,B,C);

    // D and E, one item in the objective function is (D*u+E)^T*R*(D*u+E)
	Matrix D(n_controls,n_controls);
	Vector E(n_controls);

	D.setToZero();
	D.diag().setTo(1.0);
	E.setToZero();	

	Array<Vector> D_array;
	Array<Vector> E_array;
	D_array.setSize(_lookahead_number);
	E_array.setSize(_lookahead_number);

	for(int i=0;i<_lookahead_number;i++)
	{
		D_array[i] = D.diag();
		E_array[i] = E;
	}

	_mpcSolver->setDiagDandEAarray(D_array,E_array);
	_mpcSolver->setDiagDandE(D.diag(),E);

    // currently we don't penalize the terminal time error
	Matrix P(n_y,n_y);
	P.setToZero();

    // penalty matrix, objective function item: (y-\hat y)^T * Q* (y- \hat y)
	Matrix Q(n_y,n_y);
	Q.setToZero();

	int n_tracking = n_y - numExtraObservations;
	//Q.updBlock(0,0,n_tracking,n_tracking).diag().setTo(_tracking_penalty);	//5.0e1	
	Q.updBlock(0,0,n_tracking/2,n_tracking/2).diag().setTo(_coord_tracking_penalty);	//5.0e1	
	Q.updBlock(n_tracking/2,n_tracking/2,n_tracking/2,n_tracking/2).diag().setTo(_speed_tracking_penalty);
	
	Matrix R(n_controls,n_controls);
	R.setToZero();
	
	R.updBlock(0,0,_numMuscles,_numMuscles).diag().setTo(_activation_penalty);

	if(_numNonMuscleActuators>0)
	{
		for(int i=0;i<_numNonMuscleActuators;i++)
			R.set(_numMuscles+i,_numMuscles+i,_virtual_actuator_weight_array[i]);		
	}

    // objective function time: (\dot y - \hat\dot y)^T * Qd * (\dot y - \hat\dot y)

	Matrix Qd(n_y,n_y);
	Qd.setToZero();

    Qd.updBlock(n_y/2,n_y/2,n_y/2,n_y/2).diag().setTo(_PD_penalty);		//1.0e-3

    // objective function time: \dot u^T * Sd * \dot u
	Matrix Sd(n_controls,n_controls);
	Sd.setToZero();

    // do not track the coordinates or speeds that are constrained
	{
		const Array<std::string>& labelArray = _coordData->getColumnLabels();

		_CoordsSelection.resize(_numCoords,_n_u);
		_CoordsSelection.setToZero();			

		for(int i=0;i<_numCoords;i++)
		{
			int idx_coord = i;
						
			//if(getModel().getCoordinateSet()[idx_coord].getDefaultLocked() || getModel().getCoordinateSet()[idx_coord].getDefaultIsPrescribed())
			if(getModel().getCoordinateSet()[idx_coord].isConstrained(initState))
			{
				P.set(i,i,0.0);
				P.set(_numCoords+i,_numCoords+i,0.0);
				Q.set(i,i,0.0);
				Q.set(_numCoords+i,_numCoords+i,0.0);
				Qd.set(i,i,0.0);
				Qd.set(_numCoords+i,_numCoords+i,0.0);
			}
			else
			{
				_CoordsSelection.set(i,idx_coord,1.0);
			}
		}

	}


	if(_numNonMuscleActuators>0)
	{
		for(int i=0;i<_numNonMuscleActuators;i++)
		{
			for(int j=0;j<_n_u;j++)
			{
				if(_mapping_nonMuscleActuators.get(j,i)>0)
				{
					if(getModel().getCoordinateSet().get(j).isConstrained(initState))
					//if(getModel().getCoordinateSet()[j].getDefaultLocked() || getModel().getCoordinateSet()[j].getDefaultIsPrescribed())
					{						

						R.set(_numMuscles+i,_numMuscles+i,0.0);
						Sd.set(_numMuscles+i,_numMuscles+i,0.0);
						_mapping_nonMuscleActuators.set(j,i,0.0);
					}
					break;
				}
			}
		}
	}

	_mpcSolver->setDiagPQR(P.diag(),Q.diag(),R.diag());

	_mpcSolver->setDiagQd(Qd.diag());
	_mpcSolver->setDiagSd(Sd.diag());
}

void ROCINForceController::connectToModel(Model& model)
{

//	setUpInternalModel(model);

	Super::connectToModel(model);

	setNumControls(getActuatorSet().getSize());
	
	_n_u = getModel().getNumSpeeds();
	_n_q = getModel().getNumCoordinates();


	// check coupled coordinates


	const ForceSet& forceSet = getModel().getForceSet();
	_numMuscles = 0;
	for(int i=0;i<forceSet.getSize();i++)
	{
		PathActuator* m = dynamic_cast<PathActuator*>(&forceSet.get(i));

		if(m!=NULL)
			_numMuscles++;
	}
	

}

void ROCINForceController::addEventHandlerToSystem()
{
	ROCINControlEventHandler* ROCINHandler = new ROCINControlEventHandler(this);
	getModel().updMultibodySystem().updDefaultSubsystem().addEventHandler(ROCINHandler);

}

void ROCINForceController::computeOptimalInitialState(SimTK::State& initState, double tiReal)
{
    // solve an optimization problem that minimize the joint acceleration error and muscle activations

	_internalSubsys->assignTimeAndQU(tiReal,initState);
	_internalModel->getMultibodySystem().realize(initState,SimTK::Stage::Velocity);

	if(initState.getZ().size()==0)
		return;


	StateInitializer stateInitializer(_internalModel,_internalController,initState);

	Array<double> UDotRef;
	if(_uSet != NULL)
		_uSet->evaluate(UDotRef,1,tiReal);
	else
		_qSet->evaluate(UDotRef,2,tiReal);


	Vector UDotRef_vec(_n_u);
	for(int i=0;i<_n_u;i++)
		UDotRef_vec[i] = UDotRef[i];

	int nActs = getActuatorSet().getSize();

	Vector w_actuator_forces(nActs), lower_bd_forces(nActs), upper_bd_forces(nActs);
	

	Vector maxIsoForces(_numMuscles);

	Vector acts_min(_numMuscles),acts_max(_numMuscles);
	Vector lm_min(_numMuscles), lm_max(_numMuscles);
	acts_min.setTo(0.02);
	acts_max.setTo(1.0);
	
	int idx_musc = 0;
	for(int i=0;i<nActs;i++)
	{
		Actuator& act = getActuatorSet().get(i);
		Muscle* musc = dynamic_cast<Muscle*>(&act);
		if(musc != NULL)
		{
			double h = musc->getOptimalFiberLength()*sin(musc->getPennationAngleAtOptimalFiberLength());
			double lmt = musc->getLength(initState);
			double lst = musc->getTendonSlackLength();
			lm_min[idx_musc] = h;
			lm_max[idx_musc] = sqrt(h*h+(lmt-lst)*(lmt-lst));


			maxIsoForces[idx_musc] = musc->getMaxIsometricForce();
			w_actuator_forces[idx_musc] = 1.0;//1.0/(musc->getMaxIsometricForce()*musc->getMaxIsometricForce());

			idx_musc++;
		}
	}

	for(int i=0;i<_numNonMuscleActuators;i++)
		w_actuator_forces[_numMuscles+i] = _virtual_actuator_weight_array[i];

	


	lower_bd_forces.updBlock(0,0,_numMuscles,1).setTo(0.0);
	Vector maxMuscleForces(_numMuscles);
	evaluateMuscleForceBasedOnFiberLength(*_internalModel,initState,NULL,lm_min,maxMuscleForces);
	upper_bd_forces.updBlock(0,0,_numMuscles,1) = maxMuscleForces.elementwiseDivide(maxIsoForces);

	lower_bd_forces.updBlock(_numMuscles,0,_numNonMuscleActuators,1) = _virtual_actuator_min_control_array;
	upper_bd_forces.updBlock(_numMuscles,0,_numNonMuscleActuators,1) = _virtual_actuator_max_control_array;

	Vector optForces(_numMuscles+_numNonMuscleActuators);
	optForces.updBlock(0,0,_numMuscles,1) = maxIsoForces;
	optForces.updBlock(_numMuscles,0,_numNonMuscleActuators,1).setTo(1.0);


	Matrix A_mbs;
	Vector B_mbs;
	computeMBSDynamicsNumerically(initState,A_mbs,B_mbs);

	stateInitializer.setUDotRef(UDotRef_vec);
	stateInitializer.setActuatorForceWeights(w_actuator_forces);
	stateInitializer.setAccWeights(1.0e3);
	stateInitializer.setActuatorForceBounds(lower_bd_forces,upper_bd_forces);
	stateInitializer.setActivationBounds(acts_min,acts_max);
	stateInitializer.setFiberLengthBounds(lm_min,lm_max);
	stateInitializer.setMBSDynamicsAndOptForces(A_mbs,B_mbs,optForces);

	initState = stateInitializer.getOptimizedInitState();	


	PrintVector(initState.getZ(),"initState.Z",std::cout);

}

SimTK::State& ROCINForceController::initSystem(Model& aModel)
{
	setCheckTargetTime(true);	
	aModel.buildSystem();
	addEventHandlerToSystem();
	SimTK::State& initState = aModel.initializeState();


	initMPCSolver(0);//_numMuscles



	//initialize control Set

	int nActs = getActuatorSet().getSize();
	for(int i=0;i<nActs;i++)
	{
		Actuator& act = getActuatorSet().get(i);

		ControlLinear* control = new ControlLinear();
		control->setName(act.getName()+".excitation");

		{
			control->setUseSteps(true);
		}

		_controlSet.adoptAndAppend(control);
	}

	_internalSubsys->setCoordinateTrajectories(_qSet);
	_internalSubsys->setSpeedTrajectories(_uSet);

	_momentArmSolver = new MomentArmSolver(aModel);

	return initState;
}

SimTK::State& ROCINForceController::setInitStateFromFile(Model& aModel, const std::string& aStateFileName)
{



	SimTK::State& initState = initSystem(aModel);
    Storage tmp(aStateFileName);
    Storage* stateStorage = new Storage();
    aModel.formStateStorage(tmp, *stateStorage);

	
	const Array<std::string>& strs_labels = stateStorage->getColumnLabels();

	int n_states = strs_labels.size()-1;

	Vector initStateVector(n_states);
	double initTime;
	stateStorage->getTime(0,initTime);
	stateStorage->getData(0,n_states,initStateVector);

	for(int i=1;i<=n_states;i++)
		aModel.setStateVariable(initState,strs_labels[i],initStateVector[i-1]);

	initState.updTime() = initTime;

    delete stateStorage;

	aModel.getMultibodySystem().realize(initState,SimTK::Stage::Velocity);

	return initState;
}

SimTK::State& ROCINForceController::setInitState(Model& aModel, double tiReal)
{
	SimTK::State& initState = initSystem(aModel);

	computeOptimalInitialState(initState,tiReal);

	return initState;
}

// using forward difference method to compute the partial derivatives
void ROCINForceController::computeGeneralizedForceDynamicsPartialsFwdDiff(const SimTK::State& s, const Vector& genForces, Matrix& dUDotdQ, Matrix& dUDotdU)
{
	dUDotdQ.resize(_n_u,_n_q);
	dUDotdU.resize(_n_u,_n_u);


	SimTK::State s_cpy = s;

	const SimTK::SimbodyMatterSubsystem& smss = getModel().getMatterSubsystem();		
	const SimTK::MultibodySystem& mbs = getModel().getMultibodySystem();

	const Set<Actuator>& fSet = getActuatorSet();
	int nf = fSet.getSize();

	for(int i=0;i<nf;i++)
	{

        ScalarActuator& act = (ScalarActuator&)fSet.get(i);
        act.overrideActuation(s_cpy, true);
        act.setOverrideActuation(s_cpy, 0.0);
	}

	mbs.realize(s_cpy,SimTK::Stage::Dynamics);

	const SimTK::Vector_<SimTK::SpatialVec>& bodyForces = mbs.getRigidBodyForces(s_cpy,SimTK::Stage::Dynamics);

	SimTK::Vector_<SimTK::SpatialVec> bodyAccs;

	Vector uDot_0(_n_u);

	smss.calcAcceleration(s_cpy,genForces,bodyForces,uDot_0,bodyAccs);

	Vector uDot(_n_u);

	const Vector& q_s = s.getQ();
	const Vector& u_s = s.getU();

	double delta = 1e-6;

	for(int i=0;i<_n_q;i++)
	{
		s_cpy.updQ()[i] = q_s[i]+delta;
		mbs.realize(s_cpy,SimTK::Stage::Dynamics);
		smss.calcAcceleration(s_cpy,genForces,bodyForces,uDot,bodyAccs);
		dUDotdQ.updCol(i) = (uDot-uDot_0)/delta;
		s_cpy.updQ()[i] = q_s[i];
	}

	for(int i=0;i<_n_u;i++)
	{
		s_cpy.updU()[i] = u_s[i]+delta;
		mbs.realize(s_cpy,SimTK::Stage::Dynamics);
		smss.calcAcceleration(s_cpy,genForces,bodyForces,uDot,bodyAccs);
		dUDotdU.updCol(i) = (uDot-uDot_0)/delta;
		s_cpy.updU()[i] = u_s[i];
	}

}

void ROCINForceController::computeDynamicsPartialsFwdDiff(const SimTK::State& s, const Vector& actuatorForces,Matrix& dUDotdQ, Matrix& dUDotdU)
{
	SimTK::State s_cpy = s;

	const SimTK::SimbodyMatterSubsystem& smss = getModel().getMatterSubsystem();		
	const SimTK::MultibodySystem& mbs = getModel().getMultibodySystem();
	

	const Set<Actuator>& fSet = getActuatorSet();
	int nf = fSet.getSize();

	for(int i=0;i<nf;i++)
	{

        ScalarActuator& act = (ScalarActuator&)fSet.get(i);
        act.overrideActuation(s_cpy, true);
        act.setOverrideActuation(s_cpy, actuatorForces[i]);
	}

	mbs.realize(s_cpy,SimTK::Stage::Acceleration);
	Vector udot_s = s_cpy.getUDot();

	

	const Vector& q_s = s.getQ();
	const Vector& u_s = s.getU();

	double delta = 1e-6;

	int n_q = q_s.size();
	int n_u = u_s.size();

	dUDotdQ.resize(n_u,n_q);
	dUDotdU.resize(n_u,n_u);

	for(int i=0;i<n_q;i++)
	{
		s_cpy.updQ()[i] = q_s[i]+delta;
		mbs.realize(s_cpy,SimTK::Stage::Acceleration);
		dUDotdQ.updCol(i) = (s_cpy.getUDot()-udot_s)/delta;
		s_cpy.updQ()[i] = q_s[i];
	}

	for(int i=0;i<n_u;i++)
	{
		s_cpy.updU()[i] = u_s[i]+delta;
		mbs.realize(s_cpy,SimTK::Stage::Acceleration);
		dUDotdU.updCol(i) = (s_cpy.getUDot()-udot_s)/delta;
		s_cpy.updU()[i] = u_s[i];
	}
	

}

// use central difference methods to compute the partial derivatives
void ROCINForceController::computeGeneralizedForceDynamicsPartialsCentralDiff(const SimTK::State& s, const Vector& genForces, Matrix& dUDotdQ, Matrix& dUDotdU)
{
	dUDotdQ.resize(_n_u,_n_q);
	dUDotdU.resize(_n_u,_n_u);

	SimTK::State s_cpy = s;

	const SimTK::SimbodyMatterSubsystem& smss = getModel().getMatterSubsystem();		
	const SimTK::MultibodySystem& mbs = getModel().getMultibodySystem();
	

	const Set<Actuator>& fSet = getActuatorSet();
	int nf = fSet.getSize();

	for(int i=0;i<nf;i++)
	{
        ScalarActuator& act = (ScalarActuator&)fSet.get(i);
        act.overrideActuation(s_cpy, true);
        act.setOverrideActuation(s_cpy, 0.0);
	}

	mbs.realize(s_cpy,SimTK::Stage::Dynamics);

	const SimTK::Vector_<SimTK::SpatialVec>& bodyForces = mbs.getRigidBodyForces(s_cpy,SimTK::Stage::Dynamics);

	SimTK::Vector_<SimTK::SpatialVec> bodyAccs;

	const Vector& q_s = s.getQ();
	const Vector& u_s = s.getU();

	double delta = 1e-5;

	Vector udot_plus(_n_u), udot_minus(_n_u);

	for(int i=0;i<_n_q;i++)
	{
		s_cpy.updQ()[i] = q_s[i]+delta;
		mbs.realize(s_cpy,SimTK::Stage::Dynamics);
		smss.calcAcceleration(s_cpy,genForces,bodyForces,udot_plus,bodyAccs);

		s_cpy.updQ()[i] = q_s[i]-delta;
		mbs.realize(s_cpy,SimTK::Stage::Dynamics);
		smss.calcAcceleration(s_cpy,genForces,bodyForces,udot_minus,bodyAccs);
		
		dUDotdQ.updCol(i) = (udot_plus-udot_minus)/(2.0*delta);
		s_cpy.updQ()[i] = q_s[i];
	}

	for(int i=0;i<_n_u;i++)
	{
		s_cpy.updU()[i] = u_s[i]+delta;
		mbs.realize(s_cpy,SimTK::Stage::Dynamics);
		smss.calcAcceleration(s_cpy,genForces,bodyForces,udot_plus,bodyAccs);
		

		s_cpy.updU()[i] = u_s[i]-delta;
		mbs.realize(s_cpy,SimTK::Stage::Dynamics);
		smss.calcAcceleration(s_cpy,genForces,bodyForces,udot_minus,bodyAccs);

		dUDotdU.updCol(i) = (udot_plus-udot_minus)/(2.0*delta);
		s_cpy.updU()[i] = u_s[i];
	}


}

void ROCINForceController::computeDynamicsPartialsCentralDiff(const SimTK::State& s, const Vector& actuatorForces,Matrix& dUDotdQ, Matrix& dUDotdU)
{
	SimTK::State s_cpy = s;

	const SimTK::SimbodyMatterSubsystem& smss = getModel().getMatterSubsystem();		
	const SimTK::MultibodySystem& mbs = getModel().getMultibodySystem();

	const Set<Actuator>& fSet = getActuatorSet();
	int nf = fSet.getSize();

	for(int i=0;i<nf;i++)
	{
        ScalarActuator& act = (ScalarActuator&)fSet.get(i);
        act.overrideActuation(s_cpy, true);
        act.setOverrideActuation(s_cpy, actuatorForces[i]);
	}
	
	
	const Vector& q_s = s.getQ();
	const Vector& u_s = s.getU();

	int n_q = q_s.size();
	int n_u = u_s.size();

	dUDotdQ.resize(n_u,n_q);
	dUDotdU.resize(n_u,n_u);

	double delta = 1e-5;

	Vector udot_plus(n_u), udot_minus(n_u);

	for(int i=0;i<n_q;i++)
	{
		s_cpy.updQ()[i] = q_s[i]+delta;
		mbs.realize(s_cpy,SimTK::Stage::Acceleration);
		udot_plus = s_cpy.getUDot();

		s_cpy.updQ()[i] = q_s[i]-delta;
		mbs.realize(s_cpy,SimTK::Stage::Acceleration);
		udot_minus = s_cpy.getUDot();
		
		dUDotdQ.updCol(i) = (udot_plus-udot_minus)/(2.0*delta);
		s_cpy.updQ()[i] = q_s[i];
	}

	for(int i=0;i<n_u;i++)
	{
		s_cpy.updU()[i] = u_s[i]+delta;
		mbs.realize(s_cpy,SimTK::Stage::Acceleration);
		udot_plus = s_cpy.getUDot();

		s_cpy.updU()[i] = u_s[i]-delta;
		mbs.realize(s_cpy,SimTK::Stage::Acceleration);
		udot_minus = s_cpy.getUDot();

		dUDotdU.updCol(i) = (udot_plus-udot_minus)/(2.0*delta);
		s_cpy.updU()[i] = u_s[i];
	}

}

// \dot U (i.e. \ddot q) = B*f_actuator + A 
void ROCINForceController::computeMBSDynamicsNumerically(const SimTK::State& s, Matrix& B, Vector& A)
{
	B.resize(_n_u,_numMuscles+_numNonMuscleActuators);
	A.resize(_n_u);

	SimTK::State s_cpy = s;

	const SimTK::SimbodyMatterSubsystem& smss = getModel().getMatterSubsystem();		
	const SimTK::MultibodySystem& mbs = getModel().getMultibodySystem();

	const Set<Actuator>& fSet = getActuatorSet();
	int nf = fSet.getSize();


	for(int i=0;i<nf;i++)
	{

        ScalarActuator& act = (ScalarActuator&)fSet.get(i);
        act.overrideActuation(s_cpy, true);
        act.setOverrideActuation(s_cpy, 0.0);
	}

	mbs.realize(s_cpy,SimTK::Stage::Acceleration);

	A = s_cpy.getUDot();

	for(int i=0;i<nf;i++)
	{
        ScalarActuator& act = (ScalarActuator&)fSet.get(i);
        act.setOverrideActuation(s_cpy, 1.0);



		mbs.realize(s_cpy,SimTK::Stage::Acceleration);
		B.updCol(i) = s_cpy.getUDot()-A;

        act.setOverrideActuation(s_cpy, 0.0);
	}

}

// \tau = C_muscles * fm
void ROCINForceController::computeMuscleJacobian(const SimTK::State& s, Matrix& C_muscles)
{
	const SimTK::SimbodyMatterSubsystem& smss = getModel().getMatterSubsystem();		
	const SimTK::MultibodySystem& mbs = getModel().getMultibodySystem();
	mbs.realize(s,SimTK::Stage::Velocity);

	const ForceSet& forceSet = getModel().getForceSet();


	//Muscle Jacobian:
	C_muscles.resize(_n_u, _numMuscles);	
	
	
	int col_idx=0;
	for(int i=0;i<forceSet.getSize();i++)
	{
		PathActuator* m = dynamic_cast<PathActuator*>(&forceSet.get(i));
		if(m!=NULL)
		{			
			//compute Muscle Jacobian using OpenSim API
			for(int j=0;j<_n_u;j++)
			{
				Coordinate& coord = getModel().getCoordinateSet().get(j);

				//check whether the coordinate is constrained
				if(coord.isDependent(s))
					C_muscles.set(j,col_idx,0.0);				
				else
					C_muscles.set(j,col_idx,_momentArmSolver->solve(s,coord,m->getGeometryPath()));

			}
								
			col_idx++;
			
		}

	}

}

void ROCINForceController::computeMBSDynamics(const SimTK::State& s, Matrix& B, Vector& A)
{
	//q_ddot = A+B*f; //f is the vector of actuator forces

	B.resize(_n_u,_numMuscles+_numNonMuscleActuators);
	A.resize(_n_u);

	double tiReal = s.getTime();

	const SimTK::SimbodyMatterSubsystem& smss = getModel().getMatterSubsystem();		
	const SimTK::MultibodySystem& mbs = getModel().getMultibodySystem();
	mbs.realize(s,SimTK::Stage::Velocity);


	//Inverse Mass Matrix 
	Matrix MInv(_n_u,_n_u);
	smss.calcMInv(s,MInv);


	Vector CminusE(_n_u);
	const SimTK::Vector_<SimTK::SpatialVec>& appliedGravityOnBody = getModel().getGravityForce().getBodyForces(s);
	SimTK::Vector_<SimTK::SpatialVec> totalExternalForces = appliedGravityOnBody;

	// Figure out number of muscles
	const ForceSet& forceSet = getModel().getForceSet();


	//Muscle Jacobian:
	Matrix C_muscles(_n_u, _numMuscles);	
	
	
	int col_idx=0;
	for(int i=0;i<forceSet.getSize();i++)
	{
		PathActuator* m = dynamic_cast<PathActuator*>(&forceSet.get(i));
		if(m!=NULL)
		{			
			//compute Muscle Jacobian using OpenSim API
			for(int j=0;j<_n_u;j++)
			{
				Coordinate& coord = getModel().getCoordinateSet().get(j);

				//check whether the coordinate is constrained
				if(coord.isDependent(s))
					C_muscles.set(j,col_idx,0.0);				
				else
					C_muscles.set(j,col_idx,_momentArmSolver->solve(s,coord,m->getGeometryPath()));

			}
								
			col_idx++;
			
		}
		else
		{
			ExternalForce* force_ext = dynamic_cast<ExternalForce*>(&forceSet.get(i));

			double time_force = tiReal;
			if(time_force > getFinalTime())
				time_force = getFinalTime();

			if(force_ext != NULL)
			{
				const SimbodyEngine& engine = getModel().getSimbodyEngine();

				const Body& forceExpressedInBody = getModel().getBodySet().get(force_ext->getForceExpressedInBodyName());
				const Body& pointExpressedInBody = getModel().getBodySet().get(force_ext->getPointExpressedInBodyName());
				const Body& appliedToBody = getModel().getBodySet().get(force_ext->getAppliedToBodyName());
				const Body& groundBody = engine.getGroundBody();

				if(force_ext->appliesForce())
				{
					Vec3 f_val = force_ext->getForceAtTime(time_force);
					engine.transform(s,forceExpressedInBody,f_val,groundBody,f_val);
					Vec3 point(0);
					if(force_ext->specifiesPoint())
					{
						point = force_ext->getPointAtTime(time_force);
						engine.transformPosition(s,pointExpressedInBody,point,appliedToBody,point);
					}

                    smss.addInStationForce(s, appliedToBody.getMobilizedBodyIndex(), point, f_val, totalExternalForces);

				}

				if(force_ext->appliesTorque())
				{
					Vec3 torq_val = force_ext->getTorqueAtTime(time_force);
					engine.transform(s,forceExpressedInBody,torq_val,groundBody,torq_val);
                    smss.addInBodyTorque(s, appliedToBody.getMobilizedBodyIndex(), torq_val, totalExternalForces);
				}
			}

		}

	}

	smss.calcResidualForce(s,Vector(),totalExternalForces,Vector(),Vector(),CminusE);

	Vector bias;
	smss.calcBiasForAccelerationConstraints(s,bias);



	if(bias.size()>0)
	{
		Matrix G_lambda;
		smss.calcG(s,G_lambda);

		Matrix GMInv = G_lambda*MInv;


		Matrix GMInvGT = GMInv*G_lambda.transpose();
		SimTK::FactorQTZ qtz(GMInvGT);
		Matrix invGMInvGT;
		qtz.inverse(invGMInvGT);


		Vector A_lambda = invGMInvGT*(bias-GMInv*CminusE);
		Matrix B_lambda = invGMInvGT*GMInv;

		Matrix I_minus_GTB(_n_u,_n_u);
		I_minus_GTB.setToZero();
		I_minus_GTB.diag().setTo(1.0);

		I_minus_GTB -= G_lambda.transpose()*B_lambda;

	
		//q_ddot = A+B*f;


		//A = MInv*(I_minus_GTB*C_muscles*P_passive-CminusE-G_lambda.transpose()*A_lambda);
		//B.updBlock(0,0,_n_u,_numMuscles) = MInv*I_minus_GTB*C_muscles*N_active_mat;
		//B.updBlock(0,_numMuscles,_n_u,_numNonMuscleActuators) = MInv*I_minus_GTB*_mapping_nonMuscleActuators;

		B.updBlock(0,0,_n_u,_numMuscles) =  MInv*I_minus_GTB*C_muscles;
		A.setToZero();
		A -= MInv*(CminusE+G_lambda.transpose()*A_lambda);
		B.updBlock(0,_numMuscles,_n_u,_numNonMuscleActuators) =  MInv*I_minus_GTB*_mapping_nonMuscleActuators;

	}
	else
	{

		//A = MInv*(C_muscles*P_passive-CminusE);
		//B.updBlock(0,0,_n_u,_numMuscles) = MInv*C_muscles*N_active_mat;
		//B.updBlock(0,_numMuscles,_n_u,_numNonMuscleActuators) = MInv*_mapping_nonMuscleActuators;

		B.updBlock(0,0,_n_u,_numMuscles) = MInv*C_muscles;
		A.setToZero();
		A -= MInv*CminusE;
		B.updBlock(0,_numMuscles,_n_u,_numNonMuscleActuators) = MInv*_mapping_nonMuscleActuators;

	}

}

void ROCINForceController::computeMBSGeneralizedForceDynamics(const SimTK::State& s, Matrix& B, Vector& A)
{
	//q_ddot = A+B*f; //f is the vector of generalized forces

	B.resize(_n_u,_n_u);
	A.resize(_n_u);

	double tiReal = s.getTime();

	const SimTK::SimbodyMatterSubsystem& smss = getModel().getMatterSubsystem();		
	const SimTK::MultibodySystem& mbs = getModel().getMultibodySystem();
	mbs.realize(s,SimTK::Stage::Velocity);


	//Inverse Mass Matrix 
	Matrix MInv(_n_u,_n_u);
	smss.calcMInv(s,MInv);


	Vector CminusE(_n_u);
	const SimTK::Vector_<SimTK::SpatialVec>& appliedGravityOnBody = getModel().getGravityForce().getBodyForces(s);
	SimTK::Vector_<SimTK::SpatialVec> totalExternalForces = appliedGravityOnBody;

	// Figure out number of muscles
	const ForceSet& forceSet = getModel().getForceSet();
	
	
	int col_idx=0;
	for(int i=0;i<forceSet.getSize();i++)
	{
		ExternalForce* force_ext = dynamic_cast<ExternalForce*>(&forceSet.get(i));

		if(force_ext != NULL)
		{
			double time_force = tiReal;
			if(time_force > getFinalTime())
				time_force = getFinalTime();

			if(force_ext != NULL)
			{
				const SimbodyEngine& engine = getModel().getSimbodyEngine();

				const Body& forceExpressedInBody = getModel().getBodySet().get(force_ext->getForceExpressedInBodyName());
				const Body& pointExpressedInBody = getModel().getBodySet().get(force_ext->getPointExpressedInBodyName());
				const Body& appliedToBody = getModel().getBodySet().get(force_ext->getAppliedToBodyName());
				const Body& groundBody = engine.getGroundBody();

				if(force_ext->appliesForce())
				{
					Vec3 f_val = force_ext->getForceAtTime(time_force);
					engine.transform(s,forceExpressedInBody,f_val,groundBody,f_val);
					Vec3 point(0);
					if(force_ext->specifiesPoint())
					{
						point = force_ext->getPointAtTime(time_force);
						engine.transformPosition(s,pointExpressedInBody,point,appliedToBody,point);
					}
                    smss.addInStationForce(s, appliedToBody.getMobilizedBodyIndex(), point, f_val, totalExternalForces);

				}

				if(force_ext->appliesTorque())
				{
					Vec3 torq_val = force_ext->getTorqueAtTime(time_force);
					engine.transform(s,forceExpressedInBody,torq_val,groundBody,torq_val);
                    smss.addInBodyTorque(s, appliedToBody.getMobilizedBodyIndex(), torq_val, totalExternalForces);
				}
			}

		}

	}

	smss.calcResidualForce(s,Vector(),totalExternalForces,Vector(),Vector(),CminusE);

	Vector bias;
	smss.calcBiasForAccelerationConstraints(s,bias);

	if(bias.size()>0)
	{
		Matrix G_lambda;
		smss.calcG(s,G_lambda);

		Matrix GMInv = G_lambda*MInv;


		Matrix GMInvGT = GMInv*G_lambda.transpose();
		SimTK::FactorQTZ qtz(GMInvGT);
		Matrix invGMInvGT;
		qtz.inverse(invGMInvGT);


		Vector A_lambda = invGMInvGT*(bias-GMInv*CminusE);
		Matrix B_lambda = invGMInvGT*GMInv;

		Matrix I_minus_GTB(_n_u,_n_u);
		I_minus_GTB.setToZero();
		I_minus_GTB.diag().setTo(1.0);

		I_minus_GTB -= G_lambda.transpose()*B_lambda;

	
		//q_ddot = A+B*f;


		//A = MInv*(I_minus_GTB*C_muscles*P_passive-CminusE-G_lambda.transpose()*A_lambda);
		//B.updBlock(0,0,_n_u,_numMuscles) = MInv*I_minus_GTB*C_muscles*N_active_mat;
		//B.updBlock(0,_numMuscles,_n_u,_numNonMuscleActuators) = MInv*I_minus_GTB*_mapping_nonMuscleActuators;

		B = MInv*I_minus_GTB;
		A.setToZero();
		A -= MInv*(CminusE+G_lambda.transpose()*A_lambda);


	}
	else
	{

		//A = MInv*(C_muscles*P_passive-CminusE);
		//B.updBlock(0,0,_n_u,_numMuscles) = MInv*C_muscles*N_active_mat;
		//B.updBlock(0,_numMuscles,_n_u,_numNonMuscleActuators) = MInv*_mapping_nonMuscleActuators;

		B.updBlock(0,0,_n_u,_numMuscles) = MInv;
		A.setToZero();
		A -= MInv*CminusE;
	}

}

void ROCINForceController::computeMBSGeneralizedForceDynamicsNumerically(const SimTK::State& s, Matrix& B, Vector& A)
{
	B.resize(_n_u,_n_u);
	A.resize(_n_u);
	B.setToZero();
	A.setToZero();

	SimTK::State s_cpy = s;

	const SimTK::SimbodyMatterSubsystem& smss = getModel().getMatterSubsystem();		
	const SimTK::MultibodySystem& mbs = getModel().getMultibodySystem();

	const Set<Actuator>& fSet = getActuatorSet();
	int nf = fSet.getSize();


	for(int i=0;i<nf;i++)
	{
        ScalarActuator& act = (ScalarActuator&)fSet.get(i);
        act.overrideActuation(s_cpy, true);
        act.setOverrideActuation(s_cpy, 0.0);
	}

	mbs.realize(s_cpy,SimTK::Stage::Dynamics);
	const SimTK::Vector_<SimTK::SpatialVec>& bodyForces = mbs.getRigidBodyForces(s_cpy,SimTK::Stage::Dynamics);




	Vector genForces(_n_u);
	genForces.setToZero();


	SimTK::Vector_<SimTK::SpatialVec> bodyAccs;	

	smss.calcAcceleration(s_cpy,genForces,bodyForces,A,bodyAccs);

	Vector uDot(_n_u);
	uDot.setToZero();

	for(int i=0;i<_n_u;i++)
	{
		genForces[i] = 1.0;
		smss.calcAcceleration(s_cpy,genForces,bodyForces,uDot,bodyAccs);
		B.updCol(i) = uDot-A;
		genForces[i] = 0.0;
	}

}

void ROCINForceController::computeMBSStateDynamics(const SimTK::State& s, const Vector& genForces, const Matrix& muscleJacobian, bool useTaylorExpansion, Matrix& A, Matrix& B, Vector& C)
{
	double tiReal = s.getTime();

	const SimTK::SimbodyMatterSubsystem& smss = getModel().getMatterSubsystem();		
	const SimTK::MultibodySystem& mbs = getModel().getMultibodySystem();
	mbs.realize(s,SimTK::Stage::Velocity);

	int n_x = _n_q+_n_u;
	int n_forces = _numMuscles+_numNonMuscleActuators;

	A.resize(n_x,n_x);
	B.resize(n_x,n_forces);
	C.resize(n_x);

	A.setToZero();
	B.setToZero();
	C.setToZero();

	Vector A_genf;
	Matrix B_genf;
    computeMBSGeneralizedForceDynamicsNumerically(s, B_genf, A_genf);

	Matrix N_qdot;
	computeNqdot(s,N_qdot);
	A.updBlock(0,_n_q,_n_q,_n_u) = N_qdot;

	if(useTaylorExpansion)
	{
		Matrix dUDotdQ, dUDotdU;
		computeGeneralizedForceDynamicsPartialsCentralDiff(s,genForces,dUDotdQ,dUDotdU);

		A.updBlock(_n_q,0,_n_u,_n_q) = dUDotdQ;
		A.updBlock(_n_q,_n_q,_n_u,_n_u) = dUDotdU;

		const Vector& cur_q = s.getQ();
		const Vector& cur_u = s.getU();

		C.updBlock(_n_q,0,_n_u,1) = A_genf-dUDotdQ*cur_q-dUDotdU*cur_u;

	}
	else
		C.updBlock(_n_q,0,_n_u,1) = A_genf;

	B.updBlock(0,0,_n_q,n_forces).setToZero();
	
	B.updBlock(_n_q,0,_n_u,_numMuscles) = B_genf*muscleJacobian;
	B.updBlock(_n_q,_numMuscles,_n_u,_numNonMuscleActuators) = B_genf*_mapping_nonMuscleActuators;


}

void ROCINForceController::computeMBSStateDynamics(const SimTK::State& s, const Vector& actuatorForces, bool useTaylorExpansion, Matrix& A, Matrix& B, Vector& C)
{
	double tiReal = s.getTime();

	const SimTK::SimbodyMatterSubsystem& smss = getModel().getMatterSubsystem();		
	const SimTK::MultibodySystem& mbs = getModel().getMultibodySystem();
	mbs.realize(s,SimTK::Stage::Velocity);

	int n_x = _n_q+_n_u;
	int n_forces = _numMuscles+_numNonMuscleActuators;

	A.resize(n_x,n_x);
	B.resize(n_x,n_forces);
	C.resize(n_x);

	A.setToZero();
	B.setToZero();
	C.setToZero();

	Vector A_f;
	Matrix B_f;
	computeMBSDynamicsNumerically(s,B_f,A_f);

	Matrix N_qdot;
	computeNqdot(s,N_qdot);
	A.updBlock(0,_n_q,_n_q,_n_u) = N_qdot;

	if(useTaylorExpansion)
	{
		Matrix dUDotdQ, dUDotdU;
		computeDynamicsPartialsCentralDiff(s,actuatorForces,dUDotdQ,dUDotdU);
		A.updBlock(_n_q,0,_n_u,_n_q) = dUDotdQ;
		A.updBlock(_n_q,_n_q,_n_u,_n_u) = dUDotdU;

		const Vector& cur_q = s.getQ();
		const Vector& cur_u = s.getU();

		C.updBlock(_n_q,0,_n_u,1) = A_f-dUDotdQ*cur_q-dUDotdU*cur_u;

	}
	else
		C.updBlock(_n_q,0,_n_u,1) = A_f;

	B.updBlock(0,0,_n_q,n_forces).setToZero();
	B.updBlock(_n_q,0,_n_u,n_forces) = B_f;
}

	//x = [q u];
void ROCINForceController::combineSystemDynamics(const Matrix& A_mbs, const Matrix& B_mbs, const Vector& C_mbs, const Vector& A_fm_u, const Vector& B_fm_u, Matrix& A, Matrix& B, Vector& C)
{
	int n_x = A_mbs.nrow();
	int n_controls = B_mbs.ncol();

	A = A_mbs;
	B.resize(n_x,n_controls);
	Matrix N_active_mat(_numMuscles,_numMuscles);
	N_active_mat.setToZero();
	N_active_mat.updDiag() = A_fm_u;
	B.updBlock(0,0,n_x,_numMuscles) = B_mbs.block(0,0,n_x,_numMuscles)*N_active_mat;
	B.updBlock(0,_numMuscles,n_x,_numNonMuscleActuators) = B_mbs.block(0,_numMuscles,n_x,_numNonMuscleActuators);
	C = C_mbs+B_mbs.block(0,0,n_x,_numMuscles)*B_fm_u;

}

void ROCINForceController::computeNqdot(const SimTK::State& s,  Matrix& Nqdot)
{
	const SimTK::SimbodyMatterSubsystem& smss = getModel().getMatterSubsystem();		
	const SimTK::MultibodySystem& mbs = getModel().getMultibodySystem();
	mbs.realize(s,SimTK::Stage::Velocity);

	Nqdot.resize(_n_q,_n_u);
	Vector qdot_tmp(_n_q), u_tmp(_n_u);
	u_tmp.setToZero();
	for(int i=0;i<_n_u;i++)
	{
		u_tmp[i] = 1.0;
		smss.multiplyByN(s,false,u_tmp,qdot_tmp);
		Nqdot.updCol(i) = qdot_tmp;
		u_tmp[i] = 0.0;
	}

}

void ROCINForceController::doPrecomputation_Varying_Dynamics(const SimTK::State& s)
{

//    for (int j = 0; j < _n_u; j++)
//    {
//        Coordinate& coord = getModel().getCoordinateSet().get(j);
//        std::cout << "coord "<<j<<": "<<coord.getName().c_str() << std::endl;
//    }

    //std::cout << "s = " << s << std::endl;
//    PrintVector(s.getQ(), "s.Q", std::cout);
//    PrintVector(s.getU(), "s.U", std::cout);
//    exit(0);
    

	const SimTK::SimbodyMatterSubsystem& smss = getModel().getMatterSubsystem();		
	const SimTK::MultibodySystem& mbs = getModel().getMultibodySystem();
	mbs.realize(s,SimTK::Stage::Velocity);

	Array<Vector> G_array,D_array,E_array;
	Array<Matrix> F_array,H_array;

	int n_dynamics = (_lookahead_number<_mpcSolver->getCurWinSize())? _lookahead_number : _mpcSolver->getCurWinSize();

	if(n_dynamics == 0)
	{		
		_controlSet.setControlValues(s.getTime()+_lookahead_window,&_actuator_controls_new[0]);
		return;
	}


	//int n_dynamics = SimTK::min(_lookahead_number,_mpcSolver->getCurWinSize());


	G_array.setSize(n_dynamics);
	F_array.setSize(n_dynamics);
	H_array.setSize(n_dynamics);

	D_array.setSize(n_dynamics);
	E_array.setSize(n_dynamics);

	Vector y;	

	double dt_ahead = _mpcSolver->getDt();
	int N = _numMuscles+_numNonMuscleActuators;


	setInternalTrajectoryCorrections(s);

	Array<SimTK::State> newState_mbs_array;
	Array<SimTK::State> newState_musc_array;
	newState_mbs_array.setSize(n_dynamics);
	newState_musc_array.setSize(n_dynamics);


	

	Array<double> xmin(0.0,N), xmax(1.0,N);
	for(int i=0;i<_numMuscles;i++)
	{
		xmin[i] = _muscle_min_control_array[i];
		xmax[i] = _muscle_max_control_array[i];
	}

	for(int i=_numMuscles;i<N;i++)
	{
		xmin[i] = _virtual_actuator_min_control_array[i-_numMuscles];
		xmax[i] = _virtual_actuator_max_control_array[i-_numMuscles];
	}

	Array<double> zero(0.0,N);
	Array<double> fmin(0.0,N), fmax(0.0,N);
	Array<double> actsmin(0.0,_numMuscles), actsmax(0.0,_numMuscles);

	Array<Vector> A_fm_u_array;
	Array<Vector> B_fm_u_array;
	Array<Vector> A_act_u_array;
	Array<Vector> B_act_u_array;

	A_fm_u_array.setSize(n_dynamics);
	B_fm_u_array.setSize(n_dynamics);
	A_act_u_array.setSize(n_dynamics);
	B_act_u_array.setSize(n_dynamics);

    // get the array of future states (currently only use joint angles and joint velocities) 
	SimTK::State s_cpy = s;

	for(int i=0;i<n_dynamics;i++)
	{
		_internalSubsys->assignTimeAndQU(s.getTime()+dt_ahead*(i+1),s_cpy);

		newState_mbs_array[i] = s_cpy;

	}

    // compute the muscle relation: fm = A_fm*u+B_fm; a = A_act*u+B_act
	computeMuscleRelation(s,dt_ahead,n_dynamics,A_fm_u_array,B_fm_u_array,A_act_u_array,B_act_u_array);

	double dt_eval = dt_ahead*n_dynamics;
	
	Array<Matrix> muscleJacobian_array;
	muscleJacobian_array.setSize(n_dynamics);

	for(int i=0;i<n_dynamics;i++)
		computeMuscleJacobian(newState_mbs_array[i],muscleJacobian_array[i]);

	for(int i=0;i<n_dynamics;i++)
	{
		Matrix A_mbs, B_mbs;
		Vector C_mbs;

        // compute multibody dynamics
		computeMBSStateDynamics(newState_mbs_array[i],_gen_forces_array[i],muscleJacobian_array[i],_use_taylor_expansion,A_mbs,B_mbs,C_mbs);
        // combine multibody dynamics and muscle relation
		combineSystemDynamics(A_mbs,B_mbs,C_mbs,A_fm_u_array[i],B_fm_u_array[i],F_array[i],H_array[i],G_array[i]);

        //PrintMatrix(A_mbs, "A_mbs", std::cout);
        //PrintMatrix(B_mbs, "B_mbs", std::cout);
        //PrintVector(C_mbs, "C_mbs", std::cout);

        //PrintVector(A_fm_u_array[i], "A_fm_u", std::cout);
        //PrintVector(B_fm_u_array[i], "B_fm_u", std::cout);
        //exit(0);

		
		D_array[i].resize(N);E_array[i].resize(N);
		D_array[i].setTo(1.0);E_array[i].setToZero();

		D_array[i].updBlock(0,0,_numMuscles,1) = A_act_u_array[i];
		E_array[i].updBlock(0,0,_numMuscles,1) = B_act_u_array[i];
	}

	_mpcSolver->setABCArray(F_array,H_array,G_array);
	_mpcSolver->setDiagDandEAarray(D_array,E_array);

	getCurObservation(s,y,false,false,false);

    // use MPC solver to solve control signals
	_mpcSolver->precomputeU(s.getTime(),y);

	Vector cur_u = _mpcSolver->getCurrentU();	

	for(int i=0;i<cur_u.size();i++)
		_actuator_controls_new[i] = cur_u[i];	

	_controlSet.setControlValues(s.getTime()+_lookahead_window,&_actuator_controls_new[0]);

    // store generalized forces into an array
	const Matrix& U_array = _mpcSolver->getUArray();

	for(int i=0;i<n_dynamics;i++)
	{
		_actuator_forces_array[i].updBlock(0,0,_numMuscles,1) = A_fm_u_array[i].elementwiseMultiply(U_array.block(0,i,_numMuscles,1).getAsVector())+B_fm_u_array[i];												
		_actuator_forces_array[i].updBlock(_numMuscles,0,_numNonMuscleActuators,1) = U_array.block(_numMuscles,i,_numNonMuscleActuators,1);
	}

	for(int i=0;i<n_dynamics;i++)
		_gen_forces_array[i] = muscleJacobian_array[i]*_actuator_forces_array[i].block(0,0,_numMuscles,1).getAsVector()+_mapping_nonMuscleActuators*_actuator_forces_array[i].block(_numMuscles,0,_numNonMuscleActuators,1).getAsVector();

}

void ROCINForceController::computeMuscleRelation(const SimTK::State& s, double dt, int n_dynamics, Array<Vector>& A_fm_u_array, Array<Vector>& B_fm_u_array, Array<Vector>& A_act_u_array, Array<Vector>& B_act_u_array)
{
	A_fm_u_array.setSize(n_dynamics);
	B_fm_u_array.setSize(n_dynamics);
	A_act_u_array.setSize(n_dynamics);
	B_act_u_array.setSize(n_dynamics);

	int N = _numMuscles+_numNonMuscleActuators;

	Array<double> zero(0.0,N);
	Array<double> fmin(0.0,N), fmax(0.0,N);
	Array<double> actsmin(0.0,_numMuscles), actsmax(0.0,_numMuscles);

	Array<double> xmin(0.0,N), xmax(1.0,N);
	for(int i=0;i<_numMuscles;i++)
	{
		xmin[i] = _muscle_min_control_array[i];
		xmax[i] = _muscle_max_control_array[i];
	}

	for(int i=_numMuscles;i<N;i++)
	{
		xmin[i] = 0.0;
		xmax[i] = 0.0;
	}

	SimTK::State s_cpy_min = s;
	SimTK::State s_cpy_max = s;

	for(int i=0;i<n_dynamics;i++)
	{
		_internalSubsys->setIntegTargetTime(s.getTime()+dt*(i+1));
		evaluateInternalModelForces(s_cpy_min,&zero[0],&xmin[0],&fmin[0]);
		getInternalModelActivations(_internalSubsys->getCompleteState(),&actsmin[0]);
		s_cpy_min = _internalSubsys->getCompleteState();
		
		evaluateInternalModelForces(s_cpy_max,&zero[0],&xmax[0],&fmax[0]);
		getInternalModelActivations(_internalSubsys->getCompleteState(),&actsmax[0]);
		s_cpy_max = _internalSubsys->getCompleteState();

		solveLinearCoefs(_numMuscles,&xmin[0],&fmin[0],&xmax[0],&fmax[0],A_fm_u_array[i],B_fm_u_array[i]);
		solveLinearCoefs(_numMuscles,&xmin[0],&actsmin[0],&xmax[0],&actsmax[0],A_act_u_array[i],B_act_u_array[i]);

        //PrintArray(fmin, "fmin", std::cout);
        //PrintArray(fmax, "fmax", std::cout);
        //exit(0);

	}

}

void ROCINForceController::getCurObservation(const SimTK::State& s, Vector& y, bool includeMuscleForce, bool includeMuscleFiberLength, bool includeActivation)
{


	int n_y = _n_q+_n_u;
	if(includeMuscleForce)
		n_y += _numMuscles;

	if(includeMuscleFiberLength)
		n_y += _numMuscles;

	if(includeActivation)
		n_y += _numMuscles;

	if(_track_coords)
	{
		y.resize(n_y);
		y.updBlock(0,0,_n_q,1) = s.getQ();
		y.updBlock(_n_q,0,_n_u,1) = s.getU();

		int idx_start = _n_q+_n_u;

		const SimTK::MultibodySystem& mbs = getModel().getMultibodySystem();
		SimTK::State s_cpy = s;

		const ForceSet& forceSet = getModel().getForceSet();	

		if(includeMuscleForce)
		{
			mbs.realize(s_cpy,SimTK::Stage::Dynamics);
			int col_idx=0;
			for(int i=0;i<forceSet.getSize();i++)
			{
				PathActuator* m = dynamic_cast<PathActuator*>(&forceSet.get(i));
				if(m!=NULL)
				{                    
                    y.set(idx_start + col_idx, m->getActuation(s_cpy));                    
					col_idx++;
				}
			}

			idx_start += _numMuscles;
		}

		if(includeMuscleFiberLength)
		{
			mbs.realize(s_cpy,SimTK::Stage::Dynamics);
			int col_idx=0;
			for(int i=0;i<forceSet.getSize();i++)
			{
				Muscle* m = dynamic_cast<Muscle*>(&forceSet.get(i));
				if(m!=NULL)
				{
					y.set(idx_start+col_idx,m->getFiberLength(s_cpy));
					col_idx++;
				}
			}

			idx_start += _numMuscles;
		}

		if(includeActivation)
		{
			mbs.realize(s_cpy,SimTK::Stage::Dynamics);
			int col_idx=0;
			for(int i=0;i<forceSet.getSize();i++)
			{
				PathActuator* m = dynamic_cast<PathActuator*>(&forceSet.get(i));
				if(m!=NULL)
				{
					Muscle* musc = dynamic_cast<Muscle*>(&forceSet.get(i));
					if(musc != NULL)
						y.set(idx_start+col_idx,musc->getActivation(s_cpy));
					else
						y.set(idx_start+col_idx,m->getControl(s_cpy));

					col_idx++;
				}
			}

			idx_start += _numMuscles;

		}


	}

}

// correct the trajectory for the estimated future state, so that the state could transit 
// from the current state to the desired state
void ROCINForceController::setInternalTrajectoryCorrections(const SimTK::State& s)
{
	double dt_ahead = _mpcSolver->getDt();
	
	Array<double> qCorrections(0.0,_n_u), uCorrections(0.0,_n_u);
	Array<double> qSupposed(0.0,_n_u),uSupposed(0.0,_n_u);
	_qSet->evaluate(qSupposed,0,s.getTime());
	_qSet->evaluate(uSupposed,1,s.getTime());

	double maxDeltaq = 0.0;

	for(int i=0;i<_n_u;i++) 
	{
		qCorrections[i] = s.getQ()[i] - qSupposed[i];
		if(fabs(qCorrections[i])>maxDeltaq)
			maxDeltaq = fabs(qCorrections[i]);

		uCorrections[i] = s.getU()[i] - uSupposed[i];
	}
	_internalSubsys->setCoordinateCorrections(&qCorrections[0]);
	_internalSubsys->setSpeedCorrections(&uCorrections[0]);

	double dt_transit = maxDeltaq/5.0;

	if(dt_transit<dt_ahead)
		dt_transit = dt_ahead;

	_internalSubsys->setTransitTargetTime(s.getTime()+dt_transit,dt_transit);

}


// override the computeControls method, this simply assign the precomputed control signals to the actuators
void ROCINForceController::computeControls(const SimTK::State& s, SimTK::Vector& controls) const
{

	for(int iActuator = 0; iActuator<getActuatorSet().getSize(); iActuator++)
	{
			
		Actuator& thisActuator = getActuatorSet().get(iActuator);

		Vector thisActuatorControl(1,_controlSet[iActuator].getControlValue(s.getTime()));

		thisActuator.addInControls(thisActuatorControl,controls);

	}

	return;
}


void ROCINForceController::createExternalLoads(const std::string& aExternalLoadsFileName, Model& aModel)
{

	if(_externalLoads!=NULL)
		delete _externalLoads;
	
	std::string savedCwd = IO::getCwd();
	
	_externalLoads = new ExternalLoads(aModel,aExternalLoadsFileName);	
	_externalLoads->setMemoryOwner(false);
	_externalLoads->invokeConnectToModel(aModel);


	std::string loadKinematicsFileName = _externalLoads->getExternalLoadsModelKinematicsFileName();
	const Storage* loadKinematicsForPointTransformation = NULL;

	IO::TrimLeadingWhitespace(loadKinematicsFileName);
	Storage *temp = NULL;
	// fine if there are no kinematics as long as it was not assigned
	if(!(loadKinematicsFileName == "") && !(loadKinematicsFileName == "Unassigned")){
		temp = new Storage(loadKinematicsFileName);
		if(!temp){
			IO::chDir(savedCwd);
			throw Exception("DynamicsTool: could not find external loads kinematics file '"+loadKinematicsFileName+"'."); 
		}
	}

	// if loading the data, do whatever filtering operations are also specified
	if(temp && _externalLoads->getLowpassCutoffFrequencyForLoadKinematics() >= 0) {
		std::cout<<"\n\nLow-pass filtering coordinates data with a cutoff frequency of "<<_externalLoads->getLowpassCutoffFrequencyForLoadKinematics()<<"."<<std::endl;
		temp->pad(temp->getSize()/2);
		temp->lowpassIIR(_externalLoads->getLowpassCutoffFrequencyForLoadKinematics());
	}
	loadKinematicsForPointTransformation = temp;

	if(loadKinematicsForPointTransformation)
	{
		SimTK::State& s = aModel.initSystem();
		Storage *qStore=NULL;
		Storage *uStore=NULL;

		aModel.getSimbodyEngine().formCompleteStorages(s, *loadKinematicsForPointTransformation,qStore,uStore);
		// qStore should be in radians
		if (qStore->isInDegrees()){
			aModel.getSimbodyEngine().convertDegreesToRadians(*qStore);
		}

		double ti=0.0, tf=0.0;
		if(_track_coords)
		{
			ti = _coordData->getFirstTime();
			tf = _coordData->getLastTime();
		}
		else
		{
			ti = _markerData->getStartFrameTime();
			tf = _markerData->getLastFrameTime();
		}

		_externalLoads->transformPointsExpressedInGroundToAppliedBodies(*qStore, ti, tf);
		delete qStore;
		delete uStore;

		aModel.invalidateSystem();

	}		

	for(int i=0;i<_externalLoads->getSize();i++)
		aModel.updForceSet().adoptAndAppend(&_externalLoads->get(i));

}

void ROCINForceController::addResidualAndReserveActuators(Model& aModel, const std::string& actuatorFiles)
{
	Object::registerType(VirtualActuator());

	VirtualActuatorSet actuatorSet(actuatorFiles);

	actuatorSet.getCoordAndWeightArray(_virtual_actuator_coord_array,_virtual_actuator_weight_array);
	actuatorSet.getMinAndMaxControlArray(_virtual_actuator_min_control_array,_virtual_actuator_max_control_array);

    // add the additional actuators (currently we only support coordinate actuator) to the model
	addVirtualActuators(aModel,_virtual_actuator_coord_array);

}

void ROCINForceController::addVirtualActuators(Model& aModel, const Array<std::string>& coords_virtual)
{
	int n_u = aModel.getNumSpeeds();
	_numNonMuscleActuators = coords_virtual.size();

	_mapping_nonMuscleActuators.resize(n_u,_numNonMuscleActuators);
	_mapping_nonMuscleActuators.setToZero();

	for(int i=0;i<coords_virtual.size();i++)
	{
		std::string coord_name = coords_virtual.get(i);
		
		int idx_coord = aModel.getCoordinateSet().getIndex(coord_name,0);
		Coordinate& coord = aModel.getCoordinateSet().get(idx_coord);

		_mapping_nonMuscleActuators.set(idx_coord,i,1.0);

		CoordinateActuator* newActuator = new CoordinateActuator();
		
		newActuator->setCoordinate(&coord);

		std::string name_actuator = coord_name+"_reserve";

		newActuator->setName(name_actuator);
		newActuator->setOptimalForce(1.0);
		newActuator->setMinControl(-SimTK::Infinity);
		newActuator->setMinControl(SimTK::Infinity);		
		aModel.addForce(newActuator);

	}		
	

}

void ROCINForceController::setUpInternalModel(const Model& aModel)
{
    // make a copy of the model as the internal model
	_internalModel = new Model(aModel);
	
    SimTK::State& si = _internalModel->initSystem();

    // add a ConstantController to the internalModel so that we can easily try different control signals on the internal model
	_internalController = new ConstantController(_internalModel->getActuators().getSize());
	_internalController->setActuators(_internalModel->getActuators());
	_internalModel->addController(_internalController);

    // initialize internal system
	_internalSys = new ROCINActuatorSystem();
	_internalSubsys = new ROCINActuatorSubsystem(*_internalSys,_internalModel);
	_internalSys->realizeTopology();    	

    // initialize internal integrator and internal Manager
	_internalIntegrator = new SimTK::RungeKuttaMersonIntegrator(*_internalSys);
	_internalIntegrator->setAccuracy(1.0e-5);	//1e-5
	_internalIntegrator->setProjectInterpolatedStates(false);

	_internalManager = new Manager(*_internalModel,*_internalIntegrator);
	_internalManager->setSystem(_internalSys);
	_internalManager->setPerformAnalyses(false);
	_internalManager->setWriteToStorage(false);

    _internalModel->initSystem();
    _internalSubsys->setCompleteState(si);

 	return;

}

void ROCINForceController::getInternalModelActivations(const SimTK::State& s, double* acts)
{
	int idx = 0;
	int N = _internalModel->getActuators().getSize();
	for(int iForce = 0; iForce<N; iForce++)
	{
		PathActuator* m = dynamic_cast<PathActuator *>(&_internalModel->getActuators().get(iForce));
		if(m != NULL)
		{
			Muscle* musc = dynamic_cast<Muscle*>(m);
			if(musc!=NULL)
				acts[idx] = musc->getActivation(s);
			else
				acts[idx] = m->getControl(s);
			idx++;
		}
	}

}


void ROCINForceController::evaluateInternalModelForces(const SimTK::State& s, double* targetForces, double* controls, double* forces)
{
	double tiReal = s.getTime();
	int N = _internalModel->getActuators().getSize();

	_internalController->setControls(controls,_internalSubsys->getIntegTargetTime());

	_internalSubsys->setCompleteState(s);

	_internalManager->setInitialTime(tiReal);
	_internalManager->setFinalTime(_internalSubsys->getIntegTargetTime());
	
	SimTK::State actSysState = _internalSys->updDefaultState();
    

	//SimTK::State& actSysState = _internalSys->updDefaultState();
	_internalSubsys->updZ(actSysState) = _internalModel->getMultibodySystem().getDefaultSubsystem().getZ(s);
	actSysState.setTime(tiReal);
	
	
	


	//Vector controlVector(N);
	//for(int i=0;i<N;i++)
	//	controlVector[i] = controls[i];
	//_internalController->setControls(controlVector);




	double t_init = actSysState.getTime();

    SimTK::State s_cpy = s;
	
    _internalModel->getMultibodySystem().realize(s_cpy, SimTK::Stage::Velocity);


    //_internalModel->initSystem();
	_internalManager->integrate(actSysState,1e-6);
    //_internalManager->integrate(s_cpy, 1e-6);



    //_internalManager->integrate(*_internalState, 1e-6);
	//_internalModel->getMultibodySystem().realize(*_internalState,SimTK::Stage::Dynamics);



	for(int iForce = 0; iForce<N; iForce++)
	{

        ScalarActuator& act = (ScalarActuator&)_internalModel->getActuators().get(iForce);
        forces[iForce] = act.getActuation(_internalSubsys->getCompleteState()) - targetForces[iForce];



	}

}


void ROCINForceController::evaluateMuscleForceBasedOnActivation(const Model& aModel, const SimTK::State& s, const double* para, const Vector& activations, Vector& muscleForces)
{
	SimTK::State s_cpy = s;

	const SimTK::SimbodyMatterSubsystem& smss = aModel.getMatterSubsystem();		
	const SimTK::MultibodySystem& mbs = aModel.getMultibodySystem();

	int numMuscles = activations.size();
	Vector lm(numMuscles),lm_dot(numMuscles);
	for(int i=0;i<numMuscles;i++)
	{
		lm[i] = para[i];
		lm_dot[i] = para[numMuscles+i];
	}

	setStateFiberLength(aModel,lm,s_cpy);

	muscleForces.resize(numMuscles);

	mbs.realize(s_cpy,SimTK::Stage::Dynamics);

	int N = aModel.getActuators().getSize();
	int idx_musc = 0;
	for(int j=0;j<N;j++)
	{
		Actuator& act = aModel.getActuators().get(j);
		Thelen2003Muscle* m = dynamic_cast<Thelen2003Muscle*>(&act);

		if(m!=NULL)
		{
			muscleForces[idx_musc] = m->calcActiveFiberForceAlongTendon(activations[idx_musc],lm[idx_musc],lm_dot[idx_musc])+m->getPassiveFiberForceAlongTendon(s_cpy);
			idx_musc++;	
		}

	}

}

void ROCINForceController::rootSolveActivations(const SimTK::State& s, const Vector& muscleForces, const Vector& lm, const Vector& lm_dot, const Vector& acts_min, const Vector& acts_max, const Vector& tol_acts, Vector& activations) const
{
	int numMuscles = lm.size();
	Array<double> para(0.0,numMuscles*2);

	for(int i=0;i<numMuscles;i++)
	{
		para[i] = lm[i];
		para[numMuscles+i]=lm_dot[i];		
	}

	rootSolve(&evaluateMuscleForceBasedOnActivation,*_internalModel,s,muscleForces,acts_min,acts_max,tol_acts,activations,&para[0]);
}

void ROCINForceController::setStateActivation(const Vector& acts,SimTK::State& s) const
{
	setStateActivation(*_internalModel,acts,s);
}

void ROCINForceController::setStateActivation(const Model& aModel, const Vector& acts,SimTK::State& s)
{
	int N = aModel.getActuators().getSize();
	int idx_musc = 0;
	for(int i=0;i<N;i++)
	{
		Actuator& act = aModel.getActuators().get(i);
		Muscle* m = dynamic_cast<Muscle*>(&act);
		if(m!=NULL)
		{
			m->setActivation(s,acts[idx_musc]);
			idx_musc++;	
		}

	}

}
void ROCINForceController::setStateFiberLength(const Vector& fiberLens, SimTK::State& s) const
{
	setStateFiberLength(*_internalModel,fiberLens,s);
}

void ROCINForceController::setStateFiberLength(const Model& aModel, const Vector& fiberLens, SimTK::State& s)
{
	int N = aModel.getActuators().getSize();
	int idx_musc = 0;
	for(int i=0;i<N;i++)
	{
		Actuator& act = aModel.getActuators().get(i);
		ActivationFiberLengthMuscle* m = dynamic_cast<ActivationFiberLengthMuscle*>(&act);
		if(m!=NULL)
		{
			m->setFiberLength(s,fiberLens[idx_musc]);
			idx_musc++;	
		}

	}

}

void ROCINForceController::getMuscleForces(const Model& aModel, const SimTK::State& s, int numMuscles, Vector& muscleForces)
{
	int N = aModel.getActuators().getSize();
	muscleForces.resize(numMuscles);
	int idx_musc = 0;
	for(int i=0;i<N;i++)
	{
		Actuator& act = aModel.getActuators().get(i);
		PathActuator* m = dynamic_cast<PathActuator*>(&act);
		if(m!=NULL)
		{					
			//muscleForces[idx_musc]=m->getForce(s);
            muscleForces[idx_musc] = m->getActuation(s);
			idx_musc++;
		}
	}

}

void ROCINForceController::evaluateMuscleForceBasedOnFiberLength(const Model& aModel,const SimTK::State& s, const double* para, const Vector& fiberLength, Vector& muscleForces)
{
	SimTK::State s_cpy = s;
	const SimTK::SimbodyMatterSubsystem& smss = aModel.getMatterSubsystem();		
	const SimTK::MultibodySystem& mbs = aModel.getMultibodySystem();

	setStateFiberLength(aModel,fiberLength,s_cpy);
	mbs.realize(s_cpy,SimTK::Stage::Dynamics);
	int numMuscles = fiberLength.size();
	getMuscleForces(aModel,s_cpy,numMuscles,muscleForces);

}


void ROCINForceController::rootSolve(void(*evaluate)(const Model&, const SimTK::State&, const double*, const Vector&, Vector&),const Model& aModel,const SimTK::State& s,const Vector& y, const Vector& ax, const Vector& bx, const Vector& tol, Vector& x,const double* para)
{

	int N = ax.size();
	if(y.size() != N)
	{
		std::cout<<"Does not support different size for input and output!"<<std::endl;
		return;
	}

	Vector a = ax;
	Vector b = bx;

	Vector prev_step(N), tol_act(N),p(N),q(N),new_step(N);
	prev_step.setToZero();
	new_step.setToZero();
	tol_act.setToZero();
	p.setToZero();
	q.setToZero();

	Vector f_a_err(N), f_b_err(N);

	evaluate(aModel,s,para,a,f_a_err);
	f_a_err -= y;

	evaluate(aModel,s,para,b,f_b_err);
	f_b_err -= y;

	Vector c = a;
	Vector f_c_err = f_a_err;

	Array<int> converged(0,N);

	bool finished = false;
	int iter;
	for(iter=0;!finished;iter++)
	{
		for(int i=0;i<N;i++)
		{
			if(converged[i]) continue;

			if((f_b_err[i]>0.0 && f_c_err[i]>0.0) || (f_b_err[i]<0.0 && f_c_err[i]<0.0))
			{
				c[i] = a[i];
				f_c_err[i] = f_a_err[i];
			}

			prev_step[i] = b[i]-a[i];

			if(fabs(f_c_err[i])<fabs(f_b_err[i]))
			{
				a[i] = b[i]; b[i] = c[i]; c[i] = a[i];
				f_a_err[i] = f_b_err[i]; f_b_err[i] = f_c_err[i]; f_c_err[i] = f_a_err[i];
			}

			tol_act[i] = 2.0*DBL_EPSILON*fabs(b[i])+0.5*tol[i];
			new_step[i] = 0.5*(c[i]-b[i]);

			if(fabs(new_step[i])<=tol_act[i] || f_b_err[i] == (double)0.0)
			{
				converged[i] = iter;
				continue;
			}

			if(fabs(prev_step[i])>=tol_act[i] && fabs(f_a_err[i])>fabs(f_b_err[i]))
			{
				register double t1,cb,t2;
				cb = c[i]-b[i];

				if(a[i]==c[i])
				{
					t1 = f_b_err[i]/f_a_err[i];
					p[i] = cb*t1;
					q[i] = 1.0-t1;
				}
				else
				{
					q[i] = f_a_err[i]/f_c_err[i]; t1 = f_b_err[i]/f_c_err[i]; t2 = f_b_err[i]/f_a_err[i];
					p[i] = t2*(cb*q[i]*(q[i]-t1)-(b[i]-a[i])*(t1-1.0));
					q[i] = (q[i]-1.0)*(t1-1.0)*(t2-1.0);
				}

				if(p[i]>(double)0.0)
					q[i] = -q[i];
				else
					p[i] = -p[i];

				if(p[i]<(0.75*cb*q[i]-0.5*fabs(tol_act[i]*q[i])) && p[i]<fabs(0.5*prev_step[i]*q[i]))
					new_step[i] = p[i]/q[i];
			}

			if(fabs(new_step[i])<tol_act[i])
			{
				if(new_step[i]>(double)0.0)
					new_step[i] = tol_act[i];
				else
					new_step[i] = -tol_act[i];
			}

			a[i] = b[i]; f_a_err[i] = f_b_err[i];

			b[i] += new_step[i];

		}

		evaluate(aModel,s,para,b,f_b_err);

		f_b_err -= y;

		for(int i=0;i<N;i++)
		{
			finished = true;
			if(!converged[i])
			{
				finished = false;
				break;
			}
		}
	


	}

	x = b;

}

void ROCINForceController::rootSolveMuscleFiberLength(const SimTK::State& s, const Vector& muscleForces, const Vector& lm_min, const Vector& lm_max, const Vector& tol, Vector& muscleFiberLengthVec)
{
	rootSolve(&evaluateMuscleForceBasedOnFiberLength,*_internalModel,s,muscleForces,lm_min,lm_max,tol,muscleFiberLengthVec);
}

void ROCINForceController::getStateActivation(const SimTK::State& s, Vector& acts) const
{
	getStateActivation(*_internalModel,s,_numMuscles,acts);
}

void ROCINForceController::getStateActivation(const Model& aModel, const SimTK::State& s, int numMuscles, Vector& acts)
{
	acts.resize(numMuscles);
	int N = aModel.getActuators().getSize();
	int idx_musc = 0;
	for(int i=0;i<N;i++)
	{
		Actuator& act = aModel.getActuators().get(i);
		PathActuator* m = dynamic_cast<PathActuator*>(&act);
		if(m!=NULL)
		{
			Muscle* musc = dynamic_cast<Muscle*>(m);
			if(musc!=NULL)
				acts[idx_musc] = musc->getActivation(s);
			else
				acts[idx_musc] = m->getControl(s);

			idx_musc++;	
		}

	}

}
void ROCINForceController::getStateFiberLength(const SimTK::State& s, Vector& fiberLens) const
{
	getStateFiberLength(*_internalModel,s,_numMuscles,fiberLens);
}

void  ROCINForceController::getStateFiberLength(const Model& aModel, const SimTK::State& s, int numMuscles,Vector& fiberLens)
{
	fiberLens.resize(numMuscles);
	int N = aModel.getActuators().getSize();
	int idx_musc = 0;
	for(int i=0;i<N;i++)
	{
		Actuator& act = aModel.getActuators().get(i);
		Muscle* m = dynamic_cast<Muscle*>(&act);
		if(m!=NULL)
		{
			fiberLens[idx_musc]=m->getFiberLength(s);
			idx_musc++;	
		}

	}

}
