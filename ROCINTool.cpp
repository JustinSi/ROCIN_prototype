#include "ROCINTool.h"
#include <OpenSim/OpenSim.h>
#include "TestHelper.h"
#include "ControlEvaluationAnalysis.h"

using namespace OpenSim;


bool isFileExist(const std::string& name)
{
	std::ifstream f(name.c_str());
	if(f.good())
	{
		f.close();
		return true;
	}
	else
	{
		f.close();
		return false;
	}
}

ROCINTool::ROCINTool(const std::string& setUpFile)
{
	_controller = NULL;
	constructProperties();
	XMLDocument document(setUpFile);
	SimTK::Xml::Element myNode = document.getRootDataElement();
    //read in the control setup XML file
	readObjectFromXMLNodeOrFile(myNode,document.getDocumentVersion());

}


ROCINTool::~ROCINTool()
{
	if(_controller != NULL)
		delete _controller;
}

void ROCINTool::evalCoordTrackingErrs(const std::string& modelFile, const std::string& desired_kinematics_file, const std::string& state_file, const std::string& err_file)
{
	Model osimModel(modelFile);
	osimModel.initSystem();
	ROCINForceController::evalModelCoordTrackingErrs(osimModel,desired_kinematics_file,state_file,err_file);
}

bool ROCINTool::run()
{
    //read the model file
	std::string modelFileName = getProperty_model_file().getValue();
	if(!isFileExist(modelFileName))
	{
		std::cout<<"Model file "<<modelFileName.c_str()<<" does not exist!"<<std::endl;
		exit(1);
	}

	Model osimModel(getProperty_model_file().getValue());
	clock_t time_start = clock(), diff;


	std::string desired_kinematics_file = getProperty_desired_kinematics_file().getValue();

	if(!isFileExist(desired_kinematics_file))
	{
		std::cout<<"Desired kinematics file "<<desired_kinematics_file.c_str()<<" does not exist!"<<std::endl;
		exit(1);
	}

    //construct the controller
	_controller = new ROCINForceController(osimModel,desired_kinematics_file,true,getProperty_update_horizon().getValue(),getProperty_prediction_step().getValue());

	if(getProperty_use_taylor_expansion().getValue())
	{
		if(getProperty_use_implicit_mpc_integ().getValue() == false)
		{
			std::cout<<"You must use implicit mpc integration if you use taylor expansion!!!"<<std::endl;
			exit(1);
		}
	}

	if(getProperty_prediction_step().getValue()>1)
	{
		std::cout<<"We suggest that you use implicit mpc integration and taylor expansion if you have prediction step larger than 1!!!"<<std::endl;		
	}

    //set the penalty weight
	_controller->setCoordTrackingPenalty(getProperty_coord_tracking_penalty().getValue());
	_controller->setSpeedTrackingPenalty(getProperty_speed_tracking_penalty().getValue());
	_controller->setActivationPenalty(getProperty_activation_penalty().getValue());
	_controller->setPDPenalty(getProperty_PD_penalty().getValue());

	_controller->setUseTaylorExpansion(getProperty_use_taylor_expansion().getValue());
	_controller->setLowPassCutOffFrequency(getProperty_lowpass_cutoff_frequency().getValue());

	if(!getProperty_external_loads_file().empty())
	{
		std::string externalLoadsFileName = getProperty_external_loads_file().getValue();
		if(!isFileExist(externalLoadsFileName))
		{
			if(externalLoadsFileName != "")
			{
				std::cout<<"External loads file "<<externalLoadsFileName.c_str()<<" does not exist!"<<std::endl;
				exit(1);
			}
			std::cout<<"No external loads are added to the model."<<std::endl;
		}
        // read in the external loads file
		else
			_controller->createExternalLoads(externalLoadsFileName,osimModel);
	}

	
	if(!getProperty_additional_actuator_file().empty())
	{

		std::string additionalActuatorFileName = getProperty_additional_actuator_file().getValue();
		if(!isFileExist(additionalActuatorFileName))
		{
			if(additionalActuatorFileName != "")
			{
				std::cout<<"Additional actuator file "<<additionalActuatorFileName.c_str()<<" does not exist!"<<std::endl;
				exit(1);
			}
			std::cout<<"No additional actuaors are added to the model."<<std::endl;
		}
        //read in the additional actuators (residual and reserve actuators)
		else
			_controller->addResidualAndReserveActuators(osimModel,additionalActuatorFileName);
	}

    //set up the internal model
	_controller->setUpInternalModel(osimModel);
    //add all the actuators of the model to the controller
	_controller->setActuators(osimModel.updActuators());
    //add the controller to the model
	osimModel.addController(_controller);
	
    //add an Analysis (for control evaluation) to the model
	ControlEvaluationAnalysis* controlAnalysis = new ControlEvaluationAnalysis(&osimModel);	
	osimModel.addAnalysis(controlAnalysis);

	double initTime = getProperty_initial_time().getValue();

	SimTK::State * si_p = NULL;
	
	if(getProperty_use_init_state_file().getValue())
	{
		std::string initStateFileName = getProperty_init_state_file().getValue();
		if(!isFileExist(initStateFileName))
		{
			std::cout<<"Initial state file "<<initStateFileName.c_str()<<" does not exist!"<<std::endl;
			exit(1);
		}

        //initialize the state using a state file
		si_p = &(_controller->setInitStateFromFile(osimModel,initStateFileName));
		initTime = si_p->getTime();
	}
	else
	{
        //initialize without a state file
		si_p = &(_controller->setInitState(osimModel,initTime));
	}

    //initialize the moment arm solver, this is used when computing muscle Jacobian
	controlAnalysis->initMomentArmSolver();



	SimTK::State& si = *si_p;
	
    //set the integrator
	SimTK::RungeKuttaMersonIntegrator integrator(osimModel.getMultibodySystem());

	integrator.setAccuracy(getProperty_integrator_error_tolerance().getValue());
	integrator.setMaximumStepSize(getProperty_maximum_integrator_step_size().getValue());
	Manager manager(osimModel,integrator);

	double finalTime = getProperty_final_time().getValue();
	manager.setInitialTime(initTime);
	manager.setFinalTime(finalTime);

	_controller->setUseImplicitIntegForMPC(getProperty_use_implicit_mpc_integ().getValue());
	_controller->setPenalizeOutputDerivate(getProperty_add_PD_penalty().getValue());
	_controller->setPenalizeControlDerivative(false);
	_controller->setModelNaturalFrequency(getProperty_natural_frequency().getValue());
	_controller->setUsingVaryingDynamics(true);

    //call doControlPrecomputation (this is the core function that compute the control signals), calling this explicitly so that the control values are 
    //also specified at the starting point
	_controller->doControlPrecomputation(si);

	_controller->setTargetTime(si.getTime()+_controller->getTargetDT());

    // do simulation
	manager.integrate(si);

	diff = clock()-time_start;
	int msec = diff*1000/CLOCKS_PER_SEC;
	std::cout<<"ROCIN Running time: "<<double(msec)/1000.0<<" sec"<<std::endl;

    //write output files
	IO::makeDir(getProperty_results_directory().getValue());

	std::string control_file = getProperty_results_directory().getValue()+"/ROCIN_controls.sto";
	std::string state_file = getProperty_results_directory().getValue()+"/ROCIN_states.sto";
    //write the control signals
	osimModel.printControlStorage(control_file);
    //write the state 
	manager.getStateStorage().print(state_file);
    //write the coord tracking error
	std::string err_file = getProperty_results_directory().getValue()+"/err_coord_tracking.sto";
	_controller->evalCoordTrackingErrs(desired_kinematics_file,state_file,err_file);
	//write the control analysis result
	std::string report_file = getProperty_results_directory().getValue()+"/control_analysis_report.sto";
	controlAnalysis->getAnalysisStorage().print(report_file);



	return true;
}




void ROCINTool::constructProperties()
{
	constructProperty_model_file("");
	constructProperty_additional_actuator_file("");
	constructProperty_results_directory("");
	constructProperty_use_init_state_file(false);
	constructProperty_init_state_file("");
	constructProperty_initial_time(0.0);
	constructProperty_final_time(0.0);
	constructProperty_maximum_integrator_step_size(1.0);
	constructProperty_integrator_error_tolerance(1.0e-5);
	constructProperty_external_loads_file("");
	constructProperty_desired_kinematics_file("");
	constructProperty_lowpass_cutoff_frequency(-1.0);
	constructProperty_use_implicit_mpc_integ(true);
	constructProperty_add_PD_penalty(false);
	constructProperty_use_taylor_expansion(false);
	constructProperty_update_horizon(0.01);
	constructProperty_prediction_step(1);
	constructProperty_coord_tracking_penalty(50.0);
	constructProperty_speed_tracking_penalty(50.0);
	constructProperty_activation_penalty(1.0);
	constructProperty_PD_penalty(0.01);
	constructProperty_natural_frequency(10.0);
}