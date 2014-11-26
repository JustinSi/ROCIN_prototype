#include "ROCINTool.h"
#include <OpenSim/OpenSim.h>
#include <OpenSim/Simulation/Manager/Manager.h>
#include "TestHelper.h"
#include <time.h>

using namespace OpenSim;



int main(int argc, const char* argv[])
{
	if(argc!=2 && argc!=6)
	{
		std::cout<<"Bad arguments !"<<std::endl;
		std::cout<<"Usage: ROCIN [setUpFile] "<<std::endl;
		std::cout<<"   Or: ROCIN eval [modelFile] [file1] [file2] [errFile] "<<std::endl;
		exit(1);
	}

    // This is generally how we call ROCIN controller
	if(argc == 2)
	{
		std::string filename = argv[1];


		ROCINTool* rocinTool = new ROCINTool(filename);

		rocinTool->run();
	}
    //This is to print out coord tracking error given two files (e.g. desired kinematics file and a state file)
    else
	{
		if(strcmp(argv[1],"eval")!=0)
		{
			std::cout<<"Unrecognized option "<<argv[1]<<"!!!"<<std::endl;
			exit(1);
		}

		std::string modelFile = argv[2];
		std::string file1 = argv[3];
		std::string file2 = argv[4];
		std::string errFile = argv[5];

		ROCINTool::evalCoordTrackingErrs(modelFile,file1,file2,errFile);
	}
		
	return 0;
}