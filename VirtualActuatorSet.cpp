#include "VirtualActuatorSet.h"

using namespace OpenSim;

VirtualActuator::VirtualActuator()
{
	constructProperties();
}

VirtualActuator::VirtualActuator(const std::string& coord, double weight, double minVal, double maxVal)
{
	constructProperties();

	updProperty_coordinate() = coord;
	updProperty_weight() = weight;
	updProperty_min_control() = minVal;
	updProperty_max_control() = maxVal;
}

void VirtualActuator::constructProperties()
{
	constructProperty_coordinate("");
	constructProperty_weight(1.0);
	constructProperty_min_control(-200.0);
	constructProperty_max_control(200.0);
}


VirtualActuatorSet::~VirtualActuatorSet()
{
}

VirtualActuatorSet::VirtualActuatorSet() :
Set<VirtualActuator>()
{
}

VirtualActuatorSet::VirtualActuatorSet(const std::string& aFileName) :
	Set<VirtualActuator>(aFileName,false)
{
	updateFromXMLDocument();
}

void VirtualActuatorSet::getCoordinateArray(Array<std::string>& coord_array) const
{
	int size = getSize();
	coord_array.setSize(size);
	for(int i=0;i<size;i++)
	{
		VirtualActuator& act = get(i);
		coord_array[i] = act.getProperty_coordinate().getValue();
	}
}

void VirtualActuatorSet::getWeightArray(Array<double>& weight_array) const
{
	int size = getSize();
	weight_array.setSize(size);
	for(int i=0;i<size;i++)
	{
		VirtualActuator& act = get(i);
		weight_array[i] = act.getProperty_weight().getValue();
	}
}

void VirtualActuatorSet::getCoordAndWeightArray(Array<std::string>& coord_array, Array<double>& weight_array) const
{
	int size = getSize();
	coord_array.setSize(size);
	weight_array.setSize(size);

	for(int i=0;i<size;i++)
	{
		VirtualActuator& act = get(i);
		coord_array[i] = act.getProperty_coordinate().getValue();
		weight_array[i] = act.getProperty_weight().getValue();
	}
}

void VirtualActuatorSet::getMinAndMaxControlArray(SimTK::Vector& min_control_array, SimTK::Vector& max_control_array) const
{
	int size = getSize();
	min_control_array.resize(size);
	max_control_array.resize(size);

	for(int i=0;i<size;i++)
	{
		VirtualActuator& act = get(i);
		min_control_array[i] = act.getProperty_min_control().getValue();
		max_control_array[i] = act.getProperty_max_control().getValue();
	}
}