// Petsc Libs
#include "petscdmda.h"
#include "petscksp.h"
#include "petscpc.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>

// My Libraries
#include <libNavierStokes/Message.hpp>

using namespace std;
using namespace bitpit;

int main(int argc, char** argv)
{
	Message test;
	test.display();

	return 0;
}