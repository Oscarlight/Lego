/*
 * InOut.h
 *
 *  Created on: Oct 15, 2015
 *      Author: oscar
 */

#ifndef INOUT_H_
#define INOUT_H_

#include "lib/getpot-c++/GetPot"

class InOut {
private:

public:
	InOut();
	virtual ~InOut();

	int readInput(int argc, char* argv[]);
};

#endif /* INOUT_H_ */
