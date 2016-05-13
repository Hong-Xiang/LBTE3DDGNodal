#pragma once
#include <exception>

namespace IMATLABC {
	void iMC_readMATFile(double* data, const char* file_name, const char* variable_name, const int possible_index = 0, const bool using_index = false);

//Exceptions:
	class MATFileOpenFail : std::exception {};
}