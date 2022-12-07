#include "MatrixArithmetics.h"
using namespace arrmath;

int main(){
	
	vector<double> test1 = vectorOfZeros(10);
	if (test1 == {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }){
		std::cout << "vectorOfZeros OK." << std::endl;
	} else {
		std::cout << "vectorOfZeros FAILED." << std::endl;
	}








	return 0;
}