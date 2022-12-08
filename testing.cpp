#include "MatrixArithmetics.h"
using namespace arrmath;

void printResult(const char name[], bool test) {
	if (test) {
		std::cout << name << " OK." << std::endl;
	}
	else {
		std::cout << name << " FAILED." << std::endl;
	}
}

int main(){
	const double PI = 3.14159265358979323846;

	vector<double> test1 = vectorOfZeros<double>(10);
	printResult("vectorOfZeros", test1 == vector<double>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 });

	test1[3] = 10.0;
	double test2 = vectorMax(test1);
	printResult("vectorMax", test2 == 10.0);

	test1[9] = -10.0;
	double test3 = vectorMin(test1);
	printResult("vectorMin", test3 == -10.0);

	double test4 = vectorSpan(test1);
	printResult("vectorSpan", test4 == 20.0);

	double test5 = vectorMaxIndex(test1);
	printResult("vectorMaxIndex", test5 == 3);

	double test6 = vectorMinIndex(test1);
	printResult("vectorMinIndex", test6 == 9);

	test1[6] = 10.0;
	vector<double> test7 = vectorZeroAvg(test1);
	printResult("vectorZeroAvg", test7 == vector<double>{0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 5.0, 0.0, 0.0, -10.0 });

	vector<double> test8 = vectorElemMult(test1, vector<double>{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0});
	printResult("vectorElemMult", test8 == vector<double>{0.0, 0.0, 0.0, 30.0, 0.0, 0.0, 60.0, 0.0, 0.0, -90.0 });

	vector<double> test9 = vectorElemAdd(test1, vector<double>{0.0, 1.0, 2.0, -7.0, 4.0, 5.0, -4.0, 7.0, 8.0, 19.0});
	printResult("vectorElemAdd", test9 == vector<double>{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 });
	
	vector<double> test10 = vectorScale(test1, 0.5); 
	printResult("vectorScale", test10 == vector<double>{0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 5.0, 0.0, 0.0, -5.0 });

	matrix<double> test11 = vectorTranspose(test1);
	printResult("vectorTranspose", test11 == matrix<double>{ { 0.0 }, { 0.0 }, { 0.0 }, { 10.0 }, { 0.0 }, { 0.0 }, { 10.0 }, { 0.0 }, { 0.0 }, { -10.0 } });

	double test12 = vectorDotProduct(test1, vector<double>{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0});
	printResult("vectorDotProduct", test12 == 0.0);

	double test13 = vectorMagnitude(vector<double>{2.0, 3.0, 4.0, 5.0, 6.0, 3.0, 1.0, 0.0, 0.0, 0.0});
	printResult("vectorMagnitude", test13 == 10.0);

	vector<double> test14 = vectorNormalise(vector<double>{2.0, 3.0, 4.0, 5.0, 6.0, 3.0, 1.0, 0.0, 0.0, 0.0});
	double testSum = 0.0;
	for (size_t i{}; i < test14.size(); i++) { testSum += test14[i]; }
	printResult("vectorNormalise", (float)testSum == 2.4f);

	vector<double> test15 = vectorAverage(vector<double>{1.0, 2.0, 3.0}, vector<double>{1.0, -2.0, 1.0});
	vector<double> test16 = vectorAverage(matrix<double>{ {1.0, 2.0, 3.0}, {1.0, 1.0, 1.0}, {4.0, 3.0, 2.0} });
	printResult("vectorAverage", test15 == vector<double>{1.0, 0.0, 2.0} && test16 == vector<double>{2.0, 2.0, 2.0});

	double test17 = vectorGetDistance(vector<double>{1.0, 0.0, 0.0}, vector<double>{1.0, 0.0, 1.0});
	double test18 = vectorGetDistance(vector<double>{0.0, 0.5}, vector<double>{PI, 0.5}, CoordSystem::pol);
	double test19 = vectorGetDistance(vector<double>{0.0, 0.0, 0.5}, vector<double>{0.0, PI, 0.5}, CoordSystem::sph);
	printResult("vectorGetDistance", test17 == 1.0 && test18 == 1 && test19 == 1);
	
	vector<double> test20 = vectorFunc<double>(test1, [](double x){return 2 * x - 1; });
	printResult("vectorFunc", test20 == vector<double>{-1.0, -1.0, -1.0, 19.0, -1.0, -1.0, 19.0, -1.0, -1.0, -21.0 });


	scrt::matrixFun();

	return 0;
}