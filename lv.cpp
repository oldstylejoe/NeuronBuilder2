//Joe Snider
//1/08
//
//Calcualte the K4/K2 estimate for slopes from product_moment files.

#include <iostream>
#include <fstream>
#include <vector>

#include "../utilities/linear_regression.h"

using namespace std;

int main() {
	float labelValue;
	char junkString[1000];
	double x, y;
	char fileName[1000];

	while(cin >> fileName) {
		cerr << "Reading from " << fileName << "..." << flush;
		vector<double> dataX, dataY, errorY;
		sscanf(fileName, "length_volume_vary%f.txt", &labelValue);
		ifstream inputFile(fileName);
		while(inputFile >> y >> x) {
			dataX.push_back(log(x));
			dataY.push_back(log(y));
			errorY.push_back(0.);
		}
		inputFile.close();

		double slopeXY, slopeErrorXY, slopeYX, slopeErrorYX, offset, offsetError;
		LinearRegression(dataX, dataY, errorY, slopeXY, slopeErrorXY, offset, offsetError);
		LinearRegression(dataY, dataX, errorY, slopeYX, slopeErrorYX, offset, offsetError);
slopeYX = 1.;

		cout << labelValue << " "
			<< sqrt(slopeXY/slopeYX) << " "
			<< slopeErrorXY << "\n" << flush;

		cerr << "done\n" << flush;
	}
}
