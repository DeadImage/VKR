#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "wavelets.h"

using namespace std;

vector<double> dwt_lift_1d(vector<double> &initialValues, unsigned int numOfSteps) {
	
	vector<double> copyOfInitialValues = initialValues;
	if (numOfSteps == 0) {
		return copyOfInitialValues;
	}

	unsigned int powerOfTwo = 0;
	unsigned int steps = 0;

	size_t numOfValues = initialValues.size();
	vector<double> transformedValues(numOfValues);

	if (!isPowerOfTwo(numOfValues)) {
		fprintf(stderr, "Error - the size of array must be a power of 2.\n");
		exit(1);
	}

	steps = maxNumberOfSteps_1D(numOfValues);
	if (numOfSteps > steps) {
		fprintf(stderr, "Error - the number of steps is higher than possible.\n");
		exit(1);
	}
	
	powerOfTwo = pow(2, steps - 1);
	steps = numOfSteps;
	while (steps > 0) {
		for (unsigned int i = 0; i < powerOfTwo; i++) {
			transformedValues.at(i) = (copyOfInitialValues.at(2*i) + copyOfInitialValues.at(2*i+1)) / 2;
			transformedValues.at(i+powerOfTwo) = (copyOfInitialValues.at(2 * i) - copyOfInitialValues.at(2*i+1)) / 2;
		}

		for (unsigned int i = 0; i < powerOfTwo; i++) {
			copyOfInitialValues.at(i) = transformedValues.at(i);
		}
		steps--;
		powerOfTwo /= 2;
	}

	return transformedValues;
}

vector<vector<double>> dwt_lift_2d (vector<vector<double>>& initialValues, unsigned int numOfSteps) {

	vector<vector<double>> copyOfInitialValues = initialValues;
	unsigned int steps;
	unsigned int powerOfTwoRows;
	unsigned int powerOfTwoCols;

	if (numOfSteps == 0) {
		return copyOfInitialValues;
	}

	size_t numOfRows = initialValues.size();
	size_t numOfCols = initialValues[0].size();

	vector<vector<double>> transformedValues(numOfRows, vector<double>(numOfCols));

	if (!isPowerOfTwo(numOfRows) || !isPowerOfTwo(numOfCols)) {
		fprintf(stderr, "Error - the number of rows and columns must be powers of 2.\n");
		exit(1);
	}

	steps = maxNumberOfSteps_2D(numOfRows, numOfCols);
	if (numOfSteps > steps) {
		fprintf(stderr, "Error - the number of steps is higher than possible.\n");
		exit(1);
	}

	powerOfTwoRows = pow(2, maxNumberOfSteps_1D(numOfCols) - 1);
	powerOfTwoCols = pow(2, maxNumberOfSteps_1D(numOfRows) - 1);
	steps = numOfSteps;

	while (steps > 0) {

		for (unsigned int i = 0; i < numOfRows; i++) {
			for (unsigned int j = 0; j < powerOfTwoRows; j++) {
				transformedValues[i][j] = (copyOfInitialValues[i][j * 2] + copyOfInitialValues[i][2 * j + 1]) / 2;
				transformedValues[i][j + powerOfTwoRows] = (copyOfInitialValues[i][j * 2] - copyOfInitialValues[i][2 * j + 1]) / 2;
			}
		}

		for (unsigned int i = 0; i < numOfRows; i++) {
			for (unsigned int j = 0; j < numOfCols; j++)
				copyOfInitialValues[i][j] = transformedValues[i][j];
		}

		for (unsigned int i = 0; i < numOfCols; i++) {
			for (unsigned int j = 0; j < powerOfTwoCols; j++) {
				double temp1 = copyOfInitialValues[j * 2][i];
				double temp2 = copyOfInitialValues[2 * j + 1][i];

				transformedValues[j][i] = (temp1 + temp2) / 2;
				transformedValues[j + powerOfTwoCols][i] = (temp1 - temp2) / 2;
			}
		}

		for (unsigned int i = 0; i < powerOfTwoRows; i++) {
			for (unsigned int j = 0; j < powerOfTwoCols; j++)
				copyOfInitialValues[i][j] = transformedValues[i][j];
		}

		steps--;
		powerOfTwoRows /= 2;
		powerOfTwoCols /= 2;
		numOfRows /= 2;
		numOfCols /= 2;
	}

	return transformedValues;
}

vector<double> dwt_1d_filter(vector<double>& initialValues, unsigned int numOfSteps, double h[], double g[], short order) {
	unsigned int steps = 0;

	size_t numOfValues = initialValues.size();
	vector<double> transformedValues = initialValues;

	steps = maxNumberOfSteps_1D(numOfValues);
	if (numOfSteps > steps) {
		fprintf(stderr, "Error - the number of steps is higher than possible.\n");
		exit(1);
	}

	if (!isPowerOfTwo(numOfValues)) {
		fprintf(stderr, "Error - the size of array must be a power of 2.\n");
		exit(1);
	}

	int powerOfTwo = pow(2, steps);
	steps = numOfSteps;
	while (steps > 0) {
		vector<double> copyOfInitialValues(powerOfTwo);
		for (int i = 0; i < powerOfTwo; i++) {
			copyOfInitialValues[i] = transformedValues[i];
		}

		vector<double> tempTransformedValues = convolution(copyOfInitialValues, h, g, order);
		for (int i = 0; i < powerOfTwo; i++) {
			transformedValues[i] = tempTransformedValues[i];
		}
		powerOfTwo /= 2;
		steps--;
	}

	return transformedValues;
}

vector<vector<double>> dwt_2d_filter(vector<vector<double>>& initialValues, 
	unsigned int numOfSteps, double h[], double g[], short order) {
	
	unsigned int steps = 0;

	size_t numOfRows = initialValues.size();
	size_t numOfCols = initialValues[0].size();

	vector<vector<double>> transformedValues = initialValues;

	if (!isPowerOfTwo(numOfRows) || !isPowerOfTwo(numOfCols)) {
		fprintf(stderr, "Error - the number of rows and columns must be powers of 2.\n");
		exit(1);
	}

	steps = maxNumberOfSteps_2D(numOfRows, numOfCols);
	if (numOfSteps > steps) {
		fprintf(stderr, "Error - the number of steps is higher than possible.\n");
		exit(1);
	}

	int powerOfTwoRows = pow(2, maxNumberOfSteps_1D(numOfCols));
	int powerOfTwoCols = pow(2, maxNumberOfSteps_1D(numOfRows));
	steps = numOfSteps;

	while (steps > 0) {
		for (unsigned int i = 0; i < numOfRows; i++) {
			vector<double> tempTransformedValues(numOfCols);
			for (unsigned int j = 0; j < numOfCols; j++) {
				tempTransformedValues[j] = transformedValues[i][j];
			}
			tempTransformedValues = convolution(tempTransformedValues, h, g, order);
			for (int j = 0; j < powerOfTwoRows; j++) {
				transformedValues[i][j] = tempTransformedValues[j];
			}
		}

		for (unsigned int i = 0; i < numOfCols; i++) {
			vector<double> tempTransformedValues(numOfRows);
			for (unsigned int j = 0; j < numOfRows; j++) {
				tempTransformedValues[j] = transformedValues[j][i];
			}
			tempTransformedValues = convolution(tempTransformedValues, h, g, order);
			for (int j = 0; j < powerOfTwoCols; j++) {
				transformedValues[j][i] = tempTransformedValues[j];
			}
		}

		steps--;
		powerOfTwoRows /= 2;
		powerOfTwoCols /= 2;
		numOfRows /= 2;
		numOfCols /= 2;
	}

	return transformedValues;
}

vector<double> dwt_haar_1d(vector<double>& initialValues, unsigned int numOfSteps) {
	if (numOfSteps == 0) {
		return initialValues;
	}

	double h[] = { 0.7071, 0.7071 };
	double g[] = { 0.7071, -0.7071 };
	vector<double>transformedValues = dwt_1d_filter(initialValues, numOfSteps, h, g, 2);

	return transformedValues;
}

vector<vector<double>> dwt_haar_2d(vector<vector<double>>& initialValues, unsigned int numOfSteps) {
	if (numOfSteps == 0) {
		return initialValues;
	}

	double h[] = { 0.5, 0.5 };
	double g[] = { 0.5, -0.5 };
	vector<vector<double>>transformedValues = dwt_2d_filter(initialValues, numOfSteps, h, g, 2);

	return transformedValues;
}

int isPowerOfTwo(unsigned int x) {
	return ((x != 0) && ((x & (~x + 1)) == x));
}

unsigned int maxNumberOfSteps_1D(size_t numOfValues) {
	unsigned int powerOfTwo = 1;
	unsigned int steps = 0;
	while (powerOfTwo < numOfValues) {
		powerOfTwo = powerOfTwo * 2;
		steps++;
	}

	return steps;
}

unsigned int maxNumberOfSteps_2D(size_t numOfRows, size_t numOfCols) {
	unsigned int maxRows = maxNumberOfSteps_1D(numOfRows);
	unsigned int maxCols = maxNumberOfSteps_1D(numOfCols);

	return min(maxRows, maxCols);
}

vector<double> convolution(vector<double>& input, double h[], double g[], short order) {
	vector<double> output(input.size());
	short steps = input.size() / order;
	double* buffer = new double[order];

	for (unsigned int step = 0; step < steps; step++) {
		for (unsigned int i = 0; i < order; i++) {
			buffer[i] = input[step*order + i];
		}

		output[step] = h[0] * buffer[0];
		output[step + (steps * order / 2)] = g[0] * buffer[0];

		for (unsigned int i = order - 1; i > 0; i--) {
			output[step] += h[i] * buffer[i];
			output[step + (steps * order / 2)] += g[i] * buffer[i];
		}
	}

	return output;
}