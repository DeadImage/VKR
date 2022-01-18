#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "wavelets.h"

using namespace std;

void dwt_1d_filter(uint8_t initialValues[], int numOfValues, double h[], double g[], short order) {
	double* transformedValues = new double [numOfValues];
	for (int i = 0; i < numOfValues; i++)
		transformedValues[i] = (double)initialValues[i];

	if (!isPowerOfTwo(numOfValues)) {
		fprintf(stderr, "Error - the size of array must be a power of 2.\n");
		exit(1);
	}

	convolution(transformedValues, numOfValues, h, g, order);

	for (int i = 0; i < numOfValues; i++)
		initialValues[i] = (uint8_t)transformedValues[i];
}

void dwt_2d_filter(uint8_t initialValues[], int rows, int cols, double h[], double g[], short order) {

	double* transformedValues = new double[rows * cols];
	for (int i = 0; i < rows * cols; i++)
		transformedValues[i] = (double)initialValues[i];

	for (int i = 0; i < rows; i++) {
		convolution(&transformedValues[i * cols], cols, h, g, order);
	}

	double* temp = new double[rows];

	for (int i = 0; i < cols; i++) {
		for (int j = 0; j < rows; j++)
			temp[j] = transformedValues[i + j * cols];

		convolution(temp, rows, h, g, order);
		for (int j = 0; j < rows; j++)
			transformedValues[i + j * cols] = temp[j];
	}

	for (int i = 0; i < rows * cols; i++)
		initialValues[i] = (uint8_t)transformedValues[i];
}

void dwt_haar_1d(uint8_t initialValues[], int n, unsigned int numOfSteps) {
	if (numOfSteps > maxNumberOfSteps_1D(n)) {
		fprintf(stderr, "Error - the number of steps is higher than possible.\n");
		exit(1);
	}

	int step_length = n;

	double h[] = { 0.7071, 0.7071 };
	double g[] = { 0.7071, -0.7071 };

	for (int i = 0; i < numOfSteps; i++) {
		dwt_1d_filter(initialValues, step_length, h, g, 2);
		step_length = step_length >> 1;
	}
}

void dwt_haar_2d(uint8_t initialValues[], int rows, int cols, unsigned int numOfSteps) {
	if (numOfSteps > maxNumberOfSteps_2D(rows, cols)) {
		fprintf(stderr, "Error - the number of steps is higher than possible.\n");
		exit(1);
	}

	double h[] = { 0.5, 0.5 };
	double g[] = { 0.5, -0.5 };

	int step_rows = rows;
	int step_cols = cols;

	for (int i = 0; i < numOfSteps; i++) {
		dwt_2d_filter(initialValues, step_rows, step_cols, h, g, 2);

		step_rows = step_rows >> 1;
		step_cols = step_cols >> 1;
	}
}

int isPowerOfTwo(unsigned int x) {
	return ((x != 0) && ((x & (~x + 1)) == x));
}

void convolution(double input[], int n, double h[], double g[], short order) {
	double* output = new double[n];

	short steps = n / order;
	double* buffer = new double[order];

	for (unsigned int step = 0; step < steps; step++) {
		for (unsigned int i = 0; i < order; i++) {
			buffer[i] = input[step * order + i];
		}

		output[step] = h[0] * buffer[0];
		output[step + (steps * order / 2)] = g[0] * buffer[0];

		for (unsigned int i = order - 1; i > 0; i--) {
			output[step] += h[i] * buffer[i];
			output[step + (steps * order / 2)] += g[i] * buffer[i];
		}
	}

	for (int i = 0; i < n; i++)
		input[i] = output[i];
}