#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "wavelets.h"
using namespace std;


void dwt_1d_lifting(uint8_t a[], int n, int steps, string name) {
	if (steps > maxNumberOfSteps_1D(n)) {
		fprintf(stderr, "Error - the number of steps is higher than possible.\n");
		exit(1);
	}

	int step_length = n;

	for (int i = 0; i < steps; i++) {
		if (name == "haar")
			haar_1d_lifting(a, step_length);
		if (name == "db2")
			db2_1d_lifting(a, step_length);

		step_length = step_length >> 1;
	}
}

void dwt_2d_lifting(uint8_t a[], int rows, int cols, int steps, string name) {
	if (steps > maxNumberOfSteps_2D(rows, cols)) {
		fprintf(stderr, "Error - the number of steps is higher than possible.\n");
		exit(1);
	}

	int step_rows = rows;
	int step_cols = cols;

	for (int i = 0; i < steps; i++) {
		if (name == "haar")
			haar_2d_lifting(a, step_rows, step_cols);
		if (name == "db2")
			db2_2d_lifting(a, step_rows, step_cols);

		step_rows = step_rows >> 1;
		step_cols = step_cols >> 1;
	}

}

void haar_1d_lifting(uint8_t a[], int n) {
	int i, j;
	int half = n >> 1;
	double* temp = new double[n];

	j = 0;
	for (i = 0; i < n; i += 2) {
		temp[j] = (a[i] + a[i + 1]) / 2;
		temp[j + half] = (a[i] - a[i + 1])/2;

		j++;
	}

	for (i = 0; i < n; i++) {
		a[i] = (uint8_t)temp[i];
	}

	delete[] temp;
}

void haar_2d_lifting(uint8_t a[], int rows, int cols) {
	int i, j;
	int half_rows = rows >> 1;
	int half_cols = cols >> 1;

	uint8_t* temp = new uint8_t[rows];

	for (i = 0; i < rows; i++) {
		haar_1d_lifting(&a[i * cols], cols);
	}

	for (i = 0; i < cols; i++) {
		for (j = 0; j < rows; j++)
			temp[j] = a[i + j * cols];

		haar_1d_lifting(temp, rows);
		for (j = 0; j < rows; j++)
			 a[i + j * cols] = temp[j];
	}

	delete[] temp;
}

void db2_1d_lifting(uint8_t a[], int n) {
	if (n < 4) {
		fprintf(stderr, "Error - the size of array must be at least 4.\n");
		exit(1);
	}

	int i, j;
	int half = n >> 1;
	double* temp = new double[n];

	j = 0;
	double x1c, x1prev, x1next, x2c;
	double y1c, y1next;
	for (i = 0; i < n; i += 2) {
		if (i == 0) {
			x1c = a[i] + 1.7320508 * a[i + 1];
			x1prev = 0;
			y1c = a[i + 1] - 0.4330127 * x1c + 0.066987 * x1prev;
		}
		if (i < n - 2) {
			x1next = a[i + 2] + 1.7320508 * a[i + 3];
			y1next = a[i + 3] - 0.4330127 * x1next + 0.066987 * x1c;
		}
		else {
			x1next = 0;
			y1next = 0;
		}
		x2c = x1c - y1next;

		temp[j] = 0.366025 * x2c;
		temp[j + half] = 1.366025 * y1c;

		if (i < n - 2) {
			x1prev = x1c;
			x1c = x1next;
			y1c = y1next;
			j++;
		}
	}

	for (i = 0; i < n; i++) {
		a[i] = (uint8_t)temp[i];
	}
}

void db2_2d_lifting(uint8_t a[], int rows, int cols) {
	int i, j;
	int half_rows = rows >> 1;
	int half_cols = cols >> 1;

	uint8_t* temp = new uint8_t[rows];

	for (i = 0; i < rows; i++) {
		db2_1d_lifting(&a[i * cols], cols);
	}

	for (i = 0; i < cols; i++) {
		for (j = 0; j < rows; j++)
			temp[j] = a[i + j * cols];

		db2_1d_lifting(temp, rows);
		for (j = 0; j < rows; j++)
			a[i + j * cols] = temp[j];
	}

	delete[] temp;
}

int maxNumberOfSteps_1D(int n) {
	int powerOfTwo = 1;
	int steps = 0;
	while (powerOfTwo < n) {
		powerOfTwo = powerOfTwo * 2;
		steps++;
	}

	return steps-1;
}

int maxNumberOfSteps_2D(int numOfRows, int numOfCols) {
	int maxRows = maxNumberOfSteps_1D(numOfRows);
	int maxCols = maxNumberOfSteps_1D(numOfCols);

	return min(maxRows, maxCols);
}