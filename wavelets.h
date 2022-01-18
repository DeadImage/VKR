#include <vector>
#include <string>
using namespace std;

void dwt_1d_lifting(uint8_t a[], int n, int steps, string name);
void dwt_2d_lifting(uint8_t a[], int rows, int cols, int steps, string name);

void haar_1d_lifting(uint8_t a[], int n);
void haar_2d_lifting(uint8_t a[], int rows, int cols);

void db2_1d_lifting(uint8_t a[], int n);
void db2_2d_lifting(uint8_t a[], int rows, int cols);

int maxNumberOfSteps_1D(int n);
int maxNumberOfSteps_2D(int numOfRows, int numOfCols);

void dwt_1d_filter(uint8_t initialValues[], int numOfValues, double h[], double g[], short order);
void dwt_2d_filter(uint8_t initialValues[], int rows, int cols, double h[], double g[], short order);

void dwt_haar_1d(uint8_t initialValues[], int n, unsigned int numOfSteps);
void dwt_haar_2d(uint8_t initialValues[], int rows, int cols, unsigned int numOfSteps);

int isPowerOfTwo(unsigned int x);
void convolution(double input[], int n, double h[], double g[], short order);