#include <vector>
using namespace std;

vector<double> dwt_lift_1d(vector<double> &initialValues, unsigned int numOfSteps);
vector<vector<double>> dwt_lift_2d(vector<vector<double>>& initialValues, unsigned int numOfSteps);

vector<double> dwt_1d_filter(vector<double>& initialValues, unsigned int numOfSteps, double h[], double g[], short order);
vector<double> dwt_2d_filter(vector<double>& initialValues, unsigned int numOfSteps, double h[], double g[], short order);

vector<double> dwt_haar_1d(vector<double>& initialValues, unsigned int numOfSteps);
vector<vector<double>> dwt_haar_2d(vector<vector<double>>& initialValues, unsigned int numOfSteps);

int isPowerOfTwo(unsigned int x);
unsigned int maxNumberOfSteps_1D(size_t numOfValues);
unsigned int maxNumberOfSteps_2D(size_t numOfRows, size_t numOfCols);
vector<double> convolution(vector<double>& input, double h[], double g[], short order);
