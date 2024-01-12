#include <iostream>
#include <cmath>
#include <functional>
#include <fstream>

using namespace std;

double G = 6.674e-11;

double calculateIntegral(function<double(double)> *f1,function<double(double)> *f2, double begin, double end){
    double points[2] = {-1/sqrt(3),1/ sqrt(3)};
    double weight = 1.0;
    double scaling_factor = (end-begin)/2.0;
    double average_factor = (begin+end)/2.0;

    double sum = 0;
    for(double point : points){
        sum += weight*(*f1)(scaling_factor*point+average_factor)*(*f2)(scaling_factor*point+average_factor);
    }
    return sum*scaling_factor;
}

void createFunctions(vector<function<double(double)>> *functions, int n, double h){
    for(int i = 0; i < n; i++) {
        function<double(double)> e = [h, i](double x) -> double {
            if(x > h*(i-1) and x <= h*i) return (x/h-i)+1;
            else if(x > h*i and x < h*(i+1)) return i-x/h+1;
            else return 0.0;
        };
        functions->push_back(e);
    }
}

void createDerivatives(vector<function<double(double)>> *derivatives,int n, double h){
    for(int i = 0; i < n; i++) {
        function<double(double)> ePrim = [h, i](double x) -> double {
            if(x >= h*(i-1) and x <= h*i)return 1/h;
            else if(x >= h*i and x <= h*(i+1))return -1/h;
            else return 0.0;
        };
        derivatives->push_back(ePrim);
    }
}

double calculateB(vector<function<double(double)>> *derivatives ,int e1, int e2, double h){
    double lower_boundB = max(h * (e1-1),0.0), upper_boundB = min(h * (e1+1),3.0);
    return (-1) * calculateIntegral(&(*derivatives)[e1], &(*derivatives)[e2],lower_boundB,upper_boundB);
}

double calculateL(vector<function<double(double)>> *functions, vector<function<double(double)>> *derivatives, int e, double h){
    double lower_boundL = max(h * (e-1),1.0), upper_boundL = min(h * (e+1),2.0);
    double lower_boundB = max(h * (e-1),0.0), upper_boundB = min(h * (e+1),3.0);
    function<double(double)> one = [](double x) { return 1.0; };
    function<double(double)> shiftD = [](double x) -> double {return - 1.0 / 3.0;};
    return ((calculateIntegral(&(*functions)[e],&one,lower_boundL,upper_boundL)
    + calculateIntegral(&(*derivatives)[e],&shiftD, lower_boundB, upper_boundB)))
    * 4.0*M_PI*G;
}

void eliminate(vector<vector<double>>& matrix, vector<double>& resultVector, int k) {
    int n = matrix.size();
    for (int i = k + 1; i < n; ++i) {
        double factor = matrix[i][k] / matrix[k][k];
        resultVector[i] -= factor * resultVector[k];
        for (int j = k; j < n; ++j) {
            matrix[i][j] -= factor * matrix[k][j];
        }
    }
}

vector<double> solveUpperTriangular(const vector<vector<double>>& matrix, const vector<double>& resultVector) {
    int n = matrix.size();
    vector<double> x(n, 0.0);

    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) {
            sum += matrix[i][j] * x[j];
        }
        x[i] = (resultVector[i] - sum) / matrix[i][i];
    }

    return x;
}

void gaussianElimination(vector<vector<double>>& matrix, vector<double>& resultVector) {
    int n = matrix.size();
    for (int k = 0; k < n - 1; ++k) {
        eliminate(matrix, resultVector, k);
    }
}

int main() {
    int n, begin=0, end=3;
    cin >> n;
    double h = (end - begin) / (double) n;
    vector<vector<double>> matrix(n-1);
    for(auto& row : matrix){
        row.resize(n-1, 0);
    }
    vector<double> resultVector(n,0);

    vector<function<double(double)>> functions;
    vector<function<double(double)>> derivatives;
    createDerivatives(&derivatives,n,h);
    createFunctions(&functions,n,h);

    for(int i=1;i<n;++i){
        for(int j=1;j<n;++j){
            if(i==j)matrix[i-1][j-1]=calculateB(&derivatives,i,j,h);
            else if(i+1==j)matrix[i-1][j-1]=matrix[j-1][i-1]=calculateB(&derivatives,i,j,h);
        }
        resultVector[i] = calculateL(&functions,&derivatives,i,h);
    }

    gaussianElimination(matrix, resultVector);
    vector<double> solution = solveUpperTriangular(matrix, resultVector);
    solution.push_back(0);

    ofstream dataFile("Gravitational_Potential.txt");

    double x = 0.0;
    for (int i = 0; i <= n; i++) {
        double y = solution[i];
        dataFile << x << " " << y + (5.0 - x / 3.0) << endl;
        x+=3.0/(double)n;
    }

    dataFile.close();

    system("gnuplot -persist -e \"plot 'Gravitational_Potential.txt' with lines\"");

    return 0;
}
