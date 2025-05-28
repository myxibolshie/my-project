#ifndef PCA_H
#define PCA_H


#include <iostream>
#include <vector>
#include <string>
#include <opencv2/opencv.hpp>
#include <Eigen/Dense>



using namespace std;
using namespace cv;
using namespace Eigen;


struct MyMatrix {
	vector<vector<double>> matrix;
	int row;
	int  col;
	MyMatrix();
        MyMatrix(int rows, int cols);
        MyMatrix(const vector<vector<double>>& data);
	MyMatrix(initializer_list<initializer_list<double>> init);
	

	int size() const;
	int cols() const;
	

	MyMatrix transposed();
	MyMatrix operator*(const MyMatrix &other);
	MyMatrix operator/(const int n);
	pair<MyMatrix, MyMatrix>  rotation_method(int n);
	double find_angle(MyMatrix &data, int i, int j);
	vector<double> find_max( MyMatrix& vec);
	MyMatrix covariation();
	MyMatrix standardization();
	double average(vector<double>& vec);
	double standart_deviation(vector<double>& vec);
	MyMatrix& operator=(const MyMatrix& _matrixp);
	void print();
		
	
};

	void main_components(MyMatrix& main_matrix, MyMatrix& eig_values,vector<string>& names);
	pair<vector<string>, MyMatrix>  writeMatrixFromCSV(const string& filename,const int row,const int col);
	void writeMatrixToCSV(const  vector<string>& names ,const MyMatrix& vec, const string& filename);
	void hello_world();
	vector<VectorXd> imageToBlocks(const Mat& image, int BlockSize);
	Mat blocksToImage( const vector<VectorXd>& blocks, int rows, int cols, int blockSize);
	VectorXd mean(X.cols());
	pair<VectorXd, MatrixXd> jacobi_eigensolver(const MatrixXd& A);
	MatrixXd eye(int n);
	
	

#endif //PCA_H


