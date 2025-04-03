#ifndef PCA_H
#define PCA_H


#include <iostream>
#include <vector>
#include <string>
#include <opencv2/opencv.hpp>



using namespace std;
using namespace cv;


struct Matrix {
	vector<vector<double>> matrix;
	size_t row;
	size_t col;
	Matrix();
        Matrix(size_t rows, size_t cols);
        Matrix(const vector<vector<double>>& data);
	Matrix(initializer_list<initializer_list<double>> init);
	

	size_t size() const;
	size_t cols() const;
	

	Matrix transposed();
	Matrix operator*(const Matrix &other);
	pair<Matrix, Matrix>  rotation_method(int n);
	double find_angle(Matrix &data, int i, int j);
	vector<double> find_max( Matrix& vec);
	Matrix covariation();
	Matrix standardization();
	double average(vector<double>& vec);
	double standart_deviation(vector<double>& vec);
	Matrix& operator=(const Matrix& _matrixp);
	void print();
	Matrix center(Matrix& Mat, int a);

	
};

	void main_components(Matrix& main_matrix, Matrix& eig_values,vector<string>& names);
	pair<vector<string>, Matrix>  writeMatrixFromCSV(const string& filename,const int row,const int col);
	void writeMatrixToCSV(const  vector<string>& names ,const Matrix& vec, const string& filename);
	void hello_world();
	vector<double> blockToVector(Mat& block);
	Matrix imageToBlocks(Mat& image, int BlockSize);
	Mat blocksToImage( Matrix& blocks, int rows, int cols, int blockSize);
	
	

#endif //PCA_H


