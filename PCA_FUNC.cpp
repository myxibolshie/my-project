#include "PCA.h"
#include <fstream>
#include <sstream>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <typeinfo>



Matrix::Matrix(size_t rows, size_t cols) {
		row = rows;
		col = cols;
		matrix.resize(rows, vector<double>(cols, 0.0));
}

Matrix::Matrix(const vector<vector<double>>& data) {
		matrix = data;
		row = matrix.size();
		col = matrix[0].size();
}

Matrix::Matrix(initializer_list<initializer_list<double>> init) {
        	for (auto row : init) {
            		matrix.push_back(vector<double>(row));
        	}
		row = matrix.size();
		col = matrix[0].size();
}

size_t Matrix::size() const {
        	return matrix.size();
}

size_t Matrix::cols() const {
        	return matrix[0].size();
}

Matrix Matrix::transposed()  {
	        Matrix transposed(col, row);
		for (int i = 0; i < row; i++) {
		           for(int j = 0; j < col; j++) {
				   transposed.matrix[j][i] = matrix[i][j];
			   }
		}
		return  transposed;
}



Matrix Matrix::operator*(const Matrix &other ) {
		Matrix new_vec(row, other.cols());
		for(size_t i = 0; i < row; i++) {
			for(size_t j = 0; j < other.cols(); j++) {
				new_vec.matrix[i][j] = 0;
				for(size_t k = 0; k < col; k++) {
					new_vec.matrix[i][j] +=  matrix[i][k] * other.matrix[k][j];
				}
			}
		}
		return new_vec;
}

pair<Matrix, Matrix>  Matrix::rotation_method(int param) {
		Matrix eig_vectors(param, param);
		for (size_t i = 0; i < param; i++) {
			eig_vectors.matrix[i][i] = 1;
		}
		Matrix vec = *this;

		double epsilon = 0.1;
		while(true) {
			auto max_index = find_max(vec);
			int i = max_index[0];
			int j = max_index[1];


	        	if (abs(vec.matrix[i][j]) < epsilon) break;
	        	Matrix R(vec.size(), vec.size());
	        	for (size_t k = 0; k < vec.size(); k++) {
			 	R.matrix[k][k] = 1;
			}
                	double c = cos(find_angle(vec, i, j));
			double s = sin(find_angle(vec, i, j));
			R.matrix[i][i] = c;
			R.matrix[j][j] = c;
			R.matrix[i][j] = -s;
			R.matrix[j][i] = s;
			Matrix R_tr = R.transposed();
			Matrix C = R_tr * vec;
			vec =  C * R;
			eig_vectors = eig_vectors * R;

		}
		return {vec, eig_vectors };
}

double Matrix::find_angle(Matrix &data, int i, int j) {
		double angle;
        	if (data.matrix[i][i] == data.matrix[j][j]) {
			angle = M_PI / 4;
		}
		else {
			angle = 0.5 * atan2( 2*data.matrix[i][j], data.matrix[i][i] - data.matrix[j][j]);
		}
		return angle;
}

vector<double> Matrix::find_max( Matrix& vec) {
		double max = -10;
		double i_max;
		double j_max;
		for(size_t i = 0; i < vec.size(); i++) {
			for(size_t j = i+1; j < vec.size(); j++) {
				if (abs(vec.matrix[i][j]) > max) {
			       		max = abs(vec.matrix[i][j]);
		               		i_max = i;
		               		j_max = j;
		        	}
	        	}
        	}
        	vector<double> v = {i_max, j_max};
		return v;
}

Matrix Matrix::covariation() {
		Matrix new_vec(row, row);
		vector<double> average_value;
		for(size_t i = 0; i < row; i++) {
			average_value.push_back(average(matrix[i]));
		}
		for(size_t i = 0; i < row; i++) {
			for(size_t j = 0; j < row; j++) {
				double cov = 0 ;
				for(size_t k = 0; k < col; k++) {
					cov += matrix[i][k] * matrix[j][k];
				}
				cov = (cov / col) - (average_value[i] * average_value[j]);
				new_vec.matrix[i][j] = cov;
			}
		}
		return new_vec;
}

Matrix Matrix::standardization() {
		Matrix new_vec(row, col);
		for(size_t i = 0; i < row; i++){
			for(size_t j = 0; j < col; j++) {
				if (abs(matrix[i][j] -  average(matrix[i])) < 1e-10) {
					new_vec.matrix[i][j] = 0;
				}
				else {
					new_vec.matrix[i][j] = (matrix[i][j] - average(matrix[i])) / standart_deviation(matrix[i]);
				}
			}
		}
		return new_vec;
}


double Matrix::average( vector<double>& vec ) {
		double average = 0;
		double i = 0;
		for (auto vec_element: vec) {
			average += vec_element;
			i += 1;
		}
		average = average / i;
        	return average;
}


Matrix Matrix::center(Matrix& Mat, int a = 0){
	Matrix Mat_Tr = Mat.transposed();
	
	for(int i = 0; i < Mat_Tr.size(); i++){
		double num = average(Mat_Tr.matrix[i]);
		for(int j = 0; j < Mat_Tr.cols(); j++){
			if (a==0){
				Mat_Tr.matrix[i][j] -= num;
			}
			else{
				Mat_Tr.matrix[i][j] += num;
} 
		}
	}
	Matrix Mat_Tr_Tr = Mat_Tr.transposed();
	return Mat_Tr_Tr;
}

double Matrix::standart_deviation( vector<double>& vec) {
		double star_dev = 0;
		double averag = average(vec);
		double size = vec.size();
		for (double vec_el: vec) {
			double s;
			s = vec_el - averag;
			star_dev += pow(s, 2);

		}
		star_dev = star_dev / size;
		star_dev = pow(star_dev, 0.5);
		return star_dev;
}

Matrix& Matrix::operator=(const Matrix& _matrixp) {
	    if (this != &_matrixp) {
		    for (int i = 0; i < _matrixp.size(); i++) {
			    for (int j = 0; j < _matrixp.cols(); j++) {
				    matrix[i][j] = _matrixp.matrix[i][j];
			    }
		   }

	}
	return *this;
}

void Matrix::print() {
		for(size_t i = 0; i < matrix.size(); i++) {
			for(size_t j = 0; j < matrix[0].size(); j ++) {
				cout << matrix[i][j] << " ";
			}
			cout << endl;
		}
}

pair<vector<string>, Matrix>  writeMatrixFromCSV(const string& filename,const int row,const int col) {
		ifstream file(filename);
		string line;
		if (!file.is_open()) {
			cerr << "Error" << endl;

		}
		vector<string> names;
		for(size_t i = 0; i < col; i++) {
			string s;
			getline(file,s,',');
			names.push_back(s);
			cout << names[i] << endl;
		}

		int v = 0;
		Matrix data(row, col);
		for (size_t i = 0; i < row; i++) {
			for (size_t j = 0; j < col; j++) {
				string line;
				cout << line << endl;
			//	cout << typeid(line).name() << endl;
				getline(file,line,',');
			//	cout << typeid(line).name() << endl;
			//	cout << typeid(stod(line)).name() << endl;
			//	cout << line << endl;
			//	try {
				data.matrix[i][j] = stod(line);
			//	} catch ( const invalid_argument& e ) {
			//		v += 1;
			//		cout << "Error value " << v <<" " <<  line<<  endl;
			//	}
			}
		}
		return {names, data};
}

void writeMatrixToCSV(const  vector<string>& names ,const Matrix& vec, const string& filename) {
		ofstream file(filename);
		if (!file.is_open()) {
        		cerr << "Ошибка открытия файла!" << endl;
        		return;
    		}

    		for(auto name: names) {
	    		file << name;
	    		if ( name != names.back() ) {
		    		file << ",";
		   	}
    		}
    		file << endl;

    		for (int i = 0; i < vec.size(); i++) {
        		for (int j = 0; j < vec.cols(); j++) {
            			file << vec.matrix[i][j];
            			if (j < vec.cols() - 1) {
                			file << ",";
            			}
        		}
        		file << endl;
    		}

    		file.close();
}

void hello_world() {
	cout << "Hello, you are welcomed by a program that uses Principal Component Analysis" << endl;
	cout << "The introduction was created for the correct use of the program" << endl;
	cout << "The csv file type that is needed for the program to work correctly will now be shown." << endl;
	cout << endl;
	Matrix  S = {{1,2,3,4,45,45,89,12} ,{23,45,67,89,54,23,89,43}, {45,56,67,78,12,18,65,87}};
	vector<string> s = {"Length", "Width", "Height"};
	Matrix S_tr = S.transposed();
	for(auto el: s){
		cout << el << ",";
	}
	cout << endl;
	for(size_t i = 0; i < 8; i++) {
		for(size_t j = 0; j < 3; j++) {
			cout << S_tr.matrix[i][j] << ",";
		}
		cout << endl;
	}
	cout << "As you can see, the first line of the csv file shows the system parameters and under each parameter there is a data set that corresponds to it." << endl;
	cout << "It is also important to note that each value is separated by a comma and each line ends with a comma." << endl;
	cout << "It should be understood that, for example, the first 3 parameters of the system correspond to something. In our example, the first 3 parameters indicate the length, width                  and height of the cube" << endl;

}


void main_components(Matrix& main_matrix, Matrix& eig_values, vector<string>& names) {
    double sum = 0;
    vector<double> eig_value;

    // Суммируем собственные значения и сохраняем их в вектор
    for(size_t i = 0; i < eig_values.size(); i++) {
        sum += eig_values.matrix[i][i];
        eig_value.push_back(eig_values.matrix[i][i]);
    }

    vector<double> dispersion;
    for(size_t i = 0; i < eig_values.size(); i++) {
        dispersion.push_back((eig_value[i] / sum) * 100);
    }

    string word;
    cout << "Do you want to delete some components?(Y/N) ";
    cin >> word;

    if (word == "Y") {
        for(size_t i = 0; i < eig_values.size(); i++) {
            cout << "If you delete component number " << i + 1 << ", you will lose " << dispersion[i] << " percent of the information." << endl;
        }

        cout << "Which component do you want to remove? " << endl;
        int s;
        cout << "If you don't want to delete anything else, write the number 0: ";
        cin >> s;


        while (s != 0) {

            if (s > 0 && s <= static_cast<int>(main_matrix.matrix.size())) {
                main_matrix.matrix.erase(main_matrix.matrix.begin() + (s - 1));
		main_matrix.row -= 1;
		if (names.size() != 0){
                names.erase(names.begin() + (s - 1));
		}
            } 
	    else {
                cout << "Invalid component number. Please enter a valid number." << endl;
            }


            cout << "Which component do you want to remove? (0 to stop): ";
            cin >> s;
        }
    }

}


Matrix imageToBlocks(Mat& image, int BlockSize){
	vector<vector<double>>  blocks;
	for (size_t i=0; i < image.rows; i += BlockSize){
		for (size_t j=0; j < image.cols; j += BlockSize){
			Mat block = image(Rect(j, i, BlockSize, BlockSize)).clone();
			vector<double> v = blockToVector(block);
			blocks.push_back(v);
			



		}
	}
	Matrix blocksMatrix(blocks);
	return blocksMatrix;


}

vector<double>  blockToVector(Mat& block){
	vector<double> v;
	
	if (block.empty() || block.type() != CV_64F){
		cout << "block is empty or type is not CV_64F"<<endl;
		return v;
	}
	
	
	for (int i = 0;i < block.rows; i += 1){
		for (int j = 0; j < block.cols; j +=1){
			v.push_back(block.at<double>(i,j));

		}
	}
	return v;
}



Mat blocksToImage(Matrix& blocks, int rows, int cols, int blockSize) {
    Mat image(rows, cols, CV_64F, Scalar(0));
    int index = 0;
    for (int i = 0; i < rows; i += blockSize) {
        for (int j = 0; j < cols; j += blockSize) {
            int actualHeight = min(blockSize, rows - i);
            int actualWidth = min(blockSize, cols - j);
            
            
            if (blocks.matrix[index].size() != blockSize*blockSize) {
                cerr << "Ошибка размера блока!" << endl;
                exit(1);
            }
            
            Mat fullBlock(blockSize, blockSize, CV_64F, const_cast<double*>(blocks.matrix[index].data()));
            Mat targetROI = image(Rect(j, i, actualWidth, actualHeight));
            fullBlock(Rect(0, 0, actualWidth, actualHeight)).copyTo(targetROI);
            index++;
        }
    }
    return image;
}

Matrix Matrix::leftCols(int k) const { 
    

    Matrix result(size(), k); 

    
    for (int i = 0; i < size(); ++i) {
        for (int j = 0; j < k; ++j) {
            result.matrix[i][j] = (*this).matrix[i][j]; 
        }
    }

    return result;
}









