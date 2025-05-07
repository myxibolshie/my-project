#include "PCA.h"
#include <fstream>
#include <sstream>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <typeinfo>



MyMatrix::MyMatrix(int rows, int cols) {
		row = rows;
		col = cols;
		matrix.resize(rows, vector<double>(cols, 0.0));
}

MyMatrix::MyMatrix(const vector<vector<double>>& data) {
		matrix = data;
		row = matrix.size();
		col = matrix[0].size();
}

MyMatrix::MyMatrix(initializer_list<initializer_list<double>> init) {
        	for (auto row : init) {
            		matrix.push_back(vector<double>(row));
        	}
		row = matrix.size();
		col = matrix[0].size();
}

int MyMatrix::size() const {
        	return matrix.size();
}

int MyMatrix::cols() const {
        	return matrix[0].size();
}

MyMatrix MyMatrix::transposed()  {
	        MyMatrix transposed(col, row);
		for (int i = 0; i < row; i++) {
		           for(int j = 0; j < col; j++) {
				   transposed.matrix[j][i] = matrix[i][j];
			   }
		}
		return  transposed;
}



MyMatrix MyMatrix::operator*(const MyMatrix &other ) {
		MyMatrix new_vec(row, other.cols());
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

MyMatrix MyMatrix::operator/(const int  n) {
	MyMatrix new_vec(row, col);
	for(size_t i = 0; i < row; i) {
		for(size_t j = 0; j < col; j++) {
			new_vec.matrix[i][j] = matrix[i][j] / n;
		}
	}
	return new_vec;
}

pair<MyMatrix, MyMatrix>  MyMatrix::rotation_method(int param) {
		MyMatrix eig_vectors(param, param);
		for (size_t i = 0; i < param; i++) {
			eig_vectors.matrix[i][i] = 1;
		}
		MyMatrix vec = *this;

		double epsilon = 0.001;
		while(true) {
			auto max_index = find_max(vec);
			int i = max_index[0];
			int j = max_index[1];


	        	if (abs(vec.matrix[i][j]) < epsilon) break;
	        	MyMatrix R(vec.size(), vec.size());
	        	for (size_t k = 0; k < vec.size(); k++) {
			 	R.matrix[k][k] = 1;
			}
                	double c = cos(find_angle(vec, i, j));
			double s = sin(find_angle(vec, i, j));
			R.matrix[i][i] = c;
			R.matrix[j][j] = c;
			R.matrix[i][j] = -s;
			R.matrix[j][i] = s;
			MyMatrix R_tr = R.transposed();
			MyMatrix C = R_tr * vec;
			vec =  C * R;
			eig_vectors = eig_vectors * R;

		}
		return {vec, eig_vectors };
}

double MyMatrix::find_angle(MyMatrix &data, int i, int j) {
		double angle;
        	if (data.matrix[i][i] == data.matrix[j][j]) {
			angle = M_PI / 4;
		}
		else {
			angle = 0.5 * atan2( 2*data.matrix[i][j], data.matrix[i][i] - data.matrix[j][j]);
		}
		return angle;
}

vector<double> MyMatrix::find_max( MyMatrix& vec) {
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

MyMatrix MyMatrix::covariation() {
		MyMatrix new_vec(row, row);
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
MyMatrix MyMatrix::standardization() {
		MyMatrix new_vec(row, col);
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


double MyMatrix::average( vector<double>& vec ) {
		double average = 0;
		double i = 0;
		for (auto vec_element: vec) {
			average += vec_element;
			i += 1;
		}
		average = average / i;
        	return average;
}



double MyMatrix::standart_deviation( vector<double>& vec) {
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

MyMatrix& MyMatrix::operator=(const MyMatrix& _matrixp) {
	    if (this != &_matrixp) {
		    for (int i = 0; i < _matrixp.size(); i++) {
			    for (int j = 0; j < _matrixp.cols(); j++) {
				    matrix[i][j] = _matrixp.matrix[i][j];
			    }
		   }

	}
	return *this;
}

void MyMatrix::print() {
		for(size_t i = 0; i < matrix.size(); i++) {
			for(size_t j = 0; j < matrix[0].size(); j ++) {
				cout << matrix[i][j] << " ";
			}
			cout << endl;
		}
}

pair<vector<string>, MyMatrix>  writeMatrixFromCSV(const string& filename,const int row,const int col) {
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
		MyMatrix data(row, col);
		for (size_t i = 0; i < row; i++) {
			for (size_t j = 0; j < col; j++) {
				string line;
				cout << line << endl;
				getline(file,line,',');
				data.matrix[i][j] = stod(line);
			}
		}
		return {names, data};
}

void writeMatrixToCSV(const  vector<string>& names ,const MyMatrix& vec, const string& filename) {
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
	MyMatrix  S = {{1,2,3,4,45,45,89,12} ,{23,45,67,89,54,23,89,43}, {45,56,67,78,12,18,65,87}};
	vector<string> s = {"Length", "Width", "Height"};
	MyMatrix S_tr = S.transposed();
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


void main_components(MyMatrix& main_matrix, MyMatrix& eig_values, vector<string>& names) {
    double sum = 0;
    vector<double> eig_value;

    
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



vector<VectorXd> imageToBlocks(const Mat& image, int blockSize) {
    vector<VectorXd> blocks;
    for (int i = 0; i < image.rows; i += blockSize) {
        for (int j = 0; j < image.cols; j += blockSize) {
            // Рассчитываем актуальные размеры блока
            int actualHeight = min(blockSize, image.rows - i);
            int actualWidth = min(blockSize, image.cols - j);
            
            Mat block = image(Rect(j, i, actualWidth, actualHeight)).clone();
            
                        if (actualHeight < blockSize || actualWidth < blockSize) {
                Mat paddedBlock = Mat::zeros(blockSize, blockSize, CV_64F);
                block.copyTo(paddedBlock(Rect(0, 0, actualWidth, actualHeight)));
                block = paddedBlock;
            }
            
            VectorXd vec(Map<VectorXd>(block.ptr<double>(), blockSize * blockSize));
            blocks.push_back(vec);
        }
    }
    return blocks;
}

Mat blocksToImage(const vector<VectorXd>& blocks, int rows, int cols, int blockSize) {
    Mat image(rows, cols, CV_64F, Scalar(0));
    int index = 0;
    for (int i = 0; i < rows; i += blockSize) {
        for (int j = 0; j < cols; j += blockSize) {
            int actualHeight = min(blockSize, rows - i);
            int actualWidth = min(blockSize, cols - j);
            
            
            if (blocks[index].size() != blockSize*blockSize) {
                cerr << "Ошибка размера блока!" << endl;
                exit(1);
            }
            
            Mat fullBlock(blockSize, blockSize, CV_64F, const_cast<double*>(blocks[index].data()));
            Mat targetROI = image(Rect(j, i, actualWidth, actualHeight));
            fullBlock(Rect(0, 0, actualWidth, actualHeight)).copyTo(targetROI);
            index++;
        }
    }
    return image;
}










