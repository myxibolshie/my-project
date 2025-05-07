#include "PCA.h"



int main() {
    string the_purpose;
    cout << "The purpose(photo or data): ";
    cin >> the_purpose;
    if ( the_purpose == "data") {
    	hello_world();
    	vector<string> name;
    	int number_of_parametres;
    	int parametres;
    	string name_of_file;
    	cout << "name_of_file: " ;
    	cin >> name_of_file;
  
    	cout << "The number of parametres: ";
    	cin >> number_of_parametres;
    	cout << "parametres: ";
    	cin >> parametres;
    	MyMatrix data(number_of_parametres, parametres);    
   
    	pair<vector<string>, MyMatrix> Pair = writeMatrixFromCSV(name_of_file, number_of_parametres, parametres);
    	name = Pair.first;
    	data = Pair.second;


    	MyMatrix data_T = data.transposed();
    	MyMatrix data_TR = data.transposed();

   	MyMatrix data_T_standardization =  data_T.standardization();

    	MyMatrix cov_data = data_T_standardization.covariation();

    	pair<MyMatrix, MyMatrix> EIG = cov_data.rotation_method(cov_data.size());
    	MyMatrix eigenvalues = EIG.first;
    	MyMatrix eigenvectors = EIG.second;
	

    	MyMatrix eigenvectors_T = eigenvectors.transposed();
   	MyMatrix final_matrix = eigenvectors_T * data_TR;
   	
  	main_components(final_matrix, eigenvalues, name);
  
 	MyMatrix  final_matrix_tr = final_matrix.transposed();
    	writeMatrixToCSV(name, final_matrix_tr, "output.csv");
        }
    else {
	   
    Mat image = imread("image.jpg", IMREAD_GRAYSCALE);
    if (image.empty()) {
        cerr << "sdОшибка загрузки изображения!" << endl;
        return 1;
    }
    
    const int blockSize = 16;    
    const int k = 32;           
    
    
    Mat imageDouble;
    image.convertTo(imageDouble, CV_64F);
    normalize(imageDouble, imageDouble, 0, 1, NORM_MINMAX);

 
    vector<VectorXd> blocks = imageToBlocks(imageDouble, blockSize);
	
    
    MatrixXd X(blocks.size(), blockSize*blockSize);
    for (size_t i = 0; i < blocks.size(); ++i) {
        X.row(i) = blocks[i];
    }


    
    VectorXd mean = X.colwise().mean();
    X.rowwise() -= mean.transpose();


   
    MatrixXd cov = (X.adjoint() * X) / (X.rows() - 1);
   
    
    SelfAdjointEigenSolver<MatrixXd> solver(cov);
    
    
    MatrixXd eigenvectors = solver.eigenvectors();
    
    VectorXd eigenvalues = solver.eigenvalues();
    
    
    eigenvectors = eigenvectors.rowwise().reverse().eval();
    eigenvalues = eigenvalues.reverse().eval();

  
    MatrixXd V_k = eigenvectors.leftCols(k);

    
    MatrixXd Y = X * V_k;
    MatrixXd X_reconstructed = Y * V_k.transpose();
    X_reconstructed.rowwise() += mean.transpose();

   
    vector<VectorXd> reconstructedBlocks;
    for (int i = 0; i < X_reconstructed.rows(); ++i) {
        reconstructedBlocks.push_back(X_reconstructed.row(i));
    }
    
    Mat reconstructedImage = blocksToImage(reconstructedBlocks, image.rows, image.cols, blockSize);
    
  
    normalize(reconstructedImage, reconstructedImage, 0, 255, NORM_MINMAX);
    reconstructedImage.convertTo(reconstructedImage, CV_8U);
    imwrite("denoised_image.jpg", reconstructedImage);

    cout << "Обработка завершена. Результат: denoised_image.jpg" << endl; 
	}    

	return 0;
}
