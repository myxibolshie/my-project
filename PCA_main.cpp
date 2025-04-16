







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
    	Matrix data(number_of_parametres, parametres);    
   
    	pair<vector<string>, Matrix> Pair = writeMatrixFromCSV(name_of_file, number_of_parametres, parametres);
    	name = Pair.first;
    	data = Pair.second;


    	Matrix data_T = data.transposed();
    	Matrix data_TR = data.transposed();

   	Matrix data_T_standardization =  data_T.standardization();

    	Matrix cov_data = data_T_standardization.covariation();

    	pair<Matrix, Matrix> EIG = cov_data.rotation_method(cov_data.size());
    	Matrix eigenvalues = EIG.first;
    	Matrix eigenvectors = EIG.second;



    	Matrix eigenvectors_T = eigenvectors.transposed();
   	Matrix final_matrix = eigenvectors_T * data_TR;
   	final_matrix.print();
  	main_components(final_matrix, eigenvalues, name);
  
 	Matrix  final_matrix_tr = final_matrix.transposed();
    	writeMatrixToCSV(name, final_matrix_tr, "output1.csv");
        }
    else {
	Mat image = imread("chisto (1).jpg",IMREAD_GRAYSCALE);
	const int blockSize = 16;
	const int k = 64;

	Mat imageDouble;
	image.convertTo(imageDouble,CV_64F);
	normalize(imageDouble, imageDouble, 0, 1, NORM_MINMAX);
	
	
	Matrix blocks = imageToBlocks(imageDouble, blockSize);
	Matrix blocks_cnt = blocks.center(blocks, 0);
		
//	b.print();
        

	Matrix blocks_cov = blocks_cnt.covariation();
//	blocks_cov.print();
//	cout << blocks_cov.size() << endl;
//	cout << blocks_cov.cols()<< endl;	
	pair<Matrix, Matrix> EIG = blocks_cov.rotation_method(blocks_cov.size());
	Matrix eigenvalues = EIG.first;
	Matrix eigenvectors = EIG.second;
	Matrix eigenvectors_new = eigenvectors.leftCols(k);
//	eigenvalues.print();
	Matrix Y = blocks_cnt * eigenvectors_new;
        Matrix X_final = Y * eigenvectors_new.transposed();	
	X_final = X_final.center(X_final,1);
	



//	main_components(final_matrix, eigenvalues, name);
	Mat construktImage = blocksToImage(X_final, image.rows, image.cols, blockSize);
	normalize(construktImage,construktImage,0,255,NORM_MINMAX);
	construktImage.convertTo(construktImage,CV_8U);
	imwrite("reconstruct_image_2.jpg",construktImage);
	cout << "The method is done!" << endl;

	
	
	
	
		


	}    
return 0;
}
