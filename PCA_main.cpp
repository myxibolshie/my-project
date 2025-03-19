







#include "PCA.h"





int main() {
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
   // Matrix data_T;
    pair<vector<string>, Matrix> Pair = writeMatrixFromCSV(name_of_file, number_of_parametres, parametres);
    name = Pair.first;
    data = Pair.second;
//    name.print();
 //  data.print();
    Matrix data_T = data.transposed();
    Matrix data_TR = data.transposed();
//    data_T.print();
   Matrix data_T_standardization =  data_T.standardization();
//    data_T_standardization.print();
    Matrix cov_data = data_T_standardization.covariation();
//    cov_data.print();
    pair<Matrix, Matrix> EIG = cov_data.rotation_method(64);
    Matrix eigenvalues = EIG.first;
    Matrix eigenvectors = EIG.second;
  //  eigenvalues.print();
 //   cout << endl;
 //   eigenvectors.print();
    Matrix eigenvectors_T = eigenvectors.transposed();
    Matrix final_matrix = eigenvectors_T * data_TR;
   final_matrix.print();
  main_components(final_matrix, eigenvalues, name);
  
 Matrix  final_matrix_tr = final_matrix.transposed();
    writeMatrixToCSV(name, final_matrix_tr, "output1.csv");

return 0;
}
