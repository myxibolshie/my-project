#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace Eigen;
using namespace cv;

// ������� ��� ��������� ����������� �� ����� � ��������� ������
vector<VectorXd> imageToBlocks(const Mat& image, int blockSize) {
    vector<VectorXd> blocks;
    for (int i = 0; i < image.rows; i += blockSize) {
        for (int j = 0; j < image.cols; j += blockSize) {
            // ������������ ���������� ������� �����
            int actualHeight = min(blockSize, image.rows - i);
            int actualWidth = min(blockSize, image.cols - j);
            
            Mat block = image(Rect(j, i, actualWidth, actualHeight)).clone();
            
            // ���� ���� ������ ������� ������� - ��������� ������
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

// ������� ��� ������ ����������� � ���������� ������
Mat blocksToImage(const vector<VectorXd>& blocks, int rows, int cols, int blockSize) {
    Mat image(rows, cols, CV_64F, Scalar(0));
    int index = 0;
    for (int i = 0; i < rows; i += blockSize) {
        for (int j = 0; j < cols; j += blockSize) {
            int actualHeight = min(blockSize, rows - i);
            int actualWidth = min(blockSize, cols - j);
            
            // �������� ������������ ������� ������
            if (blocks[index].size() != blockSize*blockSize) {
                cerr << "������ ������� �����!" << endl;
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

int main() {
    // �������� � ������������� �����������
    Mat image = imread("chisto (1).jpg", IMREAD_GRAYSCALE);
    if (image.empty()) {
        cerr << "������ �������� �����������!" << endl;
        return 1;
    }
    
    const int blockSize = 16;    // ������ �����
    const int k = 32;           // ����������� ����� ��������� ��� ������� ���������� ����
    
    // ����������� � double � ������������
    Mat imageDouble;
    image.convertTo(imageDouble, CV_64F);
    normalize(imageDouble, imageDouble, 0, 1, NORM_MINMAX);

    // ��������� �� �����
    vector<VectorXd> blocks = imageToBlocks(imageDouble, blockSize);
	cout << blocks << endl;
    // �������� ������� ������
/*    MatrixXd X(blocks.size(), blockSize*blockSize);
    for (size_t i = 0; i < blocks.size(); ++i) {
        X.row(i) = blocks[i];
    }

    // ������������� ������
    VectorXd mean = X.colwise().mean();
    X.rowwise() -= mean.transpose();

    // ���������� �������������� �������
    MatrixXd cov = (X.adjoint() * X) / (X.rows() - 1);

    // ���������� ����������� ��������/��������
    SelfAdjointEigenSolver<MatrixXd> solver(cov);
    
    // ���������� ��������� �� �������� ����������� ��������
    MatrixXd eigenvectors = solver.eigenvectors();
    VectorXd eigenvalues = solver.eigenvalues();
    
    // ����������� ������� ��� ���������� �� ��������
    eigenvectors = eigenvectors.rowwise().reverse().eval();
    eigenvalues = eigenvalues.reverse().eval();

    // ����� ������� ���������
    MatrixXd V_k = eigenvectors.leftCols(k);

    // �������� � ��������������
    MatrixXd Y = X * V_k;
    MatrixXd X_reconstructed = Y * V_k.transpose();
    X_reconstructed.rowwise() += mean.transpose();

    // ������ �����������
    vector<VectorXd> reconstructedBlocks;
    for (int i = 0; i < X_reconstructed.rows(); ++i) {
        reconstructedBlocks.push_back(X_reconstructed.row(i));
    }
    
    Mat reconstructedImage = blocksToImage(reconstructedBlocks, image.rows, image.cols, blockSize);
    
    // ������������� � ����������
    normalize(reconstructedImage, reconstructedImage, 0, 255, NORM_MINMAX);
    reconstructedImage.convertTo(reconstructedImage, CV_8U);
    imwrite("denoised_image.jpg", reconstructedImage);

    cout << "��������� ���������. ���������: denoised_image.jpg" << endl; */
    
    return 0;
}
