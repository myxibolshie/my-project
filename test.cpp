#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace Eigen;
using namespace cv;

// Функция для разбиения изображения на блоки с проверкой границ
vector<VectorXd> imageToBlocks(const Mat& image, int blockSize) {
    vector<VectorXd> blocks;
    for (int i = 0; i < image.rows; i += blockSize) {
        for (int j = 0; j < image.cols; j += blockSize) {
            // Рассчитываем актуальные размеры блока
            int actualHeight = min(blockSize, image.rows - i);
            int actualWidth = min(blockSize, image.cols - j);
            
            Mat block = image(Rect(j, i, actualWidth, actualHeight)).clone();
            
            // Если блок меньше нужного размера - дополняем нулями
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

// Функция для сборки изображения с обработкой границ
Mat blocksToImage(const vector<VectorXd>& blocks, int rows, int cols, int blockSize) {
    Mat image(rows, cols, CV_64F, Scalar(0));
    int index = 0;
    for (int i = 0; i < rows; i += blockSize) {
        for (int j = 0; j < cols; j += blockSize) {
            int actualHeight = min(blockSize, rows - i);
            int actualWidth = min(blockSize, cols - j);
            
            // Проверка соответствия размера данных
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

int main() {
    // Загрузка и предобработка изображения
    Mat image = imread("chisto (1).jpg", IMREAD_GRAYSCALE);
    if (image.empty()) {
        cerr << "Ошибка загрузки изображения!" << endl;
        return 1;
    }
    
    const int blockSize = 16;    // Размер блока
    const int k = 32;           // Уменьшенное число компонент для лучшего подавления шума
    
    // Конвертация в double и нормализация
    Mat imageDouble;
    image.convertTo(imageDouble, CV_64F);
    normalize(imageDouble, imageDouble, 0, 1, NORM_MINMAX);

    // Разбиение на блоки
    vector<VectorXd> blocks = imageToBlocks(imageDouble, blockSize);
	cout << blocks << endl;
    // Создание матрицы данных
/*    MatrixXd X(blocks.size(), blockSize*blockSize);
    for (size_t i = 0; i < blocks.size(); ++i) {
        X.row(i) = blocks[i];
    }

    // Центрирование данных
    VectorXd mean = X.colwise().mean();
    X.rowwise() -= mean.transpose();

    // Вычисление ковариационной матрицы
    MatrixXd cov = (X.adjoint() * X) / (X.rows() - 1);

    // Вычисление собственных значений/векторов
    SelfAdjointEigenSolver<MatrixXd> solver(cov);
    
    // Сортировка компонент по убыванию собственных значений
    MatrixXd eigenvectors = solver.eigenvectors();
    VectorXd eigenvalues = solver.eigenvalues();
    
    // Инвертируем порядок для сортировки по убыванию
    eigenvectors = eigenvectors.rowwise().reverse().eval();
    eigenvalues = eigenvalues.reverse().eval();

    // Выбор главных компонент
    MatrixXd V_k = eigenvectors.leftCols(k);

    // Проекция и восстановление
    MatrixXd Y = X * V_k;
    MatrixXd X_reconstructed = Y * V_k.transpose();
    X_reconstructed.rowwise() += mean.transpose();

    // Сборка изображения
    vector<VectorXd> reconstructedBlocks;
    for (int i = 0; i < X_reconstructed.rows(); ++i) {
        reconstructedBlocks.push_back(X_reconstructed.row(i));
    }
    
    Mat reconstructedImage = blocksToImage(reconstructedBlocks, image.rows, image.cols, blockSize);
    
    // Постобработка и сохранение
    normalize(reconstructedImage, reconstructedImage, 0, 255, NORM_MINMAX);
    reconstructedImage.convertTo(reconstructedImage, CV_8U);
    imwrite("denoised_image.jpg", reconstructedImage);

    cout << "Обработка завершена. Результат: denoised_image.jpg" << endl; */
    
    return 0;
}
