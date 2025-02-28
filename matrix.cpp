#include <iostream>
#include <stdexcept>
#include <vector>
#include <windows.h>


class Matrix {
private:
    int rows_, cols_;
    std::vector<std::vector<double>> matrix_;

public:
    Matrix() : rows_(3), cols_(3), matrix_(3, std::vector<double>(3, 0)) {}
    Matrix(int rows, int cols) : rows_(rows), cols_(cols), matrix_(rows, std::vector<double>(cols, 0)) {}
    Matrix(const Matrix& other) : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {}
    Matrix(Matrix&& other) noexcept : rows_(other.rows_), cols_(other.cols_), matrix_(std::move(other.matrix_)) {
        other.rows_ = 0;
        other.cols_ = 0;
    }
    ~Matrix() = default;


    int getRows() const { return rows_; }
    int getCols() const { return cols_; }

    void print(){
        std::cout<<"\n";
        for(auto & i : matrix_){
            for(double j : i){
             std::cout<< j<<" ";
            }
            std::cout<<"\n";
        }

    }


    bool EqMatrix(const Matrix& other) const {
        if (rows_ != other.rows_ || cols_ != other.cols_) return false;
        for (int i = 0; i < rows_; i++) {
            for (int j = 0; j < cols_; j++) {
                if (matrix_[i][j] != other.matrix_[i][j]) return false;
            }
        }
        return true;
    }

    void SumMatrix(const Matrix& other) {
        if (rows_ != other.rows_ || cols_ != other.cols_) {
            throw std::invalid_argument("Матрицы разного размера!!!!!!!");
        }
        for (int i = 0; i < rows_; ++i) {
            for (int j = 0; j < cols_; ++j) {
                matrix_[i][j] += other.matrix_[i][j];
            }
        }
    }

    void SubMatrix(const Matrix& other) {
        if (rows_ != other.rows_ || cols_ != other.cols_) {
            throw std::invalid_argument("Матрицы разного размера!!!!!!!");
        }
        for (int i = 0; i < rows_; ++i) {
            for (int j = 0; j < cols_; ++j) {
                matrix_[i][j] -= other.matrix_[i][j];
            }
        }
    }

    void MulNumber(const double num) {
        for (int i = 0; i < rows_; ++i) {
            for (int j = 0; j < cols_; ++j) {
                matrix_[i][j] *= num;
            }
        }
    }

    void MulMatrix(const Matrix& other) {
        if (cols_ != other.rows_) {
            throw std::invalid_argument("");
        }
        Matrix result(rows_, other.cols_);
        for (int i = 0; i < rows_; ++i) {
            for (int j = 0; j < other.cols_; ++j) {
                for (int k = 0; k < cols_; ++k) {
                    result.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
                }
            }
        }
        *this = result;
    }

    Matrix Transpose() const {
        Matrix result(cols_, rows_);
        for (int i = 0; i < rows_; ++i) {
            for (int j = 0; j < cols_; ++j) {
                result.matrix_[j][i] = matrix_[i][j];
            }
        }
        return result;
    }

    double Determinant() const {
        if (rows_ != cols_) {
            throw std::invalid_argument("Матрица не квадратная!");
        }

        // Копируем матрицу, чтобы не изменять оригинальную
        std::vector<std::vector<double>> mat = matrix_;
        double det = 1.0;

        for (int i = 0; i < rows_; ++i) {
            // Поиск ведущего элемента (pivot)
            int pivot = i;
            for (int j = i + 1; j < rows_; ++j) {
                if (std::abs(mat[j][i]) > std::abs(mat[pivot][i])) {
                    pivot = j;
                }
            }

            if (pivot != i) {
                std::swap(mat[i], mat[pivot]);
                det *= -1; // Меняем знак определителя при перестановке строк
            }

            if (mat[i][i] == 0) {
                return 0; // Если ведущий элемент нулевой, определитель равен 0
            }

            det *= mat[i][i]; // Умножаем определитель на ведущий элемент

            for (int j = i + 1; j < rows_; ++j) {
                double factor = mat[j][i] / mat[i][i];
                for (int k = i + 1; k < rows_; ++k) {
                    mat[j][k] -= factor * mat[i][k];
                }
            }
        }

        return det;
    }

    Matrix CalcComplements() const {
        if (rows_ != cols_) {
            throw std::invalid_argument("Матрица не квадратная!!!!");
        }
        Matrix result(rows_, cols_);
        for (int i = 0; i < rows_; ++i) {
            for (int j = 0; j < cols_; ++j) {
                Matrix subMat(rows_ - 1, cols_ - 1);
                int subRow = 0;
                for (int k = 0; k < rows_; ++k) {
                    if (k == i) continue;
                    int subCol = 0;
                    for (int l = 0; l < cols_; ++l) {
                        if (l == j) continue;
                        subMat(subRow, subCol) = matrix_[k][l];
                        subCol++;
                    }
                    subRow++;
                }
                double det = subMat.Determinant();
                result(i, j) = ((i + j) % 2 == 0 ? 1 : -1) * det;
            }
        }
        return result;
    }

    Matrix InverseMatrix() const {
        if (rows_ != cols_) {
            throw std::invalid_argument("Матрица не квадратная!");
        }

        double det = Determinant();
        if (det == 0) {
            throw std::invalid_argument("Определитель равен нулю, обратной матрицы не существует!");
        }

        // Вычисляем матрицу алгебраических дополнений
        Matrix complements = CalcComplements();
        // Транспонируем её (получаем присоединённую матрицу)
        Matrix adjugate = complements.Transpose();
        // Делим каждый элемент на определитель
        adjugate.MulNumber(1.0 / det);

        return adjugate;
    }

    Matrix operator+(const Matrix& other) const {
        Matrix result(*this);
        result.SumMatrix(other);
        return result;
    }

    Matrix operator-(const Matrix& other) const {
        Matrix result(*this);
        result.SubMatrix(other);
        return result;
    }

    Matrix operator*(const Matrix& other) const {
        Matrix result(*this);
        result.MulMatrix(other);
        return result;
    }

    Matrix operator*(const double num) const {
        Matrix result(*this);
        result.MulNumber(num);
        return result;
    }

    bool operator==(const Matrix& other) const {
        return EqMatrix(other);
    }

    Matrix& operator=(const Matrix& other) {
        if (this != &other) {
            rows_ = other.rows_;
            cols_ = other.cols_;
            matrix_ = other.matrix_;
        }
        return *this;
    }

    Matrix& operator+=(const Matrix& other) {
        SumMatrix(other);
        return *this;
    }

    Matrix& operator-=(const Matrix& other) {
        SubMatrix(other);
        return *this;
    }

    Matrix& operator*=(const Matrix& other) {
        MulMatrix(other);
        return *this;
    }

    Matrix& operator*=(const double num) {
        MulNumber(num);
        return *this;
    }

    double& operator()(int i, int j) {
        if (i < 0 || i >= rows_ || j < 0 || j >= cols_) {
            throw std::out_of_range("Некорректный индекс!!!!");
        }
        return matrix_[i][j];
    }

    const double& operator()(int i, int j) const {
        if (i < 0 || i >= rows_ || j < 0 || j >= cols_) {
            throw std::out_of_range("Некорректный индекс!!!!");
        }
        return matrix_[i][j];
    }
};

int main() {
    SetConsoleOutputCP(CP_UTF8);

    Matrix a(3, 3);
    a(0, 0) = 1; a(0, 1) = 2; a(0, 2) = 5;
    a(1, 0) = 3; a(1, 1) = 4; a(1, 2) = 7;
    a(2, 0) = 2; a(2, 1) = 6; a(2, 2) = 4;

    Matrix b(2, 2);
    b(0, 0) = 5; b(0, 1) = 6;
    b(1, 0) = 7; b(1, 1) = 8;

    Matrix c(2, 2);
    b(0, 0) = 1; b(0, 1) = 3;
    b(1, 0) = 8; b(1, 1) = 8;

    Matrix m = c + b;

    m.print();
    a.print();



    std::cout<<"\n"<<a.Determinant()<<std::endl;
    std::cout<<b.Determinant()<<std::endl;
    a.CalcComplements().print();
    a.InverseMatrix().print();

    a = a.Transpose();
    a.print();


    return 0;
}