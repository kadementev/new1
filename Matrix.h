#pragma once

#define MATRIX_EXPORTS 1

#ifdef MATRIX_EXPORTS
#define MATRIX_API __declspec(dllexport)
#else
#define MATRIX_API __declspec(dllimport)
#endif



#include <initializer_list>
#include <iostream>
#include <vector>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <fstream>
#include <string>
using namespace std;

class MATRIX_API Matrix {
public:
    Matrix();
    Matrix(int rows, int columns);
    Matrix(const initializer_list<const initializer_list<double>> &m);
    void random_fill(double min, double max);
    void set(double num, int row, int column);
    Matrix admult(Matrix B);
    double trace();
    double det();
    double norm();
    double max1();
    double sum();
    void copy(Matrix m);
    void get_bin(ifstream& in);
    void write_bin(ofstream& out);
    Matrix inv();
    Matrix minor(int row, int column);
    double angle(Matrix B);
    double dot(Matrix B);
    double get(int row, int column) const;
    int rows() const;
    int columns() const;
    int rang();
    void add_matrix(Matrix m);
    Matrix slice(int start_row, int start_column, int end_row, int end_column);
    virtual ~Matrix();
    Matrix T();
    Matrix operator+(Matrix B);
    Matrix operator+(double B);
    Matrix operator-(Matrix B);
    Matrix operator*(double B);
    Matrix operator*(Matrix B);
    Matrix operator/(double B);
    Matrix operator/(Matrix m);
    friend ostream &operator<<(ostream& os, const Matrix& m);
    friend ifstream &operator<<(Matrix& m, ifstream& in);
    friend ofstream &operator>>(Matrix& m, ofstream& out);
    friend class PCA;
protected:
    int _rows;
    int _columns;
    vector<vector<double>*>* _data = nullptr;

};

ostream& operator<<(ostream& os, const Matrix& m);



Matrix::Matrix() {};


Matrix::Matrix(int rows, int columns) {
    _rows = rows;
    _columns = columns;
    _data = new vector<vector<double>*>;
    for (int i = 0; i < _rows; i++) {
        _data->push_back(new vector<double>);
        for (int j = 0; j < _columns; j++) {
            _data->at(i)->push_back(0);
        }
    }
};
Matrix::Matrix(const initializer_list<const initializer_list<double>>& m) {
    int i = 0;
    int j = 0;
    _rows = 0;
    _columns = m.begin()->size();
    _data = new vector<vector<double>*>;
    for (auto const& list : m) {
        _data->push_back(new vector<double>);
        _rows++;
        for (auto const& element : list) {
            _data->at(i)->push_back(element);
            j++;
        }
        if (j != _columns) {
            cerr << "Definition error\nWrong matrix row size" << endl;
            std::exit(1);
        }
        j = 0;
        i++;
    }
}
ostream& operator<<(ostream& os, const Matrix& m) {
    if (m._data == nullptr) {
        cerr << "Memory Error\nMatrix not specified" << endl;
        std::exit(1);
    }
    double num;
    for (int i = 0; i < m._rows; i++) {
        for (int j = 0; j < m._columns; j++) {
            if (m.get(i, j) > -pow(10, -15) && m.get(i, j) < pow(10, -15)) {
                num = 0;
            }
            else {
                num = m.get(i, j);
            }
            os << setprecision(5) << left << setw(8) << num << " ";
        }
        os << endl;
    }
    os << endl;
    return os;
}
void Matrix::random_fill(double min, double max) {
    if (_data == nullptr) {
        cerr << "Memory Error\nMatrix not specified" << endl;
        std::exit(1);
    }
    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < _columns; j++) {
            set(min + (double)rand() / RAND_MAX * (max - min), i, j);
        }
    }
}
void Matrix::set(double num, int row, int column) {
    if (_data == nullptr) {
        cerr << "Memory Error\nMatrix not specified" << endl;
        std::exit(1);
    }
    try {
        _data->at(row)->at(column) = num;
    }
    catch (const out_of_range& oor) {
        cerr << "Set error\nIndex is out of range" << endl;
        std::exit(1);
    }
}
double Matrix::get(int row, int column) const {
    if (_data == nullptr) {
        cerr << "Memory Error\nMatrix not specified" << endl;
        std::exit(1);
    }
    try {
        return _data->at(row)->at(column);
    }
    catch (const out_of_range& oor) {
        cerr << "Get error\nIndex is out of range" << endl;
        std::exit(1);
    }
}
Matrix Matrix::T() {
    if (_data == nullptr) {
        cerr << "Memory Error\nMatrix not specified" << endl;
        std::exit(1);
    }
    Matrix m_new(_columns, _rows);
    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < _columns; j++) {
            m_new.set(get(i, j), j, i);
        }
    }
    return m_new;
}
Matrix Matrix::operator+(Matrix B) {
    if (_data == nullptr || B._data == nullptr) {
        cerr << "Memory Error\nMatrix not specified" << endl;
        std::exit(1);
    }
    Matrix C(_rows, _columns);
    double num;
    if (_rows != B._rows || _columns != B._columns) {
        cerr << "Addition error. Matrix sizes do not match" << endl;
        std::exit(1);
    }
    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < _columns; j++) {
            num = get(i, j) + B.get(i, j);
            C.set(num, i, j);
        }
    }
    return C;
}
Matrix Matrix::operator+(double B) {
    if (_data == nullptr) {
        cerr << "Memory Error\nMatrix not specified" << endl;
        std::exit(1);
    }
    Matrix C(_rows, _columns);
    double num;
    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < _columns; j++) {
            num = get(i, j) + B;
            C.set(num, i, j);
        }
    }
    return C;
}
Matrix Matrix::operator -(Matrix B) {
    if (_data == nullptr || B._data == nullptr) {
        cerr << "Memory Error\nMatrix not specified" << endl;
        std::exit(1);
    }

    Matrix C(_rows, _columns);
    double num;
    if (_rows != B._rows || _columns != B._columns) {
        cerr << "Subtraction error\nMatrix sizes do not match" << endl;
        std::exit(1);
    }
    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < _columns; j++) {
            num = get(i, j) - B.get(i, j);
            C.set(num, i, j);
        }
    }
    return C;
}
Matrix Matrix::operator*(double B) {
    if (_data == nullptr) {
        cerr << "Memory Error\nMatrix not specified" << endl;
        std::exit(1);
    }

    Matrix C(_rows, _columns);
    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < _columns; j++) {
            C.set(get(i, j) * B, i, j);
        }
    }
    return C;
}
Matrix Matrix::operator/(double B) {
    if (_data == nullptr) {
        cerr << "Memory Error\nMatrix not specified" << endl;
        std::exit(1);
    }

    Matrix C(_rows, _columns);
    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < _columns; j++) {
            C.set(get(i, j) / B, i, j);
        }
    }
    return C;
}
Matrix Matrix::operator/(Matrix m) {
    if (_data == nullptr) {
        cerr << "Memory Error\nMatrix not specified" << endl;
        std::exit(1);
    }
    if (m.columns() != 1 || m._rows != 1) {
        cerr << "Error Can't divide matrix by matrix" << endl;
        std::exit(1);
    }

    Matrix C(_rows, _columns);
    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < _columns; j++) {
            C.set(get(i, j) / m.get(0, 0), i, j);
        }
    }
    return C;
}
Matrix Matrix::operator*(Matrix B) {
    if (_data == nullptr || B._data == nullptr) {
        cerr << "Memory Error\nMatrix not specified" << endl;
        std::exit(1);
    }
    Matrix C(_rows, B._columns);
    double sum;
    if (_columns != B._rows) {
        cerr << "Multiplication error\nMatrix sizes do not match" << endl;
        std::exit(1);
    }
    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < B._columns; j++) {
            sum = 0;
            for (int k = 0; k < _columns; k++) {
                sum += get(i, k) * B.get(k, j);
            }
            C.set(sum, i, j);
        }
    }
    return C;
}
Matrix Matrix::admult(Matrix B) {
    if (_data == nullptr || B._data == nullptr) {
        cerr << "Memory Error\nMatrix not specified" << endl;
        std::exit(1);
    }
    Matrix C(_rows, _columns);
    double num;
    if (_rows != B._rows || _columns != B._columns) {
        cerr << "Multiplication error\nMatrix sizes do not match" << endl;
        std::exit(1);
    }
    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < _columns; j++) {
            num = get(i, j) * B.get(i, j);
            C.set(num, i, j);
        }
    }
    return C;
}
double Matrix::trace() {
    if (_data == nullptr) {
        cerr << "Memory Error\nMatrix not specified" << endl;
        std::exit(1);
    }
    if (_rows != _columns) {
        cerr << "Trace error\nMatrix sizes do not match" << endl;
        std::exit(1);
    }
    double sum = 0;
    for (int i = 0; i < _rows; i++) {
        sum += get(i, i);
    }
    return sum;
}
double Matrix::det() {
    if (_data == nullptr) {
        cerr << "Memory Error\nMatrix not specified" << endl;
        std::exit(1);
    }
    Matrix temp;
    double* num = new double[_rows];
    double mult = 1;
    bool ok = 1;
    if (_rows != _columns) {
        cerr << "Det error\nMatrix sizes do not match" << endl;
        std::exit(1);
    }
    if (_rows == 1 && _columns == 1) {
        return get(0, 0);
    }

    temp.copy(*this);
    for (int i = 0; i < _rows; i++) { //elementary
        if (temp.get(i, i) == 0) {
            ok = 0;
            for (int j = 0; j < _rows; j++) {
                if (temp.get(j, i) != 0 && !ok) {
                    for (int k = 0; k < _rows; k++) {
                        temp.set(temp.get(j, k) + temp.get(i, k), i, k);
                        ok = 1;
                    }
                }
            }
        }

        if (!ok) {
            return 0;
        }
    }

    for (int i = 0; i < _rows - 1; i++) { //elimination
        for (int j = i + 1; j < _rows; j++) {
            for (int k = 0; k < _rows; k++) {
                if (temp.get(i, i) * temp.get(i, k) == 0) {
                    num[k] = temp.get(j, k);
                }
                else {
                    num[k] = temp.get(j, k) - temp.get(j, i) / temp.get(i, i) * temp.get(i, k);
                }
            }
            for (int k = 0; k < _rows; k++) {
                temp.set(num[k], j, k);
            }
        }
    }
    for (int k = 0; k < _rows; k++) {
        mult *= temp.get(k, k);
    }
    return mult;

}
double Matrix::dot(Matrix B) {
    if (_data == nullptr || B._data == nullptr) {
        cerr << "Memory Error\nMatrix not specified" << endl;
        std::exit(1);
    }

    Matrix temp_A;
    Matrix temp_B;
    double sum = 0;
    if (_columns == 1) {
        temp_A = T();
    }
    else if (_rows != 1) {
        cerr << "Dot error\nArgument is not a vector" << endl;
        std::exit(1);
    }
    else {
        temp_A = *this;
    }
    if (B._columns == 1) {
        temp_B = B.T();
    }
    else if (B._rows != 1) {
        cerr << "Dot error\nArgument is not a vector" << endl;
        std::exit(1);
    }
    else {
        temp_B = B;
    }
    if (temp_A._columns != temp_B._columns) {
        cerr << "Dot error\nVectors sizes do not match" << endl;
        std::exit(1);
    }
    for (int i = 0; i < temp_A._columns; i++) {
        sum += temp_A.get(0, i) * temp_B.get(0, i);
    }
    return sum;
}
double Matrix::norm() {
    if (_data == nullptr) {
        cerr << "Memory Error\nMatrix not specified" << endl;
        std::exit(1);
    }
    double sum = 0;
    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < _columns; j++) {
            sum += get(i, j) * get(i, j);
        }
    }
    return sqrt(sum);
}
double Matrix::max1() {
    if (_data == nullptr) {
        cerr << "Memory Error\nMatrix not specified" << endl;
        std::exit(1);
    }

    double max = get(0, 0);
    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < _columns; j++) {
            if (get(i, j) > max) {
                max = get(i, j);
            }
        }
    }
    return max;
}
double Matrix::sum() {
    double sum = 0;
    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < _columns; j++) {
            sum += get(i, j);
        }
    }
    return sum;
}
double Matrix::angle(Matrix B) {
    if (_data == nullptr || B._data == nullptr) {
        cerr << "Memory Error\nMatrix not specified" << endl;
        std::exit(1);
    }
    Matrix temp_A;
    Matrix temp_B;
    if (_columns == 1) {
        temp_A = T();
    }
    else if (_rows != 1) {
        cerr << "Dot error\nArgument is not a vector" << endl;
        std::exit(1);
    }
    else {
        temp_A = *this;
    }
    if (B._columns == 1) {
        temp_B = B.T();
    }
    else if (B._rows != 1) {
        cerr << "Dot error\nArgument is not a vector" << endl;
        std::exit(1);
    }
    else {
        temp_B = B;
    }
    return acos(temp_A.dot(temp_B) / (temp_A.norm() * temp_B.norm()));
}
int Matrix::rang() {
    if (_data == nullptr) {
        cerr << "Memory Error\nMatrix not specified" << endl;
        std::exit(1);
    }
    Matrix temp(_rows, _columns);
    double* num = new double[_columns];
    int rang = 0;
    bool skip;
    bool ok = 1;
    temp.copy(*this);
    for (int i = 0; i < _rows; i++) { //elementary
        if (temp.get(i, i) == 0) {
            ok = 0;
            for (int j = 0; j < _rows; j++) {
                if (temp.get(j, i) != 0 && !ok) {
                    for (int k = 0; k < _rows; k++) {
                        temp.set(temp.get(j, k) + temp.get(i, k), i, k);
                        ok = 1;
                    }
                }
            }
        }

        if (!ok) {
            return 0;
        }
    }
    for (int i = 0; i < _rows - 1; i++) {
        for (int j = i + 1; j < _rows; j++) {
            for (int k = 0; k < _columns; k++) {
                if (temp.get(i, i) * temp.get(i, k) == 0) {
                    num[k] = get(j, k);
                }
                else {
                    num[k] = get(j, k) - get(j, i) / get(i, i) * get(i, k);
                }
            }
            for (int k = 0; k < _columns; k++) {
                temp.set(num[k], j, k);
            }
        }
    }
    for (int i = 0; i < _rows; i++) {
        skip = 0;
        for (int j = 0; j < _columns; j++) {
            if (temp.get(i, j) != 0 && !skip) {
                rang++;
                skip = 1;
            }
        }
    }
    return rang;
}
Matrix Matrix::minor(int row, int column) {
    int i2 = 0;
    int j2 = 0;
    if (_data == nullptr) {
        cerr << "Memory Error\nMatrix not specified" << endl;
        std::exit(1);
    }
    if (_rows < 2 || _columns < 2) {
        cerr << "Minor error\nMatrix sizes do not match" << endl;
        std::exit(1);
    }
    Matrix m(_rows - 1, _columns - 1);
    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < _columns; j++) {
            if (i != row && j != column) {
                m.set(get(i, j), i2, j2);
                j2++;
                if (j2 >= _columns - 1) {
                    j2 = 0;
                    i2++;
                }
            }
        }
    }
    return m;
}
Matrix Matrix::inv() {
    if (_data == nullptr) {
        cerr << "Memory Error\nMatrix not specified" << endl;
        std::exit(1);
    }
    Matrix m(_rows, _columns);
    Matrix mnr;
    if (det() == 0) {
        cerr << "Inv error\nMatrix must be nonsingular" << endl;
        std::exit(1);
    }
    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < _columns; j++) {
            mnr = minor(i, j);
            m.set(pow(-1, (i + j)) * mnr.det(), i, j);
        }
    }
    m = m.T() / det();
    return m;
}
void Matrix::copy(Matrix m) {
    if (m._data == nullptr) {
        cerr << "Memory Error\nMatrix not specified" << endl;
        std::exit(1);
    }
    _rows = m._rows;
    _columns = m._columns;
    _data = new vector<vector<double>*>;
    for (int i = 0; i < _rows; i++) {
        _data->push_back(new vector<double>);
        for (int j = 0; j < _columns; j++) {
            _data->at(i)->push_back(m.get(i, j));
        }
    }
}
ifstream& operator<<(Matrix& m, ifstream& in) {
    string s;
    double num;
    bool first = 1;
    int j = 0;
    m._rows = 0;
    m._data = new vector<vector<double>*>;

    if (!in.is_open()) {
        cerr << "File not found" << endl;
        std::exit(1);
    }
    if (in >> num) {
        m._data->push_back(new vector<double>);
        m._data->at(m._rows)->push_back(num);
        j++;
    }
    else {
        cerr << "The file is damaged" << endl;
        std::exit(1);
    }

    while (in >> num) {
        m._data->at(m._rows)->push_back(num);
        j++;
        if (in.peek() == '\n') {
            m._data->push_back(new vector<double>);
            m._rows++;
            if (first) {
                m._columns = j;
                first = 0;
            }
            else if (m._columns != j) {
                cerr << "The file is damaged" << endl;
                std::exit(1);
            }
            j = 0;
        }
    }
    return in;
}
ofstream& operator>>(Matrix& m, ofstream& out) {
    if (m._data == nullptr) {
        cerr << "Memory Error\nMatrix not specified" << endl;
        std::exit(1);
    }
    if (!out.is_open()) {
        cerr << "File not found" << endl;
        std::exit(1);
    }
    for (int i = 0; i < m._rows; i++) {
        for (int j = 0; j < m._columns; j++) {
            out << m.get(i, j) << " ";
        }
        out << endl;
    }
    return out;
}
void Matrix::write_bin(ofstream& out) {
    if (!out.is_open()) {
        cerr << "File not found" << endl;
        std::exit(1);
    }
    if (_data == nullptr) {
        cerr << "Memory Error\nMatrix not specified" << endl;
        std::exit(1);
    }
    double num;
    out.write((char*)(&_rows), sizeof(int));
    out.write((char*)(&_columns), sizeof(int));
    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < _columns; j++) {
            num = get(i, j);
            out.write((char*)(&num), sizeof(double));
        }
    }
}
void Matrix::get_bin(ifstream& in) {
    if (!in.is_open()) {
        cerr << "File not found" << endl;
        std::exit(1);
    }
    double num;
    _data = new vector<vector<double>*>;
    in.read((char*)(&_rows), sizeof(int));
    in.read((char*)(&_columns), sizeof(int));
    for (int i = 0; i < _rows; i++) {
        _data->push_back(new vector<double>);
        for (int j = 0; j < _columns; j++) {
            if (!in.read((char*)(&num), sizeof(double))) {
                cerr << "The file is damaged" << endl;
                std::exit(1);
            }
            _data->at(i)->push_back(num);
        }
    }
}
int Matrix::rows() const {
    return _rows;
}
int Matrix::columns() const {
    return _columns;
}
void Matrix::add_matrix(Matrix m) {
    if (_data == nullptr) {
        copy(m);
        return;
    }
    if (_rows != m._rows) {
        cerr << "Add: Matrix sizes don't match" << endl;
        std::exit(1);
    }
    Matrix new_m(_rows, _columns + m._columns);
    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < _columns; j++) {
            new_m.set(get(i, j), i, j);
        }
        for (int j = 0; j < m._columns; j++) {
            new_m.set(m.get(i, j), i, j + _columns);
        }
    }
    copy(new_m);
}
Matrix Matrix::slice(int start_row, int start_column, int end_row, int end_column) {
    if (start_row<0 || start_row>_rows - 1 || end_row<0 || end_row>_rows - 1 ||
        start_column<0 || start_column>_columns - 1 || end_column<0 || end_column>_columns - 1 ||
        start_row > end_row || start_column > end_column) {
        cerr << "Slice: Index Error" << endl;
        std::exit(1);
    }
    Matrix new_m(end_row - start_row + 1, end_column - start_column + 1);
    for (int i = 0; i < end_row - start_row + 1; i++) {
        for (int j = 0; j < end_column - start_column + 1; j++) {
            new_m.set(get(start_row + i, start_column + j), i, j);
        }
    }
    return new_m;
}

Matrix::~Matrix() {};




class MATRIX_API DiagMatrix :public Matrix {
public:
    DiagMatrix(int N);
    void random_fill(double min, double max);
    void set(double num, int diag_c);
    DiagMatrix(const initializer_list<double>);
};


DiagMatrix::DiagMatrix(int N) :Matrix(N, N) {};
DiagMatrix::DiagMatrix(const initializer_list<double> l) :DiagMatrix(l.size()) {
    int i = 0;
    for (auto const& list : l) {
        set(list, i);
        i++;
    }
}
void DiagMatrix::random_fill(double min, double max) {
    for (int i = 0; i < _rows; i++) {
        set(min + (double)rand() / RAND_MAX * (max - min), i);
    }
}
void DiagMatrix::set(double num, int diag_c) {
    try {
        _data->at(diag_c)->at(diag_c) = num;
    }
    catch (const out_of_range& oor) {
        cerr << "Set error\nIndex is out of range" << endl;
        std::exit(1);
    }
}








class MATRIX_API DTriangleMatrix :public Matrix {
public:
    DTriangleMatrix(int N);
    void random_fill(double min, double max);
    void set(double num, int row, int column);
};


DTriangleMatrix::DTriangleMatrix(int N) :Matrix(N, N) {};
void DTriangleMatrix::random_fill(double min, double max) {
    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < _columns; j++) {
            if (i >= j) {
                set(min + (double)rand() / RAND_MAX * (max - min), i, j);
            }
        }
    }
}
void DTriangleMatrix::set(double num, int row, int column) {
    if (row < column) {
        cerr << "Set error\nIndex is out of range" << endl;
        std::exit(1);
    }
    try {
        _data->at(row)->at(column) = num;
    }
    catch (const out_of_range& oor) {
        cerr << "Set error\nIndex is out of range" << endl;
        std::exit(1);
    }
}





class MATRIX_API EyeMatrix :public Matrix {
public:
    EyeMatrix(int N);
};

EyeMatrix::EyeMatrix(int N) :Matrix(N, N) {
    for (int i = 0; i < N; i++) {
        set(1, i, i);
    }
}






class MATRIX_API SymMatrix :public Matrix {
public:
    SymMatrix(int N);
    void random_fill(double min, double max);
};

SymMatrix::SymMatrix(int N) : Matrix(N, N) {};
void SymMatrix::random_fill(double min, double max) {
    Matrix temp(_rows, _rows);
    temp.random_fill(0, 1);
    temp = (temp * temp.T() * (max - min)) / _rows + min;
    copy(temp);
}





class MATRIX_API UTriangleMatrix :public Matrix {
public:
    UTriangleMatrix(int N);
    void random_fill(double min, double max);
    void set(double num, int row, int column);
};

UTriangleMatrix::UTriangleMatrix(int N) :Matrix(N, N) {};
void UTriangleMatrix::random_fill(double min, double max) {
    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < _columns; j++) {
            if (i <= j) {
                set(min + (double)rand() / RAND_MAX * (max - min), i, j);
            }
        }
    }
}
void UTriangleMatrix::set(double num, int row, int column) {
    if (row > column) {
        cerr << "Set error\nIndex is out of range" << endl;
        std::exit(1);
    }
    try {
        _data->at(row)->at(column) = num;
    }
    catch (const out_of_range& oor) {
        cerr << "Set error\nIndex is out of range" << endl;
        std::exit(1);
    }
}






