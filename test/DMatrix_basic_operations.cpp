#include "../DMatrix.hpp"
#include <iostream>

int main() {
    Mi::DMatrix<double> A(2, 2);
    A(0, 0) = 1;
    A(0, 1) = 2;
    A(1, 0) = 3;
    A(1, 1) = 4;

    Mi::DMatrix<double> B(2, 2);
    B(0, 0) = 5;
    B(0, 1) = 6;
    B(1, 0) = 7;
    B(1, 1) = 8;

    Mi::DMatrix<double> C = A * B;

    std::cout << C << std::endl;

    return 0;
}