#include "../Matrix.hpp"
#include <iostream>

int main() {
    Mi::DMatrix<double> A(3, 3);

    A(0, 0) = 3;
    A(0, 1) = 2;
    A(0, 2) = 2;
    A(1, 0) = 2;
    A(1, 1) = -3;
    A(1, 2) = -2;
    A(2, 0) =  4;
    A(2, 1) =  2;
    A(2, 2) =  3;

    Mi::DMatrix<double> B(3, 1);

    B(0, 0) = 1;
    B(1, 0) = -8;
    B(2, 0) = 0;

    Mi::LUDecomposition<Mi::DMatrix<double>> decomp(A);

    auto ans = decomp.solve(B);

    std::cout << ans(0, 0) << std::endl;
    std::cout << ans(1, 0) << std::endl;
    std::cout << ans(2, 0) << std::endl;

    // ans(0, 0) = -1;
    // ans(1, 0) =  2;
    // ans(2, 0) =  0;

    std::cout << "===========\n";
    auto re = A*ans;
    std::cout << re(0, 0) << std::endl;
    std::cout << re(1, 0) << std::endl;
    std::cout << re(2, 0) << std::endl;

    Mi::DMatrix<double> eye(3, 3);
    eye(0, 0) = 1;
    eye(1, 1) = 1;
    eye(2, 2) = 1;

    std::cout << "Inverse\n";
    auto inv = A.inverse();
    std::cout << inv << std::endl;

    std::cout << "===========\n";

    auto re2 = inv*B;
    std::cout << re2 << std::endl;

}