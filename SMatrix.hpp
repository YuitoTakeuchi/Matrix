#pragma once
#include "MatrixBase.hpp"

namespace Mi {
template<class Element, int Rows, int Cols>
class SMatrix {
public:
    Element data[Rows * Cols];
    const int rows = Rows;
    const int cols = Cols;

    SMatrix() {}
    ~SMatrix() {}

    Element& operator()(int row, int col) {
        return data[row * Cols + col];
    }
    Element operator()(int row, int col) const {
        return data[row * Cols + col];
    }

    SMatrix& operator=(const SMatrix& rhs) {
        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                data[i * Cols + j] = rhs(i, j);
            }
        }
        return *this;
    }
    SMatrix& operator+() const {
        return *this;
    }
    SMatrix& operator-() const {
        SMatrix result;
        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                result(i, j) = -data[i * Cols + j];
            }
        }
        return result;
    }
    template<class U>
    SMatrix& operator+=(const SMatrix<U, Rows, Cols>& rhs) {
        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                data[i * Cols + j] += rhs(i, j);
            }
        }
        return *this;
    }
    template<class U>
    SMatrix& operator-=(const SMatrix<U, Rows, Cols>& rhs) {
        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                data[i * Cols + j] -= rhs(i, j);
            }
        }
        return *this;
    }
    template<class U>
    SMatrix& operator*=(const U& scalar) {
        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                data[i * Cols + j] *= scalar;
            }
        }
        return *this;
    }
    template<class U>
    SMatrix& operator/=(const U& scalar) {
        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                data[i * Cols + j] /= scalar;
            }
        }
        return *this;
    }

    template<class U>
    SMatrix& operator*(const U& scalar) const {
        SMatrix result;
        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                result(i, j) = data[i * Cols + j] * scalar;
            }
        }
        return result;
    }

    template<class U>
    SMatrix& operator/(const U& scalar) const {
        SMatrix result;
        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                result(i, j) = data[i * Cols + j] / scalar;
            }
        }
        return result;
    }

};

template<class L, class R, int Rows, int Cols>
SMatrix<decltype(L() + R()), Rows, Cols> operator+(const SMatrix<L, Rows, Cols>& lhs, const SMatrix<R, Rows, Cols>& rhs) {
    SMatrix<decltype(L() + R()), Rows, Cols> result;
    for (int i = 0; i < Rows; ++i) {
        for (int j = 0; j < Cols; ++j) {
            result(i, j) = lhs(i, j) + rhs(i, j);
        }
    }
    return result;
}

template<class L, class R, int Rows, int Cols>
SMatrix<decltype(L() - R()), Rows, Cols> operator-(const SMatrix<L, Rows, Cols>& lhs, const SMatrix<R, Rows, Cols>& rhs) {
    SMatrix<decltype(L() - R()), Rows, Cols> result;
    for (int i = 0; i < Rows; ++i) {
        for (int j = 0; j < Cols; ++j) {
            result(i, j) = lhs(i, j) - rhs(i, j);
        }
    }
    return result;
}

template<class L, class R, int LRows, int LCols, int RRows, int RCols>
SMatrix<decltype(L() * R()), LRows, RCols> operator*(const SMatrix<L, LRows, LCols>& lhs, const SMatrix<R, RRows, RCols>& rhs) {
    static_assert(LCols == RRows, "Matrix dimensions do not match");
    SMatrix<decltype(L() * R()), LRows, RCols> result;
    for (int i = 0; i < LRows; ++i) {
        for (int j = 0; j < RCols; ++j) {
            result(i, j) = 0;
            for (int k = 0; k < LCols; ++k) {
                result(i, j) += lhs(i, k) * rhs(k, j);
            }
        }
    }
    return result;
}

template<class L, class R, int LRows, int LCols, int RRows, int RCols>
SMatrix<L, LRows, LCols>& operator*=(SMatrix<L, LRows, LCols>& lhs, const SMatrix<R, RRows, RCols>& rhs) {
    lhs = lhs * rhs;
    return lhs;
}

template<class Element, int Rows, int Cols>
class LUDecomposition<SMatrix<Element, Rows, Cols>> {

};

} // namespace Mi