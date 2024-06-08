#pragma once

#include "DMatrix.hpp"
#include "SMatrix.hpp"

namespace Mi {
template <class Element, int Rows, int Cols>
DMatrix<Element> operator*(const DMatrix<Element> &lhs, const SMatrix<Element, Rows, Cols> &rhs) {
    DMatrix<Element> result(lhs.rows, rhs.cols);
    for (int i = 0; i < lhs.rows; ++i) {
        for (int j = 0; j < rhs.cols; ++j) {
            result(i, j) = 0;
            for (int k = 0; k < lhs.cols; ++k) {
                result(i, j) += lhs(i, k) * rhs(k, j);
            }
        }
    }
    return result;
}
template <class Element, int Rows, int Cols>
DMatrix<Element> operator*(const SMatrix<Element, Rows, Cols> &lhs, const DMatrix<Element> &rhs) {
    DMatrix<Element> result(lhs.rows, rhs.cols);
    for (int i = 0; i < lhs.rows; ++i) {
        for (int j = 0; j < rhs.cols; ++j) {
            result(i, j) = 0;
            for (int k = 0; k < lhs.cols; ++k) {
                result(i, j) += lhs(i, k) * rhs(k, j);
            }
        }
    }
    return result;
}
template <class Element, int Rows, int Cols>
DMatrix<Element> operator+=(const DMatrix<Element> &lhs, const SMatrix<Element, Rows, Cols> &rhs) {
    for (int i = 0; i < lhs.rows; ++i) {
        for (int j = 0; j < lhs.cols; ++j) {
            lhs(i, j) += rhs(i, j);
        }
    }
    return lhs;
}
template <class Element, int Rows, int Cols>
DMatrix<Element> operator+=(const SMatrix<Element, Rows, Cols> &lhs, const DMatrix<Element> &rhs) {
    for (int i = 0; i < lhs.rows; ++i) {
        for (int j = 0; j < lhs.cols; ++j) {
            lhs(i, j) += rhs(i, j);
        }
    }
    return lhs;
}

} // namespace Mi