#pragma once
#include "MatrixBase.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace Mi {
template<class Element> 
class DMatrix {
private:
    Element* data;
public:
    const int rows;
    const int cols;
    using type = Element;
public:
    DMatrix(int rows, int cols = 1) : rows(rows), cols(cols) {
        data = new Element[rows * cols];
    }
    DMatrix(const DMatrix& rhs) : rows(rhs.rows), cols(rhs.cols) {
        data = new Element[rows * cols];
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                data[i * cols + j] = rhs(i, j);
            }
        }
    }
    ~DMatrix() {
        delete data;
    }

    Element& operator()(int row, int col = 0) {
        return data[row * cols + col];
    };
    Element operator()(int row, int col = 0) const {
        return data[row * cols + col];
    };

    DMatrix& operator=(const DMatrix& rhs) {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                data[i * cols + j] = rhs(i, j);
            }
        }
        return *this;
    };
    DMatrix operator+() const {
        return *this;
    };
    DMatrix operator-() const {
        DMatrix result(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result(i, j) = -data[i * cols + j];
            }
        }
        return result;
    };

    DMatrix inverse() const;
    
    Transpose<DMatrix> transpose() {
        return Transpose<DMatrix>(*this);
    }
};

template<typename T, typename Stream>
inline Stream& operator<<(Stream& stream, const DMatrix<T>& mat) {
    for(int i = 0; i < mat.rows; ++i) {
        for(int j = 0; j < mat.cols; ++j) {
            stream << mat(i, j) << " ";
        }
        stream << "\n";
    }
    return stream;
}

// Matrix-Matrix operations
// +=, -=, +, -, *, *=
template<class L, class R>
DMatrix<L>& operator+=(DMatrix<L>& lhs, const DMatrix<R>& rhs) {
    for (int i = 0; i < lhs.rows; ++i) {
        for (int j = 0; j < lhs.cols; ++j) {
            lhs(i, j) += rhs(i, j);
        }
    }
    return lhs;
};

template<class L, class R>
DMatrix<L>& operator-=(DMatrix<L>& lhs, const DMatrix<R>& rhs) {
    for (int i = 0; i < lhs.rows; ++i) {
        for (int j = 0; j < lhs.cols; ++j) {
            lhs(i, j) -= rhs(i, j);
        }
    }
    return lhs;
};

template<class L, class R>
DMatrix<decltype(L() + R())> operator+(const DMatrix<L>& lhs, const DMatrix<R>& rhs) {
    DMatrix<decltype(L() + R())> result = lhs;
    result += rhs;
    return result;
};

template<class L, class R>
DMatrix<decltype(L() + R())> operator-(const DMatrix<L>& lhs, const DMatrix<R>& rhs) {
    DMatrix<decltype(L() + R())> result = lhs;
    result -= rhs;
    return result;
};


template<class L, class R>
DMatrix<decltype(L() + R())> operator*(const DMatrix<L>& lhs, const DMatrix<R>& rhs) {
    DMatrix<decltype(L() + R())> result(lhs.rows, rhs.cols);
    for (int i = 0; i < lhs.rows; ++i) {
        for (int j = 0; j < rhs.cols; ++j) {
            result(i, j) = 0;
            for (int k = 0; k < lhs.cols; ++k) {
                result(i, j) += lhs(i, k) * rhs(k, j);
            }
        }
    }
    return result;
};

template<class L, class R>
DMatrix<L>& operator*=(const DMatrix<L>& lhs, const DMatrix<R>& rhs) {
    lhs = lhs * rhs;
    return lhs;
};

// Matrix-Scalar operations
// *=, /=, *, /
template<class L, class R>
DMatrix<L>& operator*=(DMatrix<L>& lhs, const R& rhs) {
    for (int i = 0; i < lhs.rows; ++i) {
        for (int j = 0; j < lhs.cols; ++j) {
            lhs(i, j) *= rhs;
        }
    }
    return lhs;
};

template<class L, class R>
DMatrix<L>& operator/=(DMatrix<L>& lhs, const R& rhs) {
    for (int i = 0; i < lhs.rows; ++i) {
        for (int j = 0; j < lhs.cols; ++j) {
            lhs(i, j) /= rhs;
        }
    }
    return lhs;
};

template<class L, class R>
DMatrix<L> operator*(const DMatrix<L>& lhs, const R& rhs) {
    DMatrix<L> result = lhs;
    result *= rhs;
    return result;
};

template<class L, class R>
DMatrix<L> operator/(const DMatrix<L>& lhs, const R& rhs) {
    DMatrix<L> result = lhs;
    result /= rhs;
    return result;
};

template<class Element>
class LUDecomposition<DMatrix<Element>> {
private:
public:
    int* idx;
    DMatrix<Element> mat;
    Element parity;
    bool singular;
    const int dim;

public:
    LUDecomposition(DMatrix<Element> mat_in): mat(mat_in), dim(mat.rows), parity(1), singular(false) {
        idx = new int[dim];

        for(int i = 0; i < dim; ++i) {
            idx[i] = i;
        }

        Element row_scale[dim];
        
        for(int i = 0; i < dim; ++i) {
            Element max_val = 0;
            for(int j = 0; j < dim; ++j) {
                max_val = std::max(max_val, std::abs(mat(i, j)));
            }
            if(max_val == 0) {
                singular = true;
                return;
            }
            row_scale[i] = 1 / max_val;
        }

        for(int j = 0; j < dim; ++j) {
            for(int i = 0; i < j; ++i) {
                Element sum = 0.0;

                for(int k = 0; k < i; ++k) {
                    sum += mat(i, k) * mat(k, j);
                }

                mat(i, j) -= sum;
            }

            for(int i = j; i < dim; ++i) {
                Element sum = 0.0;

                for(int k = 0; k < j; ++k) {
                    sum += mat(i, k) * mat(k, j);
                }

                mat(i, j) -= sum;
            }

            Element largest_elem = 0.0;
            int argmax = j;

            for(int i = j; i < dim; ++i) {
                Element this_elem = row_scale[i] * std::abs(mat(i, j));
                if(this_elem >= largest_elem) {
                    largest_elem = this_elem;
                    argmax = i;
                }
            }

            if(j != argmax) {
                for(int k = 0; k < dim; ++k) {
                    std::swap(mat(argmax, k), mat(j, k));
                }
                parity = -parity;
                std::swap(idx[j], idx[argmax]);
                row_scale[argmax] = row_scale[j];
            }

            if(mat(j, j) == 0.0) {
                singular = true;
                return;
            }

            if(j != dim) {
                Element pivot_inv = 1.0 / mat(j, j);
                for(int i = j + 1; i < dim; ++i) {
                    mat(i, j) *= pivot_inv;
                }
            }
        }
    }

    ~LUDecomposition() {
        delete idx;
    }

    DMatrix<Element> solve(DMatrix<Element>& b);
};

template<class Element>
DMatrix<Element> LUDecomposition<DMatrix<Element>>::solve(DMatrix<Element>& b) {
    DMatrix<Element> x(dim, b.cols);
    for(int c = 0; c < b.cols; ++c) {
        DMatrix<Element> tmp(dim);

        for(int i = 0; i < dim; ++i) {
            Element sum = 0.0;

            for(int j = 0; j < i; ++j) {
                sum += mat(i, j) * tmp(idx[j]);
            }

            tmp(idx[i]) = b(idx[i], c) - sum;
        }

        for(int i = dim-1; i>= 0; --i) {
            Element sum = 0.0;

            for(int j = i + 1; j < dim; ++j) {
                sum += mat(i, j) * tmp(idx[j]);
            }

            tmp(idx[i]) = (tmp(idx[i]) - sum) / mat(i, i);
        }

        for(int i = 0; i < dim; ++i) {
            x(i, c) = tmp(idx[i]);
        }
    }

    return x;
}
template<class Element>
DMatrix<Element> DMatrix<Element>::inverse() const {
    LUDecomposition<DMatrix<Element>> decomp(*this);
    DMatrix<Element> result(rows, rows);
    DMatrix<Element> b(rows);
    for(int i = 0; i < rows; ++i) {
        b(i) = 1;
        auto res = decomp.solve(b);
        for(int j = 0; j < rows; ++j) {
            result(j, i) = res(j, 0);
        }
        b(i) = 0;
    }
    return result;
}
} // namespace Mi