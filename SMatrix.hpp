#pragma once
#include "MatrixBase.hpp"

namespace Mi {
template<class Element, int Rows, int Cols = 1>
class SMatrix {
private:
    Element data[Rows * Cols] = {};
public:
    const int rows = Rows;
    const int cols = Cols;
    using type = Element;

    SMatrix() {}
    SMatrix(const SMatrix& rhs) {
        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                data[i * Cols + j] = rhs(i, j);
            }
        }
    }
    ~SMatrix() {}

    Element& operator()(int row, int col = 0) {
        return data[row * Cols + col];
    }
    Element operator()(int row, int col = 0) const {
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

    SMatrix inverse() const;

    Transpose<SMatrix<Element, Rows, Cols>> transpose() {
        return Transpose<SMatrix<Element, Rows, Cols>>(*this);
    }
};

template<typename T, int Rows, int Cols, typename Stream>
inline Stream& operator<<(Stream& stream, const SMatrix<T, Rows, Cols> mat) {
    for (int i = 0; i < Rows; ++i) {
        for (int j = 0; j < Cols; ++j) {
            stream << mat(i, j) << " ";
        }
        stream << "\n";
    }
    return stream;
}

template<class L, class R, int Rows, int Cols>
SMatrix<L, Rows, Cols> operator+=(SMatrix<L, Rows, Cols>& lhs, const SMatrix<R, Rows, Cols>& rhs) {
    for (int i = 0; i < Rows; ++i) {
        for (int j = 0; j < Cols; ++j) {
            lhs(i, j) += rhs(i, j);
        }
    }
    return lhs;
}

template<class L, class R, int Rows, int Cols>
SMatrix<L, Rows, Cols> operator-=(SMatrix<L, Rows, Cols>& lhs, const SMatrix<R, Rows, Cols>& rhs) {
    for (int i = 0; i < Rows; ++i) {
        for (int j = 0; j < Cols; ++j) {
            lhs(i, j) -= rhs(i, j);
        }
    }
    return lhs;
}

template<class L, class R, int Rows, int Cols>
SMatrix<decltype(L() + R()), Rows, Cols> operator+(const SMatrix<L, Rows, Cols>& lhs, const SMatrix<R, Rows, Cols>& rhs) {
    SMatrix<decltype(L() + R()), Rows, Cols> result = lhs;
    result += rhs;
    return result;
}

template<class L, class R, int Rows, int Cols>
SMatrix<decltype(L() - R()), Rows, Cols> operator-(const SMatrix<L, Rows, Cols>& lhs, const SMatrix<R, Rows, Cols>& rhs) {
    SMatrix<decltype(L() - R()), Rows, Cols> result = lhs;
    result -= rhs;
    return result;
}

template<class L, class R, int LRows, int LCols, int RRows, int RCols>
SMatrix<decltype(L() * R()), LRows, RCols> operator*(const SMatrix<L, LRows, LCols>& lhs, const SMatrix<R, RRows, RCols>& rhs) {
    static_assert(LCols == RRows, "Matrix size mismatch");
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
SMatrix<L, LRows, RCols> operator*=(SMatrix<L, LRows, LCols>& lhs, const SMatrix<R, RRows, RCols>& rhs) {
    lhs = lhs * rhs;
    return lhs;
}

template<class L, class R, int Rows, int Cols>
SMatrix<L, Rows, Cols> operator*=(SMatrix<L, Rows, Cols>& lhs, const R& rhs) {
    for (int i = 0; i < Rows; ++i) {
        for (int j = 0; j < Cols; ++j) {
            lhs(i, j) *= rhs;
        }
    }
    return lhs;
}

template<class L, class R, int Rows, int Cols>
SMatrix<L, Rows, Cols> operator/=(SMatrix<L, Rows, Cols>& lhs, const R& rhs) {
    for (int i = 0; i < Rows; ++i) {
        for (int j = 0; j < Cols; ++j) {
            lhs(i, j) /= rhs;
        }
    }
    return lhs;
}

template<class L, class R, int Rows, int Cols>
SMatrix<L, Rows, Cols> operator*(const SMatrix<L, Rows, Cols>& lhs, const R& rhs) {
    SMatrix<L, Rows, Cols> result = lhs;
    result *= rhs;
    return result;
}

template<class L, class R, int Rows, int Cols>
SMatrix<L, Rows, Cols> operator/(const SMatrix<L, Rows, Cols>& lhs, const R& rhs) {
    SMatrix<L, Rows, Cols> result = lhs;
    result /= rhs;
    return result;
}

template<class L, class R, int Rows, int Cols>
SMatrix<R, Rows, Cols> operator*(const L& lhs, const SMatrix<R, Rows, Cols>& rhs) {
    SMatrix<R, Rows, Cols> result = rhs;
    result *= lhs;
    return result;
}

template<class L, class R, int Rows, int Cols>
SMatrix<R, Rows, Cols> operator/(const L& lhs, const SMatrix<R, Rows, Cols>& rhs) {
    SMatrix<R, Rows, Cols> result = rhs;
    result /= lhs;
    return result;
}

template<class Element, int Rows, int Cols>
class LUDecomposition<SMatrix<Element, Rows, Cols>> {
private:
    static_assert(Rows == Cols, "Matrix must be square");
    int idx[Rows];
    SMatrix<Element, Rows, Cols> mat;
    bool singular;
    Element parity;

public:
    LUDecomposition(SMatrix<Element, Rows, Cols> mat_in): mat(mat_in), parity(1), singular(false) {
        for (int i = 0; i < Rows; ++i) {
            idx[i] = i;
        }

        Element row_scale[Rows];

        for (int i = 0; i < Rows; ++i) {
            Element max_val = 0;
            
            for(int j = 0; j < Rows; ++j) {
                max_val = std::max(max_val, std::abs(mat(i, j)));
            }

            if(max_val == 0) {
                singular = true;
                return;
            }
            row_scale[i] = 1 / max_val;
        }

        for(int j = 0; j < Rows; ++j) {
            for(int i = 0; i < j; ++i) {
                Element sum = 0.0;

                for(int k = 0; k < i; ++k) {
                    sum += mat(i, k) * mat(k, j);   
                }

                mat(i, j) -= sum;
            }

            for(int i = j; i < Rows; ++i) {
                Element sum = 0.0;

                for(int k = 0; k < j; ++k) {
                    sum += mat(i, k) * mat(k, j);
                }

                mat(i, j) -= sum;
            }

            Element largest_elem = 0.0;
            int argmax = j;

            for(int i = j; i < Rows; ++i) {
                Element this_elem = row_scale[i] * std::abs(mat(i, j));
                if(this_elem >= largest_elem) {
                    largest_elem = this_elem;
                    argmax = i;
                }
            }

            if(j != argmax) {
                for(int k = 0; k < Rows; ++k) {
                    std::swap(mat(argmax, k), mat(j, k));
                }
                parity = -parity;
                std::swap(idx[argmax], idx[j]);
                row_scale[argmax] = row_scale[j];
            }

            if(mat(j, j) == 0.0) {
                singular = true;
                return;
            }

            if(j != Rows) {
                Element pivot_inv = 1.0 / mat(j, j);
                for(int i = j+1; i < Rows; ++i) {
                    mat(i, j) *= pivot_inv;
                }
            }
        }
    }

    SMatrix<Element, Rows> solve(SMatrix<Element, Rows>& b) const;
};

template<class Element, int Rows, int Cols>
SMatrix<Element, Rows> LUDecomposition<SMatrix<Element, Rows, Cols>>::solve(SMatrix<Element, Rows>& b) const {
    SMatrix<Element, Rows> x;
    SMatrix<Element, Rows> tmp;

    for(int i = 0; i < Rows; ++i) {
        Element sum = 0.0;

        for(int j = 0; j < i; ++j) {
            sum += mat(i, j) * tmp(idx[j]);
        }

        tmp(idx[i]) = b(idx[i]) - sum;
    }

    for(int i = Rows-1; i >= 0; --i) {
        Element sum = 0.0;

        for(int j = i + 1; j < Rows; ++j) {
            sum += mat(i, j) * tmp(idx[j]);
        }

        tmp(idx[i]) = (tmp(idx[i]) - sum) / mat(i, i);
    }

    for(int i = 0; i < Rows; ++i) {
        x(i) = tmp(idx[i]);
    }

    return x;
}

template<class Element, int Rows, int Cols>
SMatrix<Element, Rows, Cols> SMatrix<Element, Rows, Cols>::inverse() const {
    LUDecomposition<SMatrix<Element, Rows, Cols>> decomp(*this);
    SMatrix<Element, Rows, Cols> result;
    SMatrix<Element, Rows> b;

    for(int i = 0; i < Rows; ++i) {
        b(i) = 1;
        SMatrix<Element, Rows> x = decomp.solve(b);
        for(int j = 0; j < Rows; ++j) {
            result(j, i) = x(j);
        }
        b(i) = 0;
    }
    return result;
}
} // namespace Mi