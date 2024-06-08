#pragma once

namespace Mi {
template<class Mat>
class Transpose {
public:
    Mat& mat;
    Transpose(Mat& mat) : mat(mat) {}

    typename Mat::type operator()(int i, int j) const {
        return mat(j, i);
    }

    typename Mat::type& operator()(int i, int j) {
        return mat(j, i);
    }
};

template<class Mat>
class LUDecomposition {};
} // namespace Mi