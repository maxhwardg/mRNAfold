
#ifndef MWMRNA_LIB_MD_ARRAY_HPP
#define MWMRNA_LIB_MD_ARRAY_HPP

#include <array>
#include <numeric>
#include <vector>

namespace mwmrna {
/// @brief Multidimensional array represented efficiently as a 1D array. Does not assume dimensions
/// are known at compile time.
/// @tparam T Data type
/// @tparam dims Number of dimensions of the array
template <typename T, int dims>
struct MdArray {
    std::vector<T> data;
    std::array<int, dims> shape;

    MdArray(T init, std::array<int, dims> shape) : shape(shape) {
        size_t n =
            std::accumulate(shape.begin(), shape.end(), size_t(1), std::multiplies<size_t>());
        data = std::vector<T>(n, init);
    }

    MdArray() = default;
    MdArray(const MdArray&) = default;
    MdArray(MdArray&&) = default;
    MdArray& operator=(const MdArray&) = default;

    T& operator[](const std::array<int, dims>&& indices) {
        size_t idx = 0;
        for (int i = 0; i < dims; ++i) {
            idx = idx * shape[i] + indices[i];
        }
        return data[idx];
    }
};
}  // namespace mwmrna

#endif  // MWMRNA_LIB_MD_ARRAY_HPP