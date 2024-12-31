#ifndef MWMRNA_LIB_ENERGY_SEMIRING_HPP
#define MWMRNA_LIB_ENERGY_SEMIRING_HPP

#include <cmath>
#include <limits>

namespace mwmrna {

/// @brief A semiring for energy values
/// @note Operations names are chosen to match a Boltzmann semiring
template <typename T, typename EnergyT>
concept EnergySemiring = requires(T a, EnergyT b) {
    { T(b) } -> std::same_as<T>;
    { T() } -> std::same_as<T>;
    { T::Val0() } -> std::same_as<T>;
    { T::Val1() } -> std::same_as<T>;
    { a + a } -> std::same_as<T>;
    { a + b } -> std::same_as<T>;
    { b + a } -> std::same_as<T>;
    { a += a } -> std::same_as<T&>;
    { a += b } -> std::same_as<T&>;
    { (a * a) } -> std::same_as<T>;
    { (a * b) } -> std::same_as<T>;
    { (b * a) } -> std::same_as<T>;
    { a < a } -> std::same_as<T>;
    { a < b } -> std::same_as<T>;
    { b < a } -> std::same_as<T>;
    { power(a, 1) } -> std::same_as<T>;
    { a / a } -> std::same_as<T>;
    { a / b } -> std::same_as<T>;
    { b / a } -> std::same_as<T>;
    { EnergyT(a) } -> std::same_as<EnergyT>;
    { a.Zero() } -> std::same_as<EnergyT>;
    { a.One() } -> std::same_as<EnergyT>;
    { a.Add(b, b) } -> std::same_as<EnergyT>;
    { a.Multiply(b, b) } -> std::same_as<EnergyT>;
    { a.LessThan(b, b) } -> std::same_as<EnergyT>;
    { a.Power(b, 1) } -> std::same_as<EnergyT>;
    { a.Divide(b, b) } -> std::same_as<EnergyT>;
    { a.MultMany(b, b) } -> std::same_as<EnergyT>;
};

template <typename EnergyT>
class MFESemiring {
   public:
    EnergyT value;
    MFESemiring(EnergyT value) : value(value) {}
    MFESemiring() : value(Val0().value) {}
    using energy_type = EnergyT;
    static MFESemiring<EnergyT> Val1() { return 0; }
    static MFESemiring<EnergyT> Val0() {
        // Divide by 3 to ensure that the sum of two Zero()s does not overflow
        return std::numeric_limits<EnergyT>::max() / 3;
    }
    EnergyT Zero() const { return Val0(); }
    EnergyT One() const { return Val1(); }
    EnergyT Add(EnergyT a, EnergyT b) const { return std::min(a, b); }
    EnergyT Multiply(EnergyT a, EnergyT b) const { return a + b; }
    EnergyT LessThan(EnergyT a, EnergyT b) const { return a > b; }
    EnergyT Power(EnergyT a, int n) const { return a * n; }
    EnergyT Divide(EnergyT a, EnergyT b) const { return a - b; }

    template <typename... Args>
    EnergyT MultMany(Args... a) const {
        return (EnergyT(a) + ...);
    }

    inline operator EnergyT() { return value; }
};

template <typename EnergyT>
inline MFESemiring<EnergyT> operator+(MFESemiring<EnergyT> a, MFESemiring<EnergyT> b) {
    return a.Add(a.value, b.value);
}

template <typename EnergyT>
inline MFESemiring<EnergyT> operator+(MFESemiring<EnergyT> a, EnergyT b) {
    return a.Add(a.value, b);
}

template <typename EnergyT>
inline MFESemiring<EnergyT> operator+(EnergyT a, MFESemiring<EnergyT> b) {
    return b.Add(a, b.value);
}

template <typename EnergyT>
inline MFESemiring<EnergyT>& operator+=(MFESemiring<EnergyT>& a, MFESemiring<EnergyT> b) {
    a.value = a.Add(a.value, b.value);
    return a;
}

template <typename EnergyT>
inline MFESemiring<EnergyT>& operator+=(MFESemiring<EnergyT>& a, EnergyT b) {
    a.value = a.Add(a.value, b);
    return a;
}

template <typename EnergyT>
inline MFESemiring<EnergyT> operator*(MFESemiring<EnergyT> a, MFESemiring<EnergyT> b) {
    return a.Multiply(a.value, b.value);
}

template <typename EnergyT>
inline MFESemiring<EnergyT> operator*(MFESemiring<EnergyT> a, EnergyT b) {
    return a.Multiply(a.value, b);
}

template <typename EnergyT>
inline MFESemiring<EnergyT> operator*(EnergyT a, MFESemiring<EnergyT> b) {
    return b.Multiply(a, b.value);
}

template <typename EnergyT>
inline MFESemiring<EnergyT> operator<(MFESemiring<EnergyT> a, MFESemiring<EnergyT> b) {
    return a.LessThan(a.value, b.value);
}

template <typename EnergyT>
inline MFESemiring<EnergyT> operator<(MFESemiring<EnergyT> a, EnergyT b) {
    return a.LessThan(a.value, b);
}

template <typename EnergyT>
inline MFESemiring<EnergyT> operator<(EnergyT a, MFESemiring<EnergyT> b) {
    return a.LessThan(a, b.value);
}

template <typename EnergyT>
inline MFESemiring<EnergyT> power(MFESemiring<EnergyT> a, int n) {
    return a.Power(a.value, n);
}

template <typename EnergyT>
inline MFESemiring<EnergyT> operator/(MFESemiring<EnergyT> a, MFESemiring<EnergyT> b) {
    return a.Divide(a.value, b.value);
}

template <typename EnergyT>
inline MFESemiring<EnergyT> operator/(MFESemiring<EnergyT> a, EnergyT b) {
    return a.Divide(a.value, b);
}

template <typename EnergyT>
inline MFESemiring<EnergyT> operator/(EnergyT a, MFESemiring<EnergyT> b) {
    return a.Divide(a, b.value);
}

template <typename EnergyT>
class BoltzSemiring {
   public:
    EnergyT value;
    BoltzSemiring(EnergyT value) : value(value) {}
    BoltzSemiring() : value(Val0().value) {}
    using energy_type = EnergyT;
    static BoltzSemiring<EnergyT> Val1() { return 1; }
    static BoltzSemiring<EnergyT> Val0() { return 0; }
    EnergyT Add(EnergyT a, EnergyT b) const { return a + b; }
    EnergyT Multiply(EnergyT a, EnergyT b) const { return a * b; }
    EnergyT Zero() const { return Val0(); }
    EnergyT One() const { return Val1(); }
    EnergyT LessThan(EnergyT a, EnergyT b) const { return a < b; }
    EnergyT Power(EnergyT a, int n) const { return std::pow(a, n); }
    EnergyT Divide(EnergyT a, EnergyT b) const { return a / b; }

    template <typename... Args>
    EnergyT MultMany(Args... a) const {
        return (EnergyT(a) * ...);
    }
    inline operator EnergyT() { return value; }
};

template <typename EnergyT>
inline BoltzSemiring<EnergyT> operator+(BoltzSemiring<EnergyT> a, BoltzSemiring<EnergyT> b) {
    return a.Add(a.value, b.value);
}

template <typename EnergyT>
inline BoltzSemiring<EnergyT> operator+(BoltzSemiring<EnergyT> a, EnergyT b) {
    return a.Add(a.value, b);
}

template <typename EnergyT>
inline BoltzSemiring<EnergyT> operator+(EnergyT a, BoltzSemiring<EnergyT> b) {
    return b.Add(a, b.value);
}

template <typename EnergyT>
inline BoltzSemiring<EnergyT>& operator+=(BoltzSemiring<EnergyT>& a, BoltzSemiring<EnergyT> b) {
    a.value = a.Add(a.value, b.value);
    return a;
}

template <typename EnergyT>
inline BoltzSemiring<EnergyT>& operator+=(BoltzSemiring<EnergyT>& a, EnergyT b) {
    a.value = a.Add(a.value, b);
    return a;
}

template <typename EnergyT>
inline BoltzSemiring<EnergyT> operator*(BoltzSemiring<EnergyT> a, BoltzSemiring<EnergyT> b) {
    return a.Multiply(a.value, b.value);
}

template <typename EnergyT>
inline BoltzSemiring<EnergyT> operator*(BoltzSemiring<EnergyT> a, EnergyT b) {
    return a.Multiply(a.value, b);
}

template <typename EnergyT>
inline BoltzSemiring<EnergyT> operator*(EnergyT a, BoltzSemiring<EnergyT> b) {
    return b.Multiply(a, b.value);
}

template <typename EnergyT>
inline BoltzSemiring<EnergyT> operator<(BoltzSemiring<EnergyT> a, BoltzSemiring<EnergyT> b) {
    return a.LessThan(a.value, b.value);
}

template <typename EnergyT>
inline BoltzSemiring<EnergyT> operator<(BoltzSemiring<EnergyT> a, EnergyT b) {
    return a.LessThan(a.value, b);
}

template <typename EnergyT>
inline BoltzSemiring<EnergyT> operator<(EnergyT a, BoltzSemiring<EnergyT> b) {
    return a.LessThan(a, b.value);
}

template <typename EnergyT>
inline BoltzSemiring<EnergyT> power(BoltzSemiring<EnergyT> a, int n) {
    return a.Power(a.value, n);
}

template <typename EnergyT>
inline BoltzSemiring<EnergyT> operator/(BoltzSemiring<EnergyT> a, BoltzSemiring<EnergyT> b) {
    return a.Divide(a.value, b.value);
}

template <typename EnergyT>
inline BoltzSemiring<EnergyT> operator/(BoltzSemiring<EnergyT> a, EnergyT b) {
    return a.Divide(a.value, b);
}

template <typename EnergyT>
inline BoltzSemiring<EnergyT> operator/(EnergyT a, BoltzSemiring<EnergyT> b) {
    return a.Divide(a, b.value);
}

}  // namespace mwmrna

#endif  // MWMRNA_LIB_ENERGY_SEMIRING_HPP