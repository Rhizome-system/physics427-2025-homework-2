#include <tuple>
#include <cmath>
#include <limits>

using namespace std;

template <typename F>
tuple<bool, double> fixed_point_iteration(const F& f, double x0, double tolerance) {
    const int max_iter = 10000;
    double x = x0;

    for (int i = 0; i < max_iter; i++) {
        double x_next = f(x);
        if (isnan(x_next) || isinf(x_next)) return {false, x_next};
        if (fabs(f(x_next) - x_next) < tolerance) return {true, x_next};
        x = x_next;
    }
    return {false, x};
}


template <typename F>
tuple<bool, double> bisection(const F& f, double a, double b, double tolerance) {
    if (f(a) * f(b) >= 0.0) return {false, numeric_limits<double>::quiet_NaN()};

    const int max_iter = 1000;
    double fa = f(a);
    double mid = (a + b) / 2.0;

    for (int i = 0; i < max_iter; i++) {
        mid = (a + b) / 2.0;
        double fm = f(mid);

        if (isnan(fm) || isinf(fm)) return {false, fm};

        if (fabs(fm) < tolerance) return {true, mid};
        if (fa * fm < 0) {
            b = mid;
            fb = fm;
        } else {
            a = mid;
            fa = fm;
        }
    }
    return {false, mid};
}

template <typename F, typename Fprime>
tuple<bool, double> newton_raphson(const F& f, const Fprime& fprime, double x0, double tolerance) {
    const int max_iter = 1000;
    double x = x0;

    for (int i = 0; i < max_iter; i++) {
        double fx = f(x);
        double dfx = fprime(x);
        double x_next = x - fx / dfx;

        if (isnan(x_next) || isinf(x_next)) return {false, x_next};

        if (fabs(fx) < tolerance) return {true, x};

        x = x_next;
    }
    return {false, x};
}

