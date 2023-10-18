#include <vector>
#include <complex>
#include <Eigen/Core>
#include <Eigen/Sparse>

template<typename T>
class Interpolant {
public:
    Interpolant() {}
    Interpolant(const std::vector<double>& t, const std::vector<T>& y) {
        n = t.size();
        x = t;
        z.resize(n);
        b.resize(n-1);
        d.resize(n-1);

        // Compute slopes
        std::vector<double> h(n-1);
        std::vector<double> mu(n-1);
        std::vector<double> alpha(n-1);
        for (int i = 0; i < n-1; i++) {
            h[i] = x[i+1] - x[i];
            alpha[i] = 3/h[i]*(y[i+1]-y[i]) - 3/h[i-1]*(y[i]-y[i-1]);
        }

        // Compute coefficients
        std::vector<double> l(n), z(n);
        l[0] = 1;
        z[0] = 0;
        for (int i = 1; i < n-1; i++) {
            l[i] = 2*(x[i+1]-x[i-1]) - h[i-1]*mu[i-1];
            mu[i] = h[i]/l[i];
            z[i] = (alpha[i]-h[i-1]*z[i-1])/l[i];
        }
        l[n-1] = 1;
        z[n-1] = 0;
        for (int i = n-2; i >= 0; i--) {
            c[i] = z[i] - mu[i]*c[i+1];
            b[i] = (y[i+1]-y[i])/h[i] - h[i]*(c[i+1]+2*c[i])/3;
            d[i] = (c[i+1]-c[i])/(3*h[i]);
        }

        // Compute coefficients for linear interpolation
        for (int i = 0; i < n-1; i++) {
            a[i] = y[i];
            m[i] = (y[i+1]-y[i])/h[i];
        }
    }

    T evaluate(double t, int kind = 1) const {
        if (t < x[0]) {
            return T(0);
        } else if (t > x[n-1]) {
            return T(0);
        } else {
            int i = find_interval(t);
            double dx = t - x[i];
            if (kind == 0)
                return a[i] + m[i]*dx;
            return a[i] + b[i]*dx + c[i]*dx*dx + d[i]*dx*dx*dx;
        }
    }

    std::vector<T> evaluate(const std::vector<double>& t_target, int kind = 1) const {
        int N = t_target.size();
        std::vector<T> y_new(N);
        for (int i = 0; i < N; i++) {
            double t = t_target[i];
            if (t < x[0]) {
                y_new[i] = T(0);
            } else if (t > x[n-1]) {
                y_new[i] = T(0);
            } else {
                int j = find_interval(t);
                double dx = t - x[j];
                if (kind == 0)
                    y_new[i] = a[j] + m[j]*dx;
                else
                    y_new[i] = a[j] + b[j]*dx + c[j]*dx*dx + d[j]*dx*dx*dx;
            }
        }
        return y_new;
    }

private:
    int n;
    std::vector<double> x;
    std::vector<T> a;
    std::vector<T> b;
    std::vector<T> c;
    std::vector<T> z;
    std::vector<double> d;
    std::vector<double> m;

    int find_interval(double t) const {
        return std::upper_bound(x.begin(), x.end(), t) - x.begin() - 1;
    }
};