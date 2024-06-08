#nullable enable
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.IO;
using System.Text;
// ReSharper disable ArrangeRedundantParentheses

#pragma warning disable 219
#pragma warning disable CS0162
namespace AmericanOptionAlbena
{
    [SuppressMessage("ReSharper", "InconsistentNaming")]
    [SuppressMessage("ReSharper", "CommentTypo")]
    [SuppressMessage("ReSharper", "ParameterHidesMember")]
    [SuppressMessage("ReSharper", "SuggestBaseTypeForParameter")]
    public static class Program
    {
        private const int l_max_iterations = 100; // max iteration on l
        private const double Tol = 10e-10; // eps to refine eta
        private const double eps = 1e-3; // eps to refine rho
        private const double lb = 0d; // left bound
        private const double rb = 2d; // right bound
        private const double T0 = 0d; // the start time
        private const double Tn = 3d; // the finish time
        private const double T = Tn - T0; // time interval
        private const double sigma = 0.2d; // the sigma = volatility
        private const double sigma2 = sigma * sigma; // the squared sigma
        private const double r = 0.08d; // the risk-free rate
        private const double K = 100d; // the strike price
        private const double q = 0d; // the dividend rate
        private const int N = 400; // the number of space intervals
        private const int N1 = N + 1; // the number of points
        private const double b = -2d * r / sigma2;

        private const double h = (rb - lb) / N; // the space step

        private const int TauTuner = 1000; // the value to change the static time step
        private const int M = (int) (T * TauTuner); // the number of time intervals
        private const double tau0 = T / M; // the time step
        private const double alpha = 2d; // condense parameter for t
        private const double beta = 2.5d; // condense parameter for h

        private const bool print = true; // the parameter which allows to enable/disable printing of the results
#pragma warning disable CS0414 // Field is assigned but its value is never used
        private static readonly bool finite_element_enabled = false; // the parameter which allows to enable/disable the finite element at x_N
#pragma warning restore CS0414 // Field is assigned but its value is never used
        private static readonly bool enable_new_eta = true; // the parameter which allows to enable/disable a new formula for eta

        private static readonly bool s0_economic_style = true; // the parameter which allows to enable/disable graphs in economic style (S is Y-axis)

        private static readonly double[] taus = GetTaus();
        private static readonly double[] hs = GetHs();

        private static string? s;

        private static void print_to_log(string val)
        {
            s += val + Environment.NewLine;
        }

        private static void clear_log()
        {
            s = string.Empty;
        }

        private static void delete_log()
        {
            File.Delete("log.txt");
        }

        private static void flush_log()
        {
            Console.Out.Flush();
            using var writer = new StreamWriter("log.txt", true);
            writer.WriteLine(s);
            Console.WriteLine(s);
            Console.Out.Flush();
        }

        public static void Main()
        {
            delete_log();
            // _ = Run(false, false, new[]
            // {
            //     90d, 100d, 110d, 120d
            // }, finite_element_enabled);
            // flush_log();
            // clear_log();
            // _ = Run(true, false, new[]
            // {
            //     90d, 100d, 110d, 120d
            // }, finite_element_enabled);
            // flush_log();
            // clear_log();
            _ = Run(false, true, new[]
            {
                90d, 100d, 110d, 120d
            }, true);
            flush_log();
            clear_log();
            // _ = Run(false, true, new[]
            // {
            //     90d, 100d, 110d, 120d
            // }, false);
            // flush_log();
            // clear_log();
            // _ = Run(true, true, new[]
            // {
            //     90d, 100d, 110d, 120d
            // }, finite_element_enabled);
            // flush_log();
            // clear_log();
            // var name4 = Run(true, true, new[]
            // {
            //     80d, 90d, 100d, 110d, 120d
            // }, finite_element_enabled);
        }

        private static string? Run(bool refined_tau, bool refined_h, double[] S_vals, bool finite_element)
        {
            if (print)
                print_to_log(
                    $"Tn = {Tn} M = {M} N = {N} rb = {rb} sigma = {sigma} r = {r} K = {K} q = {q} h = {h} tau = {tau0} refined h = {refined_h} refined tau = {refined_tau} alpha = {alpha} beta = {beta} tau intervals = {(refined_tau ? taus.Length.ToString() : "-")} finite element enabled = {finite_element}");

            Debug.Assert(q < r);
            var u_curr = new double[N1];
            var u_prev = new double[N1];
            var lambda = new double[M + 1];
            var a = new double[M + 1];
            var s0_dash = new double[M + 1];
            s0_dash[0] = 0d;
            a[0] = (2d * r - 2d * q - sigma2 + 2d * 0d) / sigma2;
            lambda[0] = Math.Sqrt(a[0] * a[0] - 4d * b);
            var tau_sum = 0d;
            var printed1 = false;
            var printed2 = false;
            var nan_break = false;
            for (var j = 1; j <= M; j++)
            {
                if (nan_break)
                    break;
                var tau = refined_tau ? taus[j - 1] : tau0;
                tau_sum += tau;
                if (tau_sum >= 1d && print && refined_tau && !printed1)
                {
                    print_to_log($"Sum of tau on 0..1 = {tau_sum}, j = {j}");
                    printed1 = true;
                }

                if (tau_sum >= T && print && refined_tau && !printed2)
                {
                    print_to_log($"Sum of tau on 0..T = {tau_sum}, j = {j}");
                    printed2 = true;
                }

                var eta_j = new double[l_max_iterations];
                var rho_j = new double[l_max_iterations];
                rho_j[0] = s0_dash[j - 1];
                var l = 0;
                while (l < l_max_iterations - 1)
                {
                    if (CheckNaN(rho_j[l], l, j))
                    {
                        nan_break = true;
                        break;
                    }

                    var u_0j = K * (1d - Math.Exp(rho_j[l]));
                    var result = Solve(l, j, u_prev, rho_j, s0_dash, a, lambda, tau, u_0j, refined_tau, refined_h, finite_element);
                    u_curr = result.u_curr;
                    var h0 = refined_h ? hs[1] : h;
                    if (enable_new_eta)
                        eta_j[l] = K * Math.Exp(rho_j[l]) * (1d + (q * h0 / sigma2) + h0 / 2d) + (u_curr[0] - u_0j) / h0 - (h0 * r * K) / sigma2;
                    else
                        eta_j[l] = K * Math.Exp(rho_j[l]) + (u_curr[0] - u_0j) / h0;

                    // if (print)
                    //     print_to_string($"l = {l:G10} j = {j:G10} tau {tau} rho_j_l {rho_j[l]:N10} alpha_j {result.alpha_j:N5} u[1] {u_curr[1]:G10} u[0] {u_curr[0]:G10} eta_j_l = {eta_j[l]:G10} u[N] = {u_curr[N]:G10}");

                    if (Math.Abs(eta_j[l]) <= Tol)
                    {
                        s0_dash[j] = rho_j[l];
                        break;
                    }

                    l++;
                    rho_j[l] = l == 1 ? rho_j[l - 1] - eps : (eta_j[l - 1] * rho_j[l - 2] - eta_j[l - 2] * rho_j[l - 1]) / (eta_j[l - 1] - eta_j[l - 2]);
                }

                for (var i = 0; i < u_prev.Length; i++)
                    u_prev[i] = u_curr[i];
            }

            if (nan_break)
                return null;

            var s0 = GetS0(s0_dash);
            var zone = $"{nameof(s0_dash)}_h_condensed_{refined_h}_tau_condensed_{refined_tau}_finite_elem_{finite_element}_rb_{rb}_new_eta_{enable_new_eta}_beta_{beta}";
            var name =
                $"{GetPrefix(refined_tau, refined_h)}_{nameof(s0_dash)}_K={K}_N1={N1}_T={T}_h_condensed_{refined_h}_tau_condensed_{refined_tau}_finite_elem_{finite_element}_rb_{rb}_new_eta_{enable_new_eta}_beta_{beta}_is_economic_{s0_economic_style}.";
            // S0DashDirectTime(name + ".dat", zone, s0_dash, tau0, taus, refined_tau);

            zone = $"{nameof(s0)}_h_condensed_{refined_h}_tau_condensed_{refined_tau}_finite_elem_{finite_element}_rb_{rb}_new_eta_{enable_new_eta}_beta_{beta}";
            name =
                $"{GetPrefix(refined_tau, refined_h)}_s0_T-t_K={K}_N1={N1}_T={T}_h_condensed_{refined_h}_tau_condensed_{refined_tau}_finite_elem_{finite_element}_rb_{rb}_new_eta_{enable_new_eta}_beta_{beta}_is_economic_{s0_economic_style}.";
            if (s0_economic_style)
            {
                S0ReversedTimeE(name + "dat", zone, s0, T, tau0, taus, refined_tau);
                test("test_" + name + "dat", s0, T, tau0, taus, refined_tau);
            }
            else
                S0ReversedTime(name + "dat", zone, s0, T, tau0, taus, refined_tau);

            // CreateLay("graph.lay", name + "lay", name + "dat");

            // S (90,100,110,120)
            // Для заданных значений S вычислим значения x=ln(S/s0(tau)). При t=T у нас tau=0, т.е. s0(0)=84..., подставляем в эту формулу S=90, 100,...,
            // находим соответствующие им значения x и в этих точках выводим значения решения u.
            var s00 = get_s0_0(s0, refined_tau ? taus : null, T);
            print_to_log("s0(0) = " + s00.ToString(CultureInfo.InvariantCulture));
            foreach (var val in S_vals)
                print_results(val, u_curr, s00, refined_h);

            CheckS0Dash(s0_dash);
            CheckS0(s0);
            print_to_log("-----------------");
            return name;
        }

        private static double get_s0_0(double[] s0, double[]? taus, double T)
        {
            if (taus == null)
                return s0[s0.Length - 1];
            var t = T;
            var j = 0;
            while (t >= 0)
                t -= taus[j++];
            return j == s0.Length ? s0[j - 1] : s0[j];
        }

        private static void print_results(double S, double[] u, double s0, bool refined_h)
        {
            var x = Math.Log(S / s0);
            get_x_i(x, refined_h, out var x_i, out var x_im1, out var i);
            get_x_rb(rb, refined_h, out var x_rb, out var i_rb);
            // print_to_log($"calculated x = {x} found x_i = {x_i} {(refined_h ? string.Empty : $"h = {h}")}");
            string format;
            if (i < u.Length && i - 1 >= 0)
            {
                var u_x = ((u[i] - u[i - 1]) / (x_i - x_im1)) * x + (u[i - 1] * x_i - u[i] * x_im1) / (x_i - x_im1);
                format = $"S = {S,3}, rb = {rb} x_rb = {x_rb} i_rb = {i_rb} u(x_rb) = {u[i_rb]}, x_i-1 = {x_im1,4:F17} x = {x,4:F17} x_i = {x_i,4:F17} u(x_(i-1)) = {(i - 1 >= 0 ? u[i - 1] : -1d),4:F17} u(x) = {u_x,4:F17} u(x_i) = {u[i],4:F17}";
            }
            else
            {
                format = $"The index 'i' is out of 'u' boundaries: i = {i} len(u) = {u.Length}";
            }

            print_to_log(format);
            return;

            static void get_x_i(double val, bool refined, out double x_j, out double x_jm1, out int j)
            {
                x_j = 0d;
                j = 0;
                while (true)
                {
                    x_jm1 = x_j;
                    if (refined)
                    {
                        if (j < hs.Length)
                            x_j += hs[j + 1];
                        else
                            x_j += h;
                    }
                    else
                        x_j += h;

                    if (x_j >= val)
                        break;
                    j++;
                }
            }
            
            static void get_x_rb(double boundary, bool refined, out double x_j, out int j)
            {
                x_j = 0d;
                j = 0;
                while (true)
                {
                    if (refined)
                    {
                        if (j < hs.Length)
                            x_j += hs[j + 1];
                        else
                            x_j += h;
                    }
                    else
                        x_j += h;

                    if (x_j >= boundary)
                        break;
                    j++;
                }
            }
        }

        private static (double[] u_curr, double alpha_j) Solve(int l, int j, double[] u_prev, double[] rho_j, double[] s0_dash, double[] a, double[] lambda, double in_tau, double u_0j,
            bool refined_tau, bool refined_h, bool finite_element)
        {
            // alpha, a, lambda, nu
            var s_approx = (rho_j[l] - s0_dash[j - 1]) / in_tau;
            a[j] = (2d * r - 2d * q - sigma2 + 2d * s_approx) / sigma2;
            lambda[j] = Math.Sqrt(a[j] * a[j] - 4d * b);

            var alpha_j = r - q - sigma2 / 2d + s_approx;
            var f = new double[N1];
            var h_12 = refined_h ? hs[1] : h;
            var h_32 = refined_h ? hs[2] : h;
            if (alpha_j >= 0d)
                f[0] = u_prev[0] / in_tau + (u_0j * sigma2) / (h_12 * (h_12 + h_32));
            else
                f[0] = u_prev[0] / in_tau + (sigma2 / (h_12 * (h_12 + h_32)) - alpha_j / h_12) * u_0j;
            for (var i = 1; i < N1 - 1; i++)
                f[i] = u_prev[i] / in_tau;

            if (finite_element)
            {
                var nu_j = 0d;
                if (j > 1)
                {
                    var tau_j_32 = refined_tau ? taus[j - 2] : tau0;
                    var tau_j_12 = refined_tau ? taus[j - 1] : tau0;
                    var val = 2d * ((s0_dash[j - 1] - s0_dash[j - 2]) / tau_j_32 - (rho_j[l] - s0_dash[j - 1]) / tau_j_12);
                    nu_j = 2d / (lambda[j - 1] + lambda[j] + val / sigma2);
                }

                f[N1 - 1] = (0.5d + nu_j / (refined_h ? hs[N] : h)) * (u_prev[N1 - 1] / in_tau);
            }
            else
                f[N1 - 1] = 0d;

            var a0 = new double[N1]; // a - below main diagonal (indexed as [1;n-1])
            var b0 = new double[N1]; // b - main diagonal of matrix (indexed as [0;n-1])
            var c0 = new double[N1]; // c - up to main diagonal (indexed as [0;n-2])
            a0[0] = 0d;
            for (var i = 1; i < N1 - 1; ++i)
            {
                if (refined_h)
                {
                    if (alpha_j >= 0d)
                        a0[i] = -sigma2 / (hs[i + 1] * (hs[i + 1] + hs[i + 2]));
                    else
                        a0[i] = -sigma2 / (hs[i + 1] * (hs[i + 1] + hs[i + 2])) + alpha_j / hs[i + 1];
                }
                else
                {
                    if (alpha_j >= 0d)
                        a0[i] = -sigma2 / (2d * h * h);
                    else
                        a0[i] = -sigma2 / (2d * h * h) + alpha_j / h;
                }
            }

            if (finite_element)
            {
                if (refined_h)
                {
                    if (alpha_j >= 0d)
                        a0[N1 - 1] = -sigma2 / (2d * hs[N] * hs[N]);
                    else
                        a0[N1 - 1] = -sigma2 / (2d * hs[N] * hs[N]) + alpha_j / (2d * hs[N]);
                }
                else
                {
                    if (alpha_j >= 0d)
                        a0[N1 - 1] = -sigma2 / (2d * h * h);
                    else
                        a0[N1 - 1] = -sigma2 / (2d * h * h) + alpha_j / (2d * h);
                }
            }
            else
                a0[N1 - 1] = 0d;

            if (refined_h)
            {
                if (alpha_j >= 0d)
                    b0[0] = sigma2 / (hs[1] + hs[2]) * (1d / hs[1] + 1d / hs[2]) + (r + 1d / in_tau) + alpha_j / hs[2];
                else
                    b0[0] = sigma2 / (hs[1] + hs[2]) * (1d / hs[1] + 1d / hs[2]) + (r + 1d / in_tau) - alpha_j / hs[1];
            }
            else
            {
                if (alpha_j >= 0d)
                    b0[0] = sigma2 / (h * h) + (r + 1d / in_tau) + alpha_j / h;
                else
                    b0[0] = sigma2 / (h * h) + (r + 1d / in_tau) - alpha_j / h;
            }

            for (var i = 1; i < N1 - 1; ++i)
            {
                if (refined_h)
                {
                    if (alpha_j >= 0d)
                        b0[i] = sigma2 / (hs[i + 1] + hs[i + 2]) * (1d / hs[i + 1] + 1d / hs[i + 2]) + (r + 1d / in_tau) + alpha_j / hs[i + 2];
                    else
                        b0[i] = sigma2 / (hs[i + 1] + hs[i + 2]) * (1d / hs[i + 1] + 1d / hs[i + 2]) + (r + 1d / in_tau) - alpha_j / hs[i + 1];
                }
                else
                {
                    if (alpha_j >= 0d)
                        b0[i] = sigma2 / (h * h) + (r + 1d / in_tau) + alpha_j / h;
                    else
                        b0[i] = sigma2 / (h * h) + (r + 1d / in_tau) - alpha_j / h;
                }
            }

            if (finite_element)
            {
                if (refined_h)
                {
                    if (alpha_j >= 0)
                    {
                        var val1 = sigma2 / 2d * (1d / (hs[N] * hs[N]) + (a[j] + lambda[j]) / (2d * hs[N]));
                        var val2 = 0.5d * (1d / in_tau + r) + 1d / (lambda[j] * in_tau * hs[N]);
                        var val3 = alpha_j / (2d * hs[N]) * (Math.Exp(-(lambda[j] + a[j]) * hs[N] / 2d) - 1d);
                        b0[N1 - 1] = val1 + val2 - val3;
                    }
                    else
                    {
                        var val1 = sigma2 / 2d * (1d / (hs[N] * hs[N]) + (a[j] + lambda[j]) / (2d * hs[N]));
                        var val2 = 0.5d * (1d / in_tau + r) + 1d / (lambda[j] * in_tau * hs[N]);
                        var val3 = alpha_j / (2d * hs[N]);
                        b0[N1 - 1] = val1 + val2 - val3;
                    }
                }
                else
                {
                    if (alpha_j >= 0)
                    {
                        var val1 = sigma2 / 2d * (1d / (h * h) + (a[j] + lambda[j]) / (2d * h));
                        var val2 = 0.5d * (1d / in_tau + r) + 1d / (lambda[j] * in_tau * h);
                        var val3 = alpha_j / (2d * h) * (Math.Exp(-(lambda[j] + a[j]) * h / 2d) - 1d);
                        b0[N1 - 1] = val1 + val2 - val3;
                    }
                    else
                    {
                        var val1 = sigma2 / 2d * (1d / (h * h) + (a[j] + lambda[j]) / (2d * h));
                        var val2 = 0.5d * (1d / in_tau + r) + 1d / (lambda[j] * in_tau * h);
                        var val3 = alpha_j / (2d * h);
                        b0[N1 - 1] = val1 + val2 - val3;
                    }
                }
            }
            else
                b0[N1 - 1] = 1d;

            if (refined_h)
            {
                if (alpha_j >= 0d)
                    c0[0] = -sigma2 / (hs[2] * (hs[1] + hs[2])) - alpha_j / hs[2];
                else
                    c0[0] = -sigma2 / (hs[2] * (hs[1] + hs[2]));
            }
            else
            {
                if (alpha_j >= 0d)
                    c0[0] = -sigma2 / (2d * h * h) - alpha_j / h;
                else
                    c0[0] = -sigma2 / (2d * h * h);
            }

            for (var i = 1; i < N1 - 2; ++i)
            {
                // equation N - 1 = the [N1 - 3] index
                if (i == N1 - 3)
                {
                    if (finite_element)
                    {
                        if (refined_h)
                        {
                            if (alpha_j >= 0d)
                                c0[N1 - 3] = -sigma2 / (hs[N1 - 3 + 2] * (hs[N1 - 3 + 1] + hs[N1 - 3 + 2])) - alpha_j / hs[N1 - 3 + 2];
                            else
                                c0[N1 - 3] = -sigma2 / (hs[N1 - 3 + 2] * (hs[N1 - 3 + 1] + hs[N1 - 3 + 2]));
                        }
                        else
                        {
                            if (alpha_j >= 0d)
                                c0[N1 - 3] = -sigma2 / (2d * h * h) - alpha_j / h;
                            else
                                c0[N1 - 3] = -sigma2 / (2d * h * h);
                        }
                    }
                    else
                        c0[N1 - 3] = 0d;
                }
                else
                {
                    if (refined_h)
                    {
                        if (alpha_j >= 0d)
                            c0[i] = -sigma2 / (hs[i + 2] * (hs[i + 1] + hs[i + 2])) - alpha_j / hs[i + 2];
                        else
                            c0[i] = -sigma2 / (hs[i + 2] * (hs[i + 1] + hs[i + 2]));
                    }
                    else
                    {
                        if (alpha_j >= 0d)
                            c0[i] = -sigma2 / (2d * h * h) - alpha_j / h;
                        else
                            c0[i] = -sigma2 / (2d * h * h);
                    }
                }
            }

            c0[N1 - 2] = 0d;

            var u = SolveByTridiagonalMatrixAlgorithm(j, N1, a0, b0, c0, f);
            return (u, alpha_j);
        }

        private static string GetPrefix(bool tau, bool h)
        {
            if (!h && !tau)
                return new string('1', 1);
            if (h && !tau)
                return new string('2', 2);
            return !h ? new string('3', 3) : new string('4', 4);
        }

        private static double gh1(double x, double in_beta)
        {
            return x <= 0d ? 0d : Math.Pow(x, in_beta);
        }

        private static double gt1(double t, double in_alpha)
        {
            return t <= 0d ? 0d : Math.Pow(t, in_alpha);
        }

        private static double gh2(double x, double in_beta)
        {
            switch (x)
            {
                case <= 0d:
                    return 0d;
                case >= 1d:
                {
                    var d = in_beta;
                    var e = 1d - in_beta;
                    return d * x + e;
                }
                default:
                    return Math.Pow(x, in_beta);
            }
        }

        private static double gt2(double t, double in_alpha)
        {
            switch (t)
            {
                case <= 0d:
                    return 0d;
                case >= 1d:
                {
                    var d = in_alpha;
                    var e = 1d - in_alpha;
                    return d * t + e;
                }
                default:
                    return Math.Pow(t, in_alpha);
            }
        }

        private static double[] SolveByTridiagonalMatrixAlgorithm(int j, int n, IReadOnlyList<double> a, IReadOnlyList<double> b, IReadOnlyList<double> c, IReadOnlyList<double> d)
        {
            for (var i = 0; i < b.Count; i++)
                if (Math.Abs(b[i]) < Math.Abs(a[i]) + Math.Abs(c[i]))
                    throw new Exception($"There is no diagonal dominance! j={j} i={i}: {Math.Abs(a[i])} {Math.Abs(b[i])} {Math.Abs(c[i])} sum={Math.Abs(a[i] + c[i])} ");

            if (Math.Abs(b[0]) < double.Epsilon)
                throw new InvalidOperationException("c[0] == 0");
            var y = new double[n];
            var alpha_arr = new double[n];
            var beta_arr = new double[n];
            y[0] = b[0];
            alpha_arr[0] = -c[0] / y[0];
            beta_arr[0] = d[0] / y[0];
            for (var i = 1; i < n - 1; ++i)
            {
                y[i] = b[i] + a[i] * alpha_arr[i - 1];
                alpha_arr[i] = -c[i] / y[i];
                beta_arr[i] = (d[i] - a[i] * beta_arr[i - 1]) / y[i];
            }

            y[n - 1] = b[n - 1] + a[n - 1] * alpha_arr[n - 2];
            beta_arr[n - 1] = (d[n - 1] - a[n - 1] * beta_arr[n - 2]) / y[n - 1];
            var x = new double[n];
            x[n - 1] = beta_arr[n - 1];
            for (var i = n - 2; i >= 0; --i)
                x[i] = alpha_arr[i] * x[i + 1] + beta_arr[i];
            return x;
        }

        private static double[] GetS0(double[] s0_dash)
        {
            var res = new double[s0_dash.Length];
            for (var i = 0; i < res.Length; ++i)
                res[i] = K * Math.Exp(s0_dash[i]);
            return res;
        }

        [SuppressMessage("ReSharper", "UnusedMember.Local")]
        private static void PrintLowerBound(double r, double q, double sigma2, double K)
        {
            var omega = (-r + q + sigma2 / 2d) / sigma2;
            var alphaMin = omega - Math.Sqrt(omega * omega + 2d * r / sigma2);
            var reversedAlphaMin = 1d / alphaMin;
            var s0LowerBound = K / (1d - reversedAlphaMin);
            var lnS0LowerBound = Math.Log(s0LowerBound / K);
            Console.WriteLine("r: " + r);
            Console.WriteLine("q: " + q);
            Console.WriteLine("sigma2: " + sigma2);
            Console.WriteLine("sigma2/2d: " + sigma2 / 2d);
            Console.WriteLine("omega: " + omega);
            Console.WriteLine("alpha_minus: " + alphaMin);
            Console.WriteLine("1/alpha_minus: " + reversedAlphaMin);
            Console.WriteLine("S_0 lower bound (K/(1-1/alpha_minus)): " + s0LowerBound);
            Console.WriteLine("ln(S_0/K): " + lnS0LowerBound);
            Console.WriteLine("==");
        }

        // ReSharper disable once UnusedMember.Local
        private static void S0DashDirectTime(string name, string zone, double[] arr, double in_tau, double[] in_taus, bool nonUniformTau) // значения S0 до обратного преобразования
        {
            using var writer = new StreamWriter(name, false);
            writer.WriteLine("TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA');");
            writer.WriteLine("VARIABLES = S0_dash t");
            writer.WriteLine($"ZONE T='{zone}'");
            var lines = new List<string>();
            var t = 0d;
            for (var k = 0; k < in_taus.Length; k++)
            {
                lines.Add($"{arr[k]:e16} {t:e16}");
                t += nonUniformTau ? in_taus[k] : in_tau;
            }

            var sb = new StringBuilder();
            foreach (var line in lines)
                sb.AppendLine(line);
            var i = Math.Min(arr.Length - 1, lines.Count);
            writer.WriteLine($"I={i} K={1} ZONETYPE=Ordered");
            writer.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE)");
            writer.WriteLine(sb.ToString());
        }

        private static void S0ReversedTimeE(string name, string zone, double[] arr, double T, double in_tau, double[] in_taus, bool nonUniformTau) // значения S0 после обратного преобразования
        {
            using var writer = new StreamWriter(name, false);
            writer.WriteLine("TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA');");
            writer.WriteLine("VARIABLES = t S0");
            writer.WriteLine($"ZONE T='{zone}'");

            var builder = new StringBuilder();
            int I;
            if (nonUniformTau)
            {
                var j = 0;
                var t = T;
                while (t >= 0)
                {
                    builder.AppendLine($"{t:e16} {arr[j]:e16}");
                    t -= in_taus[j];
                    j++;
                }

                I = j;
            }
            else
            {
                var t = 0d;
                foreach (var val in arr)
                {
                    builder.AppendLine($"{t:e16} {val:e16}");
                    t += in_tau;
                    if (T - t < 1e-10)
                        t = T;
                }

                I = arr.Length;
            }

            writer.WriteLine($"I={I} K={1} ZONETYPE=Ordered");
            writer.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE)");
            writer.WriteLine(builder);
        }

        private static void test(string name, double[] arr, double T, double in_tau, double[] in_taus, bool nonUniformTau) // значения S0 после обратного преобразования
        {
            using var writer = new StreamWriter(name, false);
            var builder = new StringBuilder();
            if (nonUniformTau)
            {
                var j = 0;
                var t = T;
                while (t >= 0)
                {
                    builder.AppendLine($"{j} {t:e16} {arr[j]:e16}");
                    t -= in_taus[j];
                    j++;
                }
            }
            else
            {
                var t = 0d;
                foreach (var val in arr)
                {
                    builder.AppendLine($"{t:e16} {val:e16}");
                    t += in_tau;
                    if (T - t < 1e-10)
                        t = T;
                }
            }

            writer.WriteLine(builder);
        }

        private static void S0ReversedTime(string name, string zone, double[] arr, double T, double in_tau, double[] in_taus, bool nonUniformTau) // значения S0 после обратного преобразования
        {
            using var writer = new StreamWriter(name, false);
            writer.WriteLine("TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA');");
            writer.WriteLine("VARIABLES = S0 t");
            writer.WriteLine($"ZONE T='{zone}'");

            var builder = new StringBuilder();
            int I;
            if (nonUniformTau)
            {
                var j = 0;
                var t = T;
                while (t >= 0)
                {
                    builder.AppendLine($"{arr[j]:e16} {t:e16}");
                    t -= in_taus[j];
                    j++;
                }

                I = j;
            }
            else
            {
                var t = 0d;
                foreach (var val in arr)
                {
                    builder.AppendLine($"{val:e16} {t:e16}");
                    t -= in_tau;
                    if (t < 1e-10)
                        t = 0d;
                }

                I = arr.Length;
            }

            writer.WriteLine($"I={I} K={1} ZONETYPE=Ordered");
            writer.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE)");
            writer.WriteLine(builder);
        }

        private static double[] GetTaus()
        {
            var res = new double[M + 1];
            for (var j = 1; j < res.Length; j++)
                res[j - 1] = gt2(j * tau0, alpha) - gt2((j - 1) * tau0, alpha);
            return res;
        }

        private static double[] GetHs()
        {
            var res = new double[N1 + 1];
            for (var i = 0; i < res.Length - 1; i++)
                res[i + 1] = gh2((i + 1) * h, beta) - gh2(i * h, beta);
            return res;
        }

        private static bool CheckNaN(double val, int l, int j)
        {
            if (!double.IsNaN(val))
                return false;

            print_to_log($"rho is NaN! j = {j} l = {l}");
            return true;
        }

        private static void CheckS0Dash(double[] arr)
        {
            if (arr[0] != 0d)
                throw new Exception("arr[0] is not 0!");
            for (var i = 2; i < arr.Length; i++)
            {
                if (Math.Abs(arr[i - 1]) > Math.Abs(arr[i]))
                    throw new Exception($"Math.Abs(arr[{i - 1}]={arr[i - 1]}) is greater than Math.Abs(arr[{i}]={arr[i]})!");
            }
        }

        private static void CheckS0(double[] arr)
        {
            for (var i = 1; i < arr.Length; i++)
            {
                if (arr[i - 1] < arr[i])
                    throw new Exception($"arr[{i - 1}]={arr[i - 1]} is less than arr[{i}]={arr[i]}!");
            }
        }

        [SuppressMessage("ReSharper", "UnusedMember.Local")]
        private static void PrintCondensedMeshesTau()
        {
            Console.WriteLine($"M={M} Tn={Tn} T0={T0} T={T} alpha={alpha} tau0={tau0}");
            Console.WriteLine("TAU V1" + new string(' ', 10) + "TAU V2");
            var t = 0d;
            for (var j = 1; j < 110; j++)
            {
                var t1 = gt1(j * tau0, alpha) - gt1((j - 1) * tau0, alpha);
                var t2 = taus[j - 1];
                t += t2;
                Console.Write("j: " + j + " " + t1.ToString("F8") + " " + t2.ToString("F8") + " " + t.ToString("F8") + (t1.Equals(t2) ? '=' : '\0'));
                Console.WriteLine();
            }
        }

        [SuppressMessage("ReSharper", "UnusedMember.Local")]
        private static void PrintCondensedMeshesH()
        {
            Console.WriteLine($"N={N} beta={beta} h={h}");
            Console.WriteLine("H V1" + new string(' ', 12) + "H V2");
            var x = 0d;
            for (var i = 1; i < 100; i++)
            {
                var x1 = gh1(i * h, beta) - gh1((i - 1) * h, beta);
                var x2 = hs[i - 1];
                x += x2;
                Console.Write("i: " + i + " " + x1.ToString("F8") + " " + x2.ToString("F8") + " " + x.ToString("F8") + (x1.Equals(x2) ? '=' : '\0'));
                Console.WriteLine();
            }
        }

        // ReSharper disable once UnusedMember.Local
        private static void CreateLay(string layTemplate, string output, string datFile)
        {
            var content = File.ReadAllText(layTemplate);
            var replacedContent = content.Replace("{NAME}", datFile);
            File.WriteAllText(output, replacedContent);
        }
    }
}