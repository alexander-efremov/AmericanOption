using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Text;

#pragma warning disable 219
namespace AmericanOptionAlbena
{
    [SuppressMessage("ReSharper", "ConvertToCompoundAssignment")]
    [SuppressMessage("ReSharper", "SuggestBaseTypeForParameter")]
    [SuppressMessage("ReSharper", "ArrangeRedundantParentheses")]
    public static class Program
    {
        [SuppressMessage("ReSharper", "RedundantAssignment")]
        [SuppressMessage("ReSharper", "TooWideLocalVariableScope")]
        [SuppressMessage("ReSharper", "RedundantIfElseBlock")]
        public static void Main()
        {
            const int l_max_iterations = 1000; // max iteration on l
            const double Tol = 10e-10; // eps to refine eta
            const double eps = 10e-3; // eps to refine mu
            const double lb = 0d; // left bound
            const double rb = 100d; // right bound
            const double T0 = 0d; // the start time
            const double Tn = 60d; // the finish time
            const double T = Tn - T0; // time interval
            const double sigma = 0.1d; // the sigma = volatility
            const double sigma2 = sigma * sigma; // the squared sigma
            const double r = 0.1d; // the risk-free rate
            const double K = 10; // the strike price
            const double q = 0d; // the dividend rate
            const int N = 800; // the number of space intervals
            const int N1 = N + 1; // the number of points
            const double h = (rb - lb) / N; // the space step
            const int M = 10000; // the number of time intervals
            const double tau = T / M; // the time step
            var u_curr = new double[N1]; // the current solution
            var u_prev = new double[N1]; // the previous solution
            var s0_hat_deriv = new double[M + 1];
            var s0_wave = new double[M + 1];
            var lambda = new double[M + 1];
            var a = new double[M + 1];
            // var sols = new List<double[]>();
            var print = true;

            for (var j = 1; j <= M; j++)
            {
                var mu_j = new double[l_max_iterations];
                var eta_j = new double[l_max_iterations];
                var rho_j_l = 0d;
                var l = 0;
                mu_j[l] = s0_hat_deriv[j - 1];
                if (print)
                {
                    Console.WriteLine("l = {0} eta_j[l-2] = {1}", l, l > 1 ? eta_j[l - 2] : 0d);
                    Console.WriteLine("l = {0} eta_j[l-1] = {1}", l, l > 0 ? eta_j[l - 1] : 0d);
                    Console.WriteLine("l = {0} eta_j[l] = {1}", l, eta_j[l]);
                    Console.WriteLine("l = {0} mu_j[l-2] = {1}", l, 0d);
                    Console.WriteLine("l = {0} mu_j[l-1] = {1}", l, 0d);
                    Console.WriteLine("l = {0} mu_j[l] = {1}", l, mu_j[l]);
                    Console.WriteLine("----------------------------------");
                }

                while (l < l_max_iterations)
                {
                    a[j] = (2d * r - 2d * q - sigma2 + 2d * mu_j[l]) / sigma2;
                    var b = (-2d * r) / sigma2;
                    lambda[j] = Math.Sqrt(a[j] * a[j] - 4d * b);
                    Console.WriteLine("l = {0} j = {1} s0_wave[j] = {2}", l, j, rho_j_l);
                    u_curr = Solve(j, N1, r, q, tau, h, u_prev, sigma2, rho_j_l, s0_wave, mu_j[l], a[j], lambda);
                    //sols.Add(u_curr);
                    rho_j_l = s0_wave[j - 1] + tau * mu_j[l];

                    eta_j[l] = K - rho_j_l - u_curr[0];
                    //eta_j[l] = K * (1d - Math.Exp(rho_j_l)) - u_curr[0];
                    if (Math.Abs(eta_j[l]) <= Tol)
                    {
                        s0_hat_deriv[j] = mu_j[l];
                        s0_wave[j] = rho_j_l;
                        break;
                    }

                    if (++l == 1)
                    {
                        mu_j[l] = mu_j[0] + eps;
                        if (print)
                        {
                            Console.WriteLine("l = {0} eta_j[l-2] = {1}", l, 0d);
                            Console.WriteLine("l = {0} eta_j[l-1] = {1}", l, l > 0 ? eta_j[l - 1] : 0d);
                            Console.WriteLine("l = {0} eta_j[l] = {1}", l, eta_j[l]);
                            Console.WriteLine("l = {0} mu_j[l-2] = {1}", l, 0d);
                            Console.WriteLine("l = {0} mu_j[l-1] = {1}", l, mu_j[l - 1]);
                            Console.WriteLine("l = {0} mu_j[l] = {1}", l, mu_j[l]);
                            Console.WriteLine("----------------------------------");
                        }
                    }
                    else
                    {
                        if (print)
                        {
                            Console.WriteLine("l = {0} mu_j[l-2] = {1}", l, mu_j[l - 2]);
                            Console.WriteLine("l = {0} mu_j[l-1] = {1}", l, mu_j[l - 1]);
                            Console.WriteLine("l = {0} mu_j[l] = {1}", l, mu_j[l]);
                            Console.WriteLine("l = {0} eta_j[l-2] = {1}", l, eta_j[l - 2]);
                            Console.WriteLine("l = {0} eta_j[l-1] = {1}", l, eta_j[l - 1]);
                        }

                        var value = (eta_j[l - 1] * mu_j[l - 2] - eta_j[l - 2] * mu_j[l - 1]) /
                                    (eta_j[l - 1] - eta_j[l - 2]);
                        mu_j[l] = value;

                        if (print)
                        {
                            Console.WriteLine(
                                "l = {0} (eta_j[l - 1] * mu_j[l - 2] - eta_j[l - 2] * mu_j[l - 1]) = {1:g6}",
                                l, eta_j[l - 1] * mu_j[l - 2] - eta_j[l - 2] * mu_j[l - 1]);
                            Console.WriteLine("l = {0} (eta_j[l - 1] - eta_j[l - 2]) = {1:g6}", l,
                                eta_j[l - 1] - eta_j[l - 2]);
                            Console.WriteLine("l = {0} mu_j[l-2] = {1}", l, mu_j[l - 2]);
                            Console.WriteLine("l = {0} mu_j[l-1] = {1}", l, mu_j[l - 1]);
                            Console.WriteLine("l = {0} calculated mu_j[l] = {1}", l, mu_j[l]);
                            Console.WriteLine("----------------------------------");
                        }
                    }
                }

                if (l >= l_max_iterations)
                    throw new InvalidOperationException("The 'l' exceeded the maximum number of iterations.");

                for (var i = 0; i < u_prev.Length; i++)
                    u_prev[i] = u_curr[i];

                // if (Math.Abs(s0_wave[j] - s0_wave[j-1]) < Tol)
                //     break;
            }

            // var printer = new TecplotPrinter(K, rb);
            // for (var i = 0; i < sols.Count; i++)
            // {
            //     var frame = "V(S,t) on t=";
            //     if (i == 0)
            //         frame += "T";
            //     else
            //         frame += sols.Count - i + "(T/" + sols.Count + ")";
            //
            //     printer.PrintXY(
            //         Path.Combine(Directory.GetCurrentDirectory(), "numeric/numericSolution"),
            //         i,
            //         (tau * (M - i)).ToString("F8"),
            //         h.ToString("F8"),
            //         h,
            //         K,
            //         sols[i],
            //         "numeric_" + (tau * (M - i)).ToString("F8"),
            //         frame,
            //         i == 0);
            // }

            for (var i = 2; i < s0_wave.Length; i++)
                if (s0_wave[i - 1] > s0_wave[i])
                    throw new InvalidOperationException("The s0 is not convex!");

            WriteVectorToFile("S0.dat", s0_wave);
            PrintToTecplot($"{nameof(s0_wave)}_nx={s0_wave.Length - 1}_tau={tau}.dat", s0_wave, tau);
            PrintToTecplotT($"{nameof(s0_wave)}_reversed_nx={s0_wave.Length - 1}_tau={tau}.dat", s0_wave, tau, T);

            // print s0 to tecplot file
            var s0 = new double[s0_wave.Length];
            for (var i = 0; i < s0_wave.Length; i++)
                s0[i] = K * Math.Exp(s0_wave[i]);
            // PrintToTecplotT($"{nameof(s0)}_T-t_nx={s0.Length}_tau={tau}.dat", s0, tau, T);
            // PrintToTecplotT($"{nameof(s0)}_T-t_rev_nx={s0.Length}_tau={tau}.dat", s0.Reverse().ToArray(), tau, T);
            // PrintToTecplot($"{nameof(s0)}_t_nx={s0.Length}_tau={tau}.dat", s0, tau);
            // PrintToTecplot($"{nameof(s0)}_t_rev_nx={s0.Length}_tau={tau}.dat", s0.Reverse().ToArray(), tau);
        }

        private static void PrintToTecplot(string name, double[] s0, double tau)
        {
            using var writer2 = new StreamWriter(name, false);
            writer2.WriteLine(
                "TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA'\nVARIABLES = S_0 {0}\nZONE T='{1}'",
                "t_j", "SubZone");
            writer2.WriteLine($"I={s0.Length} K={1} ZONETYPE=Ordered");
            writer2.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE)");
            for (var i = 0; i < s0.Length; i++)
                writer2.WriteLine("{0:e8} {1:e8}", s0[i], tau * i);
        }

        private static void PrintToTecplotT(string name, double[] s0, double tau, double T)
        {
            using var writer2 = new StreamWriter(name, false);
            writer2.WriteLine(
                "TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA'\nVARIABLES = S {0}\nZONE T='{1}'",
                "T - t", "SubZone");
            writer2.WriteLine($"I={s0.Length - 1} K={1} ZONETYPE=Ordered");
            writer2.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE)");
            for (var i = 0; i < s0.Length; i++)
                writer2.WriteLine("{0:e8} {1:e8}", s0[i], T - tau * i);
        }

        private static double[] Solve(
            int j,
            int N1,
            double r,
            double q,
            double tau,
            double h,
            double[] u_prev,
            double sigma2,
            double rho_j_l, // to toke l-th approximation of s_dash_derivative (t_j)
            double[] s0_wave, // to toke l-th approximation of s_dash_derivative (t_j)
            double mu_j_l,
            double a_j,
            double[] lambda)
        {
            // calculate the right part
            var alpha = r - q - (sigma2 / 2d) + mu_j_l;
            Debug.Assert(alpha >= 0, "alpha>=0");
            var gamma1 = (1d / tau) - (alpha / h);
            var gamma2 = (alpha / h);
            var f = new double[N1];
            f[0] = gamma1 * u_prev[0] + gamma2 * u_prev[1] - (sigma2 * rho_j_l) / h;
            for (var i = 1; i < N1 - 1; i++)
                f[i] = gamma1 * u_prev[i] + gamma2 * u_prev[i + 1];
            var nu_j = 2d / (lambda[j - 1] + lambda[j] + (2d * (s0_wave[j - 1] - mu_j_l)) / sigma2);
            var u_prev_n_plus_1 = 0d; // u_prev[N1]; TODO: How to calculate the value at the index N1
            f[N1 - 1] = 0.5d * (gamma1 * u_prev[N1 - 1] 
                             + gamma2 * u_prev_n_plus_1) 
                             + (nu_j * u_prev[N1 - 1]) / (tau * h);

            // if (l == 1)
            //     Console.WriteLine();
            // double s = 0;
            // for (var i = 0; i < f.Length; i++)
            //     s += f[i];
            // Console.WriteLine(s);
            // Utils.Print(f, "f");

            // b - below main diagonal (indexed as [1;n-1])
            var a = new double[N1];
            a[0] = 0d;
            for (var i = 1; i < N1 - 1; ++i)
                a[i] = -sigma2 / (2d * h * h);

            a[N1 - 1] = -sigma2 / (2d * h * h); // right boundary condition

            // main diagonal of matrix (indexed as [0;n-1])
            var b = new double[N1];
            b[0] = sigma2 / (h * h) + r + (1d / tau) + (-sigma2 / (2d * h * h));
            for (var i = 1; i < N1 - 1; ++i)
                b[i] = sigma2 / (h * h) + r + (1d / tau);

            var val1 = sigma2 / (2d * h * h);
            var val2 = sigma2 * (a_j + lambda[j - 1]) / (4d * h);
            var val3 = 0.5d * ((1d / tau) + r);
            var val4 = 1d / (lambda[j - 1] * tau * h);
            b[N1 - 1] = val1 + val2 + val3 + val4; // right boundary condition

            // d - up to main diagonal (indexed as [0;n-2])
            var c = new double[N1];
            c[0] = -sigma2 / (h * h);
            for (var i = 1; i < N1 - 2; ++i)
                c[i] = -sigma2 / (2d * h * h);

            // right boundary condition
            c[N1 - 2] = 0d;
            // Utils.Print(db, "db");
            // Utils.Print(dc, "dc");
            // Utils.Print(dd, "dd");
            var u = SolveByTridiagonalMatrixAlgorithm(N1, a, b, c, f);
            // Utils.Print(u_curr, "u_curr");
            return u;
        }

        private static double[] SolveByTridiagonalMatrixAlgorithm(int n, double[] a, double[] b, double[] c, double[] d)
        {
            // https://pro-prof.com/forums/topic/sweep-method-for-solving-systems-of-linear-algebraic-equations
            for (var i = 0; i < b.Length; i++)
                if (Math.Abs(b[i]) < Math.Abs(a[i]) + Math.Abs(c[i]))
                    throw new Exception(
                        $"There is no diagonal dominance! i={i} {Math.Abs(a[i])} {Math.Abs(b[i])} {Math.Abs(c[i])} sum={Math.Abs(a[i] + c[i])} ");
            if (Math.Abs(b[0]) < double.Epsilon)
                throw new InvalidOperationException("c[0] == 0");
            var y = new double[n];
            var alpha = new double[n];
            var beta = new double[n];
            y[0] = b[0];
            alpha[0] = -c[0] / y[0];
            beta[0] = d[0] / y[0];
            for (var i = 1; i < n - 1; ++i)
            {
                y[i] = b[i] + a[i] * alpha[i - 1];
                alpha[i] = -c[i] / y[i];
                beta[i] = (d[i] - a[i] * beta[i - 1]) / y[i];
            }

            y[n - 1] = b[n - 1] + a[n - 1] * alpha[n - 2];
            beta[n - 1] = (d[n - 1] - a[n - 1] * beta[n - 2]) / y[n - 1];
            var x = new double[n];
            x[n - 1] = beta[n - 1];
            for (var i = n - 2; i >= 0; --i)
                x[i] = alpha[i] * x[i + 1] + beta[i];
            return x;
        }

        private static void WriteVectorToFile(string path, IList<double> arr)
        {
            var sb = new StringBuilder();
            foreach (var t in arr)
                sb.AppendLine(t.ToString("e10"));
            File.WriteAllText(path, sb.ToString());
        }
    }
}