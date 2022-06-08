using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Text;

#pragma warning disable 219
namespace AmericanOptionAlbena
{
    [SuppressMessage("ReSharper", "InconsistentNaming")]
    [SuppressMessage("ReSharper", "PossibleNullReferenceException")]
    [SuppressMessage("ReSharper", "CommentTypo")]
    [SuppressMessage("ReSharper", "ParameterHidesMember")]
    public static class Program
    {
        private const int l_max_iterations = 1000; // max iteration on l
        private const double Tol = 10e-10; // eps to refine eta
        private const double eps = 10e-3; // eps to refine mu
        private const double lb = 0d; // left bound
        private const double rb = 100d; // right bound
        private const double T0 = 0d; // the start time
        private const double Tn = 60d; // the finish time
        private const double T = Tn - T0; // time interval
        private const double sigma = 0.1d; // the sigma = volatility
        private const double sigma2 = sigma * sigma; // the squared sigma
        private const double r = 0.1d; // the risk-free rate
        private const double K = 10; // the strike price
        private const double q = 0.01d; // the dividend rate
        private const int N = 1000; // the number of space intervals
        private const int N1 = N + 1; // the number of points
        private const int N0 = (int)(N / 10d)*3; // the number of points before h0
        private const double h = (rb - lb) / N; // the space step
        private const int M = 10000; // the number of time intervals
        private const double tau = T / M; // the time step
        private const double alpha = 1.1d; // condense parameter for t
        private const double beta = 1.3d; // condense parameter for h

        public static void Main()
        {
            var u_curr = new double[N1]; // the current solution
            var u_prev = new double[N1]; // the previous solution
            var s0_hat_deriv = new double[M + 1];
            var s0_wave = new double[M + 1];
            var lambda = new double[M + 1];
            var a = new double[M + 1];
            var print = true;

            // s0_hat_deriv[0] = s0_wave[0] = -1; // failed on condition var alpha = r - q - (sigma2 / 2d) + mu_j_l; Debug.Assert(alpha >= 0, "alpha>=0");
            // s0_hat_deriv[j] = s0_wave[j] = -sigma2 / 2d;
            const double b = -2d * r / sigma2;
            s0_hat_deriv[0] = -eps;
            s0_wave[0] = 0d;
            a[0] = (2d * r - 2d * q - sigma2 + 2d * s0_hat_deriv[0]) /
                   sigma2; // here s0_hat_deriv[j] == mu_j[l] where l = 0
            lambda[0] = Math.Sqrt(a[0] * a[0] - 4d * b);
            for (var j = 1; j <= M; j++)
            {
                var mu_j = new double[l_max_iterations];
                var eta_j = new double[l_max_iterations];
                var l = 0;
                mu_j[l] = s0_hat_deriv[j - 1];
                while (l < l_max_iterations)
                {
                    a[j] = (2d * r - 2d * q - sigma2 + 2d * mu_j[l]) / sigma2;
                    lambda[j] = Math.Sqrt(a[j] * a[j] - 4d * b);
                    var rho_j_l = s0_wave[j - 1] + tau * mu_j[l];
                    if (double.IsNaN(rho_j_l))
                    {
                        Console.WriteLine("rho is NaN!");
                        break;
                    }

                    u_curr = Solve(j, u_prev, rho_j_l, s0_hat_deriv, mu_j[l], a, lambda);
                    // eta_j[l] = K * (1d - Math.Exp(rho_j_l)) - u_curr[0];
                    eta_j[l] = K * Math.Exp(rho_j_l) + ((u_curr[1] - u_curr[0]) / h);
                    Console.WriteLine("l = {0} j = {1} mu_j_l {7} K {2} rho_j_l {3} u[1] {4} u[0] {5}  eta_j_l = {6}",
                        l, j, K, rho_j_l, u_curr[1], u_curr[0], eta_j[l], mu_j[l]);
                    if (Math.Abs(eta_j[l]) <= Tol)
                    {
                        s0_hat_deriv[j] = mu_j[l];
                        s0_wave[j] = rho_j_l;
                        break;
                    }

                    if (++l == 1)
                        mu_j[l] = mu_j[l - 1] - eps;
                    else
                        mu_j[l] = (eta_j[l - 1] * mu_j[l - 2] - eta_j[l - 2] * mu_j[l - 1]) /
                                  (eta_j[l - 1] - eta_j[l - 2]);
                }

                WriteVectorToFile("eta_j1_l.arr", eta_j);
                Console.WriteLine("===========");

                if (l >= l_max_iterations)
                    throw new InvalidOperationException("The 'l' exceeded the maximum number of iterations.");

                for (var i = 0; i < u_prev.Length; i++)
                    u_prev[i] = u_curr[i];
            }

            WriteVectorToFile("S0.arr", s0_wave);
            WriteVectorToFile("S0_hat_deriv.arr", s0_hat_deriv);
            // значения производных
            PrintS0HatDirectTime($"{nameof(s0_hat_deriv)}_t_nx={s0_hat_deriv.Length}_tau={tau}.dat", s0_hat_deriv, tau);
            // значения S0 до обратного преобразования
            PrintS0WaveDirectTime($"{nameof(s0_wave)}_t_nx={s0_wave.Length}_tau={tau}.dat", s0_wave, tau);
            // print s0 to tecplot file
            var s0 = new double[s0_wave.Length];
            for (var i = 0; i < s0_wave.Length; i++)
                s0[i] = K * Math.Exp(s0_wave[i]);
            // значения S0 после обратного преобразования
            PrintS0ReversedTime($"{nameof(s0)}_T-t_nx={s0.Length}_tau={tau}.dat", s0, tau, T);
        }

        private static double[] Solve(int j, double[] u_prev, double rho_j_l, double[] s0_hat_deriv, double mu_j_l,
            double[] arr_a, double[] lambda)
        {
            // calculate the right part
            // var alpha = r - q - (sigma2 / 2d) + s0_wave[j - 1];
            var f = new double[N1];
            var alpha_jm1 = r - q - (sigma2 / 2d) + s0_hat_deriv[j - 1];
            var gamma1_jm1 = (1d / tau) - alpha_jm1 / h;
            var gamma2_jm1 = (alpha_jm1 / h);
            f[0] = K * (1d - Math.Exp(rho_j_l));
            for (var i = 1; i < N1 - 1; i++)
            {
                if (i < N0)
                {
                    var gamma1_jm1_0 = (1d / tau) - alpha_jm1 / GetH(i, h, beta);
                    // todo: возможно тут h_{i+1}, то есть надо GetH(i + 1, h, beta)
                    // todo: но тогда метод разваливается, rho = NaN
                    var gamma2_jm1_0 = (alpha_jm1 / GetH(i, h, beta)); 
                    f[i] = gamma1_jm1_0 * u_prev[i] + gamma2_jm1_0 * u_prev[i + 1];
                }
                else
                {
                    f[i] = gamma1_jm1 * u_prev[i] + gamma2_jm1 * u_prev[i + 1];
                }
            }

            var nu_j = 2d / (lambda[j - 1] + lambda[j] + (2d * (s0_hat_deriv[j - 1] - mu_j_l)) / sigma2);

            // formula 31 is applied to j-1 time level
            // here /*x_{n+1}-x_n= x_n + h - x_n = h*/
            var u_prev_np1 = u_prev[N1 - 1] * Math.Exp(-((lambda[j - 1] + arr_a[j - 1]) * h) / 2d);
            f[N1 - 1] = 0.5d * (gamma1_jm1 * u_prev[N1 - 1] + gamma2_jm1 * u_prev_np1) +
                        (nu_j * u_prev[N1 - 1]) / (tau * h);

            // if (l == 1)
            //     Console.WriteLine();
            // double s = 0;
            // for (var i = 0; i < f.Length; i++)
            //     s += f[i];
            // Console.WriteLine(s);
            // Utils.Print(f, "f");

            // a - below main diagonal (indexed as [1;n-1])
            var a0 = new double[N1];
            a0[0] = 0d;
            for (var i = 1; i < N1 - 1; ++i)
            {
                if (i < N0)
                {
                    a0[i] = -sigma2 / (2d * GetH(i, h, beta) * GetH(i, h, beta));
                }
                else
                {
                    a0[i] = -sigma2 / (2d * h * h);
                }
            }

            a0[N1 - 1] = -sigma2 / (2d * h * h);

            // main diagonal of matrix (indexed as [0;n-1])
            var b0 = new double[N1];
            // b0[0] = sigma2 / (h * h) + r + (1d / tau);
            b0[0] = 1d;
            for (var i = 1; i < N1 - 1; ++i)
            {
                if (i < N0)
                {
                    b0[i] = sigma2 / (GetH(i, h, beta) * GetH(i, h, beta)) + r + (1d / tau);
                }
                else
                {
                    b0[i] = sigma2 / (h * h) + r + (1d / tau);
                }
            }

            var val1 = sigma2 / (2d * h * h);
            var val2 = (sigma2 * (arr_a[j] + lambda[j])) / (4d * h);
            var val3 = 0.5d * ((1d / tau) + r);
            var val4 = 1d / (lambda[j] * tau * h);
            b0[N1 - 1] = val1 + val2 + val3 + val4;

            // c - up to main diagonal (indexed as [0;n-2])
            var c0 = new double[N1];
            // c0[0] = -sigma2 / (h * h);
            c0[0] = 0d;
            for (var i = 1; i < N1 - 2; ++i)
            {
                if (i < N0)
                {
                    c0[i] = -sigma2 / (2d * GetH(i, h, beta) * GetH(i, h, beta));
                }
                else
                {
                    c0[i] = -sigma2 / (2d * h * h);
                }
            }

            c0[N1 - 2] = 0d;

            // Utils.Print(db, "db");
            // Utils.Print(dc, "dc");
            // Utils.Print(dd, "dd");
            var u = SolveByTridiagonalMatrixAlgorithm(N1, a0, b0, c0, f);
            // Utils.Print(u_curr, "u_curr");
            return u;
        }

        private static double GetH(int i, double h, double in_beta) => Math.Pow(i * h, in_beta);

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

        [SuppressMessage("ReSharper", "UnusedMember.Local")]
        private static void PrintLowerBound(double r, double q, double sigma2, double K)
        {
            var omega = (-r + q + sigma2 / 2d) / (sigma2);
            var alpha_min = omega - Math.Sqrt(omega * omega + ((2d * r) / sigma2));
            var reversed_alpha_min = 1d / alpha_min;
            var S0LowerBound = K / (1d - reversed_alpha_min);
            var lnS0LowerBound = Math.Log(S0LowerBound / K);
            Console.WriteLine("r: " + r);
            Console.WriteLine("q: " + q);
            Console.WriteLine("sigma2: " + sigma2);
            Console.WriteLine("sigma2/2d: " + sigma2 / 2d);
            Console.WriteLine("omega: " + omega);
            Console.WriteLine("alpha_minus: " + alpha_min);
            Console.WriteLine("1/alpha_minus: " + reversed_alpha_min);
            Console.WriteLine("S_0 lower bound (K/(1-1/alpha_minus)): " + S0LowerBound);
            Console.WriteLine("ln(S_0/K): " + lnS0LowerBound);
            Console.WriteLine("==");
        }

        private static void PrintS0ReversedTime(string name, double[] arr, double tau, double T)
        {
            using var writer = new StreamWriter(name!, false);
            writer.WriteLine("TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA');");
            writer.WriteLine("VARIABLES = t S0");
            writer.WriteLine("ZONE T='SubZone'");
            writer.WriteLine($"I={arr.Length - 1} K={1} ZONETYPE=Ordered");
            writer.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE)");
            for (var i = 0; i < arr.Length; i++)
                writer.WriteLine("{0:e12} {1:e12}", T - tau * i, arr[i]);
        }

        private static void PrintS0WaveDirectTime(string name, double[] arr, double tau)
        {
            using var writer = new StreamWriter(name!, false);
            writer.WriteLine("TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA');");
            writer.WriteLine("VARIABLES = t S0W");
            writer.WriteLine("ZONE T='SubZone'");
            writer.WriteLine($"I={arr.Length - 1} K={1} ZONETYPE=Ordered");
            writer.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE)");
            for (var i = 0; i < arr.Length; i++)
                writer.WriteLine("{0:e12} {1:e12}", tau * i, arr[i]);
        }

        private static void PrintS0HatDirectTime(string name, double[] arr, double tau)
        {
            using var writer = new StreamWriter(name!, false);
            writer.WriteLine("TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA');");
            writer.WriteLine("VARIABLES = t S0_Deriv");
            writer.WriteLine("ZONE T='SubZone'");
            writer.WriteLine($"I={arr.Length - 1} K={1} ZONETYPE=Ordered");
            writer.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE)");
            for (var i = 0; i < arr.Length; i++)
                writer.WriteLine("{0:e12} {1:e12}", tau * i, arr[i]);
        }

        [SuppressMessage("ReSharper", "UnusedMember.Local")]
        private static void PrintToTecplotT(string name, double[] s0, double tau, double T)
        {
            using var writer = new StreamWriter(name!, false);
            writer.WriteLine("TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA'\nVARIABLES = S {0}\nZONE T='{1}'",
                "T - t", "SubZone");
            writer.WriteLine($"I={s0.Length - 1} K={1} ZONETYPE=Ordered");
            writer.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE)");
            for (var i = 0; i < s0.Length; i++)
                writer.WriteLine("{0:e12} {1:e12}", s0[i], T - tau * i);
        }

        private static void WriteVectorToFile(string path, IEnumerable<double> arr)
        {
            var builder = new StringBuilder();
            foreach (var a in arr!)
                builder.AppendLine(a.ToString("e10"));
            File.WriteAllText(path, builder.ToString());
        }
    }
}