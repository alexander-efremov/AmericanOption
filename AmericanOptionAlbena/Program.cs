using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Text;
using CoreLib;

#pragma warning disable 219
#pragma warning disable CS0162
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
        private const double Tn = 10d; // the finish time
        private const double T = Tn - T0; // time interval
        private const double sigma2 = sigma * sigma; // the squared sigma
        private const double sigma = 0.1d; // the sigma = volatility
        private const double r = 0.1d; // the risk-free rate
        private const double K = 10d; // the strike price
        private const double q = 0.01d; // the dividend rate
        private const int N = 1000; // the number of space intervals
        private const int N1 = N + 1; // the number of points

        private const double h = (rb - lb) / N; // the space step

        // private const int M = 10000; // the number of time intervals
        private const int M = 1000; // the number of time intervals
        private const double tau0 = T / M; // the time step
        private const double alpha = 2d; // condense parameter for t
        private const double beta = 2d; // condense parameter for h
        private const double beta_mu = 0.5d; // the parameter to update mu

        private const bool print = true; // the parameter which allows enable/disable printing of the results

        // private const bool is_time_condensed = true; // the parameter which allows enable/disable condensed time mesh
        // private const bool is_space_condensed = true; // the parameter which allows enable/disable condensed space mesh
        private static readonly bool is_condensed_h = false; // the parameter which allows enable/disable condensed meshes
        private static readonly bool is_condensed_tau = true; // the parameter which allows enable/disable condensed meshes

        private static readonly double[] taus = GetTaus();
        private static readonly double[] hs = GetHs();

        public static void Main()
        {
            if (File.Exists("log.txt"))
                File.Delete("log.txt");
            Debug.Assert(q < r);
            // PrintCondensedMeshesTau();
            // PrintCondensedMeshesH();
            // return;
            // PrintCondensedMeshes1();
            // PrintCondensedMeshes2();
            // return;
            var u_curr = new double[N1]; // the current solution
            var u_prev = new double[N1]; // the previous solution
            var s0_hat_deriv = new double[M + 1];
            var s0_wave = new double[M + 1];
            var lambda = new double[M + 1];
            var a = new double[M + 1];
            // s0_hat_deriv[0] = s0_wave[0] = -1; // failed on condition var alpha = r - q - (sigma2 / 2d) + mu_j_l; Debug.Assert(alpha >= 0, "alpha>=0");
            // s0_hat_deriv[j] = s0_wave[j] = -sigma2 / 2d;
            const double b = -2d * r / sigma2;
            s0_hat_deriv[0] = 0d;
            s0_wave[0] = 0d;
            a[0] = (2d * r - 2d * q - sigma2 + 2d * s0_hat_deriv[0]) / sigma2; // here s0_hat_deriv[j] == mu_j[l] where l = 0
            lambda[0] = Math.Sqrt(a[0] * a[0] - 4d * b);

            for (var j = 1; j <= M; j++)
            {
                var tau = is_condensed_tau ? taus[j - 1] : tau0;
                var mu_j = new double[l_max_iterations];
                var eta_j = new double[l_max_iterations];
                var l = 0;
                mu_j[l] = s0_hat_deriv[j - 1];
                while (l < l_max_iterations)
                {
                    a[j] = (2d * r - 2d * q - sigma2 + 2d * mu_j[l]) / sigma2;
                    lambda[j] = Math.Sqrt(a[j] * a[j] - 4d * b);
                    var rho_j_l = s0_wave[j - 1] + tau * mu_j[l];
                    CheckNaN(rho_j_l);
                    u_curr = Solve(j, u_prev, rho_j_l, s0_hat_deriv, mu_j[l], a, lambda, tau);
                    var h0 = is_condensed_h ? hs[0] : h;
                    eta_j[l] = K * Math.Exp(rho_j_l) + (u_curr[1] - u_curr[0]) / h0;
                    if (print)
                    {
                        var value = string.Format("l = {0:G10} j = {1:G10} K {2:N10} mu_j_l {7:G10} rho_j_l {3:N10} u[1] {4:G10} u[0] {5:G10} eta_j_l = {6:G10}", l, j, K, rho_j_l, u_curr[1],
                            u_curr[0], eta_j[l], mu_j[l]);
                        File.AppendAllText("log.txt", value + Environment.NewLine);
                        // Console.WriteLine(value);
                    }

                    if (Math.Abs(eta_j[l]) <= Tol)
                    {
                        s0_hat_deriv[j] = mu_j[l];
                        s0_wave[j] = rho_j_l;
                        break;
                    }

                    if (++l == 1)
                        mu_j[l] = beta_mu * mu_j[l - 1] - eps;
                    else
                        mu_j[l] = (eta_j[l - 1] * mu_j[l - 2] - eta_j[l - 2] * mu_j[l - 1]) / (eta_j[l - 1] - eta_j[l - 2]);
                }

                for (var i = 0; i < u_prev.Length; i++)
                    u_prev[i] = u_curr[i];
            }

            WriteVectorToFile("S0.arr", s0_wave);
            WriteVectorToFile("S0_hat_deriv.arr", s0_hat_deriv);
            var s0 = GetS0(s0_wave);
            S0HatDirectTime($"{(is_condensed_h ? "1" : "11")}_{nameof(s0_hat_deriv)}_N1={N1}_T={T}_h_condensed_{is_condensed_h}_tau_condensed_{is_condensed_tau}.dat", s0_hat_deriv, tau0);
            S0WaveDirectTime($"{(is_condensed_h ? "2" : "22")}_{nameof(s0_wave)}_N1={N1}_T={T}_h_condensed_{is_condensed_h}_tau_condensed_{is_condensed_tau}.dat", s0_wave, tau0);
            S0ReversedTime($"{(is_condensed_h ? "3" : "33")}_{nameof(s0)}_T-t_N1={N1}_T={T}_h_condensed_{is_condensed_h}_tau_condensed_{is_condensed_tau}.dat", s0, T, tau0);
        }

        private static double[] Solve(int j, double[] u_prev, double rho_j_l, double[] s0_hat_deriv, double mu_j_l, double[] a, double[] lambda, double in_tau)
        {
            var alpha_tj = r - q - sigma2 / 2d + s0_hat_deriv[j];
            var f = new double[N1];
            var h_12 = is_condensed_h ? hs[0] : h;
            var h_32 = is_condensed_h ? hs[1] : h;
            var u_prev_1 = u_prev[1];
            var u_0j = K * (1d - Math.Exp(rho_j_l));
            f[0] = u_prev_1 / in_tau + (sigma2 / (h_12 * (h_12 + h_32))) * u_0j;
            for (var i = 1; i < N1 - 1; i++)
                f[i] = u_prev[i] / in_tau;
            var h_Nm12 = is_condensed_h ? hs[N1 - 1] : h;
            var nu_j = 2d / (lambda[j - 1] + lambda[j] + 2d * (s0_hat_deriv[j - 1] - mu_j_l) / sigma2);
            f[N1 - 1] = u_prev[N1 - 1] / (2d * in_tau) + nu_j * u_prev[N1 - 1] / (in_tau * h_Nm12);

            var a0 = new double[N1]; // a - below main diagonal (indexed as [1;n-1])
            var b0 = new double[N1]; // b - main diagonal of matrix (indexed as [0;n-1])
            var c0 = new double[N1]; // c - up to main diagonal (indexed as [0;n-2])
            a0[0] = 0d;
            for (var i = 1; i < N1 - 1; ++i)
                if (is_condensed_h)
                    a0[i] = -sigma2 / (hs[i - 1] * (hs[i - 1] + hs[i]));
                else
                    a0[i] = -sigma2 / (2d * h * h);

            if (is_condensed_h)
                a0[N1 - 1] = -sigma2 / (2d * hs[N1 - 1] * hs[N1 - 1]) + alpha_tj / (2d * hs[N1 - 1]);
            else
                a0[N1 - 1] = -sigma2 / (2d * h * h) + alpha_tj / (2d * h);

            if (is_condensed_h)
                b0[0] = sigma2 / (hs[0] * (hs[0] + hs[1])) + sigma2 / (hs[1] * (hs[0] + hs[1])) + r + 1d / in_tau + alpha_tj / hs[1];
            else
                b0[0] = sigma2 / (2d * h * h) + sigma2 / (2d * h * h) + r + 1d / in_tau + alpha_tj / h;
            for (var i = 1; i < N1 - 1; ++i)
                if (is_condensed_h)
                    b0[i] = sigma2 / (hs[i - 1] * (hs[i - 1] + hs[i])) + sigma2 / (hs[i] * (hs[i - 1] + hs[i])) + r + 1d / in_tau + alpha_tj / hs[i];
                else
                    b0[i] = sigma2 / (2d * h * h) + sigma2 / (2d * h * h) + r + 1d / in_tau + alpha_tj / h;

            if (is_condensed_h)
            {
                var val1 = (sigma2 / 2d) * (1d / (hs[N1 - 1] * hs[N1 - 1]) + (a[j] + lambda[j]) / (2d * hs[N1 - 1]));
                var val2 = 0.5d * (1d / in_tau + r) + 1d / (lambda[j] * in_tau * hs[N1 - 1]);
                var val3 = alpha_tj / (2d * hs[N1 - 1]);
                b0[N1 - 1] = val1 + val2 - val3;
            }
            else
            {
                var val1 = (sigma2 / 2d) * (1d / (h * h) + (a[j] + lambda[j]) / (2d * h));
                var val2 = 0.5d * (1d / in_tau + r) + 1d / (lambda[j] * in_tau * h);
                var val3 = alpha_tj / (2d * h);
                b0[N1 - 1] = val1 + val2 - val3;
            }

            if (is_condensed_h)
                c0[0] = -sigma2 / (hs[0] * (hs[0] + hs[1])) - alpha_tj / hs[1];
            else
                c0[0] = -sigma2 / (2d * h * h) - alpha_tj / h;
            for (var i = 1; i < N1 - 2; ++i)
                if (is_condensed_h)
                    c0[i] = -sigma2 / (hs[i] * (hs[i - 1] + hs[i])) - alpha_tj / hs[i];
                else
                    c0[i] = -sigma2 / (2d * h * h) - alpha_tj / h;
            c0[N1 - 2] = 0d;

            if (j == 1)
                Utils.PrintMatrix1(N1, a0, b0, c0, $"A_{j}.txt");
            var u = SolveByTridiagonalMatrixAlgorithm(j, N1, a0, b0, c0, f);
            // Utils.Print(a0, "db");
            // Utils.Print(b0, "dc");
            // Utils.Print(c0, "dd");
            // Utils.Print(u_curr, "u_curr");
            return u;
        }

        private static double gh1(double x, double in_beta) => x <= 0d ? 0d : Math.Pow(x, in_beta);

        private static double gt1(double t, double in_alpha) => t <= 0d ? 0d : Math.Pow(t, in_alpha);

        private static double gh2(double x, double in_beta)
        {
            if (x <= 0d)
                return 0d;

            if (x >= 1d)
            {
                var d = in_beta;
                var e = 1d - in_beta;
                return d * x + e;
            }

            return Math.Pow(x, in_beta);
        }

        private static double gt2(double t, double in_alpha)
        {
            if (t <= 0d)
                return 0d;

            if (t >= 1d)
            {
                var d = in_alpha;
                var e = 1d - in_alpha;
                return d * t + e;
            }

            return Math.Pow(t, in_alpha);
        }

        private static double[] SolveByTridiagonalMatrixAlgorithm(int j, int n, IReadOnlyList<double> a, IReadOnlyList<double> b, IReadOnlyList<double> c, IReadOnlyList<double> d)
        {
            for (var i = 0; i < b.Count; i++)
                if (Math.Abs(b[i]) < Math.Abs(a[i]) + Math.Abs(c[i]))
                    throw new Exception(
                        $"There is no diagonal dominance! j={j} i={i}: {Math.Abs(a[i])} {Math.Abs(b[i])} {Math.Abs(c[i])} sum={Math.Abs(a[i] + c[i])} ");

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

        private static double[] GetS0(IReadOnlyList<double> s0_wave)
        {
            var res = new double[s0_wave.Count];
            for (var i = 0; i < res.Length; i++)
                res[i] = K * Math.Exp(s0_wave[i]);
            return res;
        }

        [SuppressMessage("ReSharper", "UnusedMember.Local")]
        private static void PrintLowerBound(double r, double q, double sigma2, double K)
        {
            var omega = (-r + q + sigma2 / 2d) / sigma2;
            var alpha_min = omega - Math.Sqrt(omega * omega + 2d * r / sigma2);
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

        private static void S0WaveDirectTime(string name, IReadOnlyList<double> arr, double in_tau) // значения S0 до обратного преобразования
        {
            using var writer = new StreamWriter(name!, false);
            writer.WriteLine("TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA');");
            writer.WriteLine("VARIABLES = S0_W t");
            writer.WriteLine("ZONE T='SubZone'");
            var lines = new List<string>();
            var t = 0d;
            var j = 0;
            while (t < T)
            {
                lines.Add($"{arr[j]:e16} {t:e16}");
                t += is_condensed_tau ? taus[j] : in_tau;
                j++;
            }

            var sb = new StringBuilder();
            foreach (var line in lines)
                sb.AppendLine(line);
            var i = Math.Min(arr.Count - 1, lines.Count);
            writer.WriteLine($"I={i} K={1} ZONETYPE=Ordered");
            writer.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE)");
            writer.WriteLine(sb.ToString());
        }

        private static void S0HatDirectTime(string name, IReadOnlyList<double> arr, double in_tau) // значения производных
        {
            using var writer = new StreamWriter(name!, false);
            writer.WriteLine("TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA');");
            writer.WriteLine("VARIABLES = S0_Deriv t");
            writer.WriteLine("ZONE T='SubZone'");
            var lines = new List<string>();
            var t = 0d;
            var j = 0;
            while (t < T)
            {
                lines.Add($"{arr[j]:e16} {t:e16}");
                t += is_condensed_tau ? taus[j] : in_tau;
                j++;
            }

            var sb = new StringBuilder();
            foreach (var line in lines)
                sb.AppendLine(line);
            var I = Math.Min(arr.Count - 1, lines.Count);
            writer.WriteLine($"I={I} K={1} ZONETYPE=Ordered");
            writer.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE)");
            writer.WriteLine(sb.ToString());
        }

        private static void S0ReversedTime(string name, IReadOnlyList<double> arr, double T, double in_tau) // значения S0 после обратного преобразования
        {
            using var writer = new StreamWriter(name!, false);
            writer.WriteLine("TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA');");
            writer.WriteLine("VARIABLES = S0 t");
            writer.WriteLine("ZONE T='SubZone'");
            var t = T;
            var j = 0;
            var s = string.Empty;
            var I = 0;
            if (is_condensed_h)
            {
                while (t >= 0)
                {
                    s += $"{arr[j]:e16} {t:e16}{Environment.NewLine}";
                    t -= taus[j];
                    j++;
                }
            }
            else
            {
                while (t > 0)
                {
                    s += $"{arr[j]:e16} {t:e16}{Environment.NewLine}";
                    t -= in_tau;
                    if (t < 10e-10) // иногда тут остается очень маленькое значение, нам оно не нужно
                        break;
                    j++;
                }

                I = j + 1;
            }

            writer.WriteLine($"I={I} K={1} ZONETYPE=Ordered");
            writer.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE)");
            writer.WriteLine(s);
        }

        private static void WriteVectorToFile(string path, IEnumerable<double> arr)
        {
            var builder = new StringBuilder();
            foreach (var a in arr!)
                builder.AppendLine(a.ToString("e10"));
            File.WriteAllText(path, builder.ToString());
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
            var res = new double[N1];
            for (var i = 1; i < res.Length; i++)
                res[i - 1] = gh2(i * h, beta) - gh2((i - 1) * h, beta);
            return res;
        }

        private static void CheckNaN(double val)
        {
            if (!double.IsNaN(val))
                return;
            Console.Out.Flush();
            throw new Exception("rho is NaN!");
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
    }
}