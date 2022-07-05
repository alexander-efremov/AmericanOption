using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Text;

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

        private const bool print = true; // the parameter which allows enable/disable printing of the results

        // private const bool is_time_condensed = true; // the parameter which allows enable/disable condensed time mesh
        // private const bool is_space_condensed = true; // the parameter which allows enable/disable condensed space mesh
        private static readonly bool is_condensed_h = false; // the parameter which allows enable/disable condensed meshes
        private static readonly bool is_condensed_tau = true; // the parameter which allows enable/disable condensed meshes
        private static readonly bool is_start_approximation_enabled = true; // the parameter which allows enable/disable aproximate solution 

        private static readonly double[] taus = GetTaus();
        private static readonly double[] hs = GetHs();

        public static void Main()
        {
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
            s0_hat_deriv[0] = -0d;
            s0_wave[0] = 0d;
            a[0] = (2d * r - 2d * q - sigma2 + 2d * s0_hat_deriv[0]) / sigma2; // here s0_hat_deriv[j] == mu_j[l] where l = 0
            lambda[0] = Math.Sqrt(a[0] * a[0] - 4d * b);
            const int j_star = 1; // the time step after that we switch to iteration algorigthm
            for (var j = 1; j <= M; j++)
            {
                var tau = is_condensed_tau ? taus[j - 1] : tau0;
                if (is_start_approximation_enabled && j <= j_star)
                {
                    var theta = tau * j;
                    var mu_j = sigma * (1d - Math.Abs(Math.Log(theta))) / (2d * Math.Sqrt(theta * Math.Abs(Math.Log(theta))) * (1d - sigma * Math.Sqrt(theta * Math.Abs(Math.Log(theta)))));
                    var rho_j_l = Math.Log(1d - sigma * Math.Sqrt(theta * Math.Abs(Math.Log(theta))));
                    CheckNaN(rho_j_l);
                    a[j] = (2d * r - 2d * q - sigma2 + 2d * mu_j) / sigma2;
                    lambda[j] = Math.Sqrt(a[j] * a[j] - 4d * b);
                    u_curr = Solve(j, u_prev, rho_j_l, s0_hat_deriv, mu_j, a, lambda, tau);
                    if (print)
                        Console.WriteLine("j* = {0} mu_j_l {5} K {1} rho_j_l {2} u[1] {3} u[0] {4}", j, K, rho_j_l, u_curr[1], u_curr[0], mu_j);
                    s0_hat_deriv[j] = mu_j;
                    s0_wave[j] = rho_j_l;
                }
                else
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
                        CheckNaN(rho_j_l);
                        u_curr = Solve(j, u_prev, rho_j_l, s0_hat_deriv, mu_j[l], a, lambda, tau);
                        var h0 = is_condensed_h ? hs[0] : h;
                        eta_j[l] = K * Math.Exp(rho_j_l) + (u_curr[1] - u_curr[0]) / h0;
                        if (print)
                            Console.WriteLine("l = {0} j = {1} mu_j_l {7} K {2} rho_j_l {3} u[1] {4} u[0] {5} eta_j_l = {6}", l, j, K, rho_j_l, u_curr[1], u_curr[0], eta_j[l], mu_j[l]);
                        if (Math.Abs(eta_j[l]) <= Tol)
                        {
                            s0_hat_deriv[j] = mu_j[l];
                            s0_wave[j] = rho_j_l;
                            break;
                        }

                        if (++l == 1)
                            mu_j[l] = mu_j[l - 1] - eps;
                        else
                            mu_j[l] = (eta_j[l - 1] * mu_j[l - 2] - eta_j[l - 2] * mu_j[l - 1]) / (eta_j[l - 1] - eta_j[l - 2]);
                    }
                }

                for (var i = 0; i < u_prev.Length; i++)
                    u_prev[i] = u_curr[i];
            }

            WriteVectorToFile("S0.arr", s0_wave);
            WriteVectorToFile("S0_hat_deriv.arr", s0_hat_deriv);
            // значения производных
            S0HatDirectTime(
                $"{(is_condensed_h ? "1" : "11")}_{nameof(s0_hat_deriv)}_N1={N1}_T={T}_h_condensed_{is_condensed_h}_tau_condensed_{is_condensed_tau}_approx_{is_start_approximation_enabled}.dat",
                s0_hat_deriv, tau0);
            // значения S0 до обратного преобразования
            S0WaveDirectTime($"{(is_condensed_h ? "2" : "22")}_{nameof(s0_wave)}_N1={N1}_T={T}_h_condensed_{is_condensed_h}_tau_condensed_{is_condensed_tau}_approx_{is_start_approximation_enabled}.dat",
                s0_wave, tau0);
            var s0 = GetS0(s0_wave);
            // значения S0 после обратного преобразования
            S0ReversedTime($"{(is_condensed_h ? "3" : "33")}_{nameof(s0)}_T-t_N1={N1}_T={T}_h_condensed_{is_condensed_h}_tau_condensed_{is_condensed_tau}_approx_{is_start_approximation_enabled}.dat", s0,
                T, tau0);
        }

        private static double[] Solve(int j, double[] u_prev, double rho_j_l, double[] s0_hat_deriv, double mu_j_l, double[] arr_a, double[] lambda, double in_tau)
        {
            // calculate the right part
            // var alpha = r - q - (sigma2 / 2d) + s0_wave[j - 1];
            var f = new double[N1];
            var alpha_jm1 = r - q - sigma2 / 2d + s0_hat_deriv[j - 1];
            f[0] = K * (1d - Math.Exp(rho_j_l));
            for (var i = 1; i < N1 - 1; i++)
            {
                var h_jp12 = is_condensed_h ? hs[i] : h;
                var gamma1_jm1_0 = 1d / in_tau - alpha_jm1 / h_jp12;
                var gamma2_jm1_0 = alpha_jm1 / h_jp12;
                f[i] = gamma1_jm1_0 * u_prev[i] + gamma2_jm1_0 * u_prev[i + 1];
            }

            // formula 31 is applied to j-1 time level
            // here /*x_{n+1}-x_n= x_n + h - x_n = h*/
            {
                var nu_j = 2d / (lambda[j - 1] + lambda[j] + 2d * (s0_hat_deriv[j - 1] - mu_j_l) / sigma2);
                var h_jp12 = is_condensed_h ? hs[N1 - 1] : h;
                var gamma1_jm1 = 1d / in_tau - alpha_jm1 / h_jp12;
                var gamma2_jm1 = alpha_jm1 / h_jp12;
                var u_prev_np1 = u_prev[N1 - 1] * Math.Exp(-((lambda[j - 1] + arr_a[j - 1]) * h_jp12) / 2d);
                f[N1 - 1] = 0.5d * (gamma1_jm1 * u_prev[N1 - 1] + gamma2_jm1 * u_prev_np1) + nu_j * u_prev[N1 - 1] / (in_tau * h_jp12);
            }

            var a0 = new double[N1];
            var b0 = new double[N1];
            var c0 = new double[N1];

            // a - below main diagonal (indexed as [1;n-1])
            a0[0] = 0d;
            for (var i = 1; i < N1 - 1; ++i)
                if (is_condensed_h)
                    a0[i] = -sigma2 / (hs[i - 1] * (hs[i - 1] * hs[i]));
                else
                    a0[i] = -sigma2 / (2d * h * h);

            if (is_condensed_h)
                a0[N1 - 1] = -sigma2 / (hs[N1 - 2] * (hs[N1 - 2] + hs[N1 - 1]));
            else
                a0[N1 - 1] = -sigma2 / (2d * h * h);

            // main diagonal of matrix (indexed as [0;n-1])
            // b0[0] = sigma2 / (h * h) + r + (1d / tau);
            b0[0] = 1d;
            for (var i = 1; i < N1 - 1; ++i)
                if (is_condensed_h)
                    b0[i] = sigma2 / (hs[i - 1] * hs[i]) + r + 1d / in_tau;
                else
                    b0[i] = sigma2 / (h * h) + r + 1d / in_tau;

            if (is_condensed_h)
            {
                // what h we use here? hp12 or hm12?
                var val1 = sigma2 / (hs[N1 - 2] * hs[N1 - 1]);
                var val2 = sigma2 * (arr_a[j] + lambda[j]) / (4d * hs[N1 - 1]);
                var val3 = 0.5d * (1d / in_tau + r);
                var val4 = 1d / (lambda[j] * in_tau * hs[N1 - 1]);
                b0[N1 - 1] = val1 + val2 + val3 + val4;
            }
            else
            {
                var val1 = sigma2 / (2d * h * h);
                var val2 = sigma2 * (arr_a[j] + lambda[j]) / (4d * h);
                var val3 = 0.5d * (1d / in_tau + r);
                var val4 = 1d / (lambda[j] * in_tau * h);
                b0[N1 - 1] = val1 + val2 + val3 + val4;
            }

            // c - up to main diagonal (indexed as [0;n-2])
            // c0[0] = -sigma2 / (h * h);
            c0[0] = 0d;
            for (var i = 1; i < N1 - 2; ++i)
                if (is_condensed_h)
                    c0[i] = -sigma2 / (hs[i] * (hs[i - 1] + hs[i]));
                else
                    c0[i] = -sigma2 / (2d * h * h);

            c0[N1 - 2] = 0d;

            // Utils.Print(db, "db");
            // Utils.Print(dc, "dc");
            // Utils.Print(dd, "dd");
            var u = SolveByTridiagonalMatrixAlgorithm(N1, a0, b0, c0, f);
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

        private static double[] SolveByTridiagonalMatrixAlgorithm(int n, IReadOnlyList<double> a,
            IReadOnlyList<double> b, IReadOnlyList<double> c, IReadOnlyList<double> d)
        {
            for (var i = 0; i < b.Count; i++)
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

        private static void S0WaveDirectTime(string name, IReadOnlyList<double> arr, double in_tau)
        {
            using var writer = new StreamWriter(name!, false);
            writer.WriteLine("TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA');");
            writer.WriteLine("VARIABLES = S0_W t");
            writer.WriteLine("ZONE T='SubZone'");
            var lines = new List<string>();
            var t = 0d;
            var j = 0;
            while (t <= T)
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

        private static void S0HatDirectTime(string name, IReadOnlyList<double> arr, double in_tau)
        {
            using var writer = new StreamWriter(name!, false);
            writer.WriteLine("TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA');");
            writer.WriteLine("VARIABLES = S0_Deriv t");
            writer.WriteLine("ZONE T='SubZone'");
            var lines = new List<string>();
            var t = 0d;
            var j = 0;
            while (t <= T)
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

        private static void S0ReversedTime(string name, IReadOnlyList<double> arr, double T, double in_tau)
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
            if (double.IsNaN(val))
                throw new Exception("rho is NaN!");
        }
    }
}