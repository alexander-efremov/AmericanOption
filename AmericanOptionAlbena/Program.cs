#pragma warning disable 219
#pragma warning disable CS0162
namespace AmericanOptionAlbena
{
    using System;
    using System.Collections.Generic;
    using System.Diagnostics.CodeAnalysis;
    using System.IO;
    using System.Text;

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
        private const int N0 = (int)(N / 100d); // the number of points before h0

        private const double h = (rb - lb) / N; // the space step

        // private const int M = 10000; // the number of time intervals
        // private const double M0 = 166; // the time step until we condense the time mesh
        private const int M = 10000; // the number of time intervals
        private const double M0 = 166; // the time step until we condense the time mesh
        private const double tau = T / M; // the time step
        private const double alpha = 2d; // condense parameter for t
        private const double beta = 2d; // condense parameter for h

        private const bool print = true; // the parameter which allows enable/disable printing of the results

        // private const bool is_time_condensed = true; // the parameter which allows enable/disable condensed time mesh
        // private const bool is_space_condensed = true; // the parameter which allows enable/disable condensed space mesh
        private const bool is_condensed_meshes = true; // the parameter which allows enable/disable condensed meshes

        public static void Main()
        {
            // PrintCondensedMeshesTau();
            // PrintCondensedMeshesH();
            // return;
            // PrintCondensedMeshes1();
            // PrintCondensedMeshes2();
            // return;
            CheckParameters();
            var u_curr = new double[N1]; // the current solution
            var u_prev = new double[N1]; // the previous solution
            var s0_hat_deriv = new double[M + 1];
            var s0_wave = new double[M + 1];
            var lambda = new double[M + 1];
            var a = new double[M + 1];

            // s0_hat_deriv[0] = s0_wave[0] = -1; // failed on condition var alpha = r - q - (sigma2 / 2d) + mu_j_l; Debug.Assert(alpha >= 0, "alpha>=0");
            // s0_hat_deriv[j] = s0_wave[j] = -sigma2 / 2d;
            const double b = -2d * r / sigma2;
            s0_hat_deriv[0] = -16;
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
                    var tau0 = tau;
                    if (is_condensed_meshes && j < M0)
                        tau0 = Get_tau2(j, tau, alpha);
                    var rho_j_l = s0_wave[j - 1] + tau0 * mu_j[l];
                    if (double.IsNaN(rho_j_l))
                    {
                        Console.WriteLine("rho is NaN!");
                        break;
                    }

                    u_curr = Solve(j, u_prev, rho_j_l, s0_hat_deriv, mu_j[l], a, lambda);
                    // eta_j[l] = K * (1d - Math.Exp(rho_j_l)) - u_curr[0];

                    var h0 = h;
                    if (is_condensed_meshes && 1 < N0)
                        h0 = Get_h2(1, h, beta);
                    eta_j[l] = K * Math.Exp(rho_j_l) + (u_curr[1] - u_curr[0]) / h0;

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
            PrintS0HatDirectTime($"{(is_condensed_meshes ? "1" : "11")}_{nameof(s0_hat_deriv)}_N1={N1}_T={T}_condensed_{is_condensed_meshes}.dat", s0_hat_deriv, tau, alpha, M0, is_condensed_meshes);
            // значения S0 до обратного преобразования
            PrintS0WaveDirectTime($"{(is_condensed_meshes ? "2" : "22")}_{nameof(s0_wave)}_N1={N1}_T={T}_condensed_{is_condensed_meshes}.dat", s0_wave, tau, alpha, M0, is_condensed_meshes);
            var s0 = GetS0(s0_wave);
            // значения S0 после обратного преобразования
            PrintS0ReversedTime($"{(is_condensed_meshes ? "3" : "33")}_{nameof(s0)}_T-t_N1={N1}_T={T}_condensed_{is_condensed_meshes}.dat", s0, T, tau, alpha, M0, is_condensed_meshes);
        }

        private static double[] Solve(int j, double[] u_prev, double rho_j_l, double[] s0_hat_deriv, double mu_j_l,
            double[] arr_a, double[] lambda)
        {
            // calculate the right part
            // var alpha = r - q - (sigma2 / 2d) + s0_wave[j - 1];
            var f = new double[N1];
            var alpha_jm1 = r - q - sigma2 / 2d + s0_hat_deriv[j - 1];
            f[0] = K * (1d - Math.Exp(rho_j_l));
            for (var i = 1; i < N1 - 1; i++)
            {
                var h0 = h;
                var tau0 = tau;
                if (is_condensed_meshes && i < N0)
                    h0 = Get_h2(i, h, beta);
                if (is_condensed_meshes && j < M0)
                    tau0 = Get_tau2(j, tau, alpha);
                var gamma1_jm1_0 = 1d / tau0 - alpha_jm1 / h0;
                var gamma2_jm1_0 = alpha_jm1 / h0;
                f[i] = gamma1_jm1_0 * u_prev[i] + gamma2_jm1_0 * u_prev[i + 1];
            }

            var nu_j = 2d / (lambda[j - 1] + lambda[j] + 2d * (s0_hat_deriv[j - 1] - mu_j_l) / sigma2);

            // formula 31 is applied to j-1 time level
            // here /*x_{n+1}-x_n= x_n + h - x_n = h*/
            var gamma1_jm1 = 1d / tau - alpha_jm1 / h;
            var gamma2_jm1 = alpha_jm1 / h;
            var u_prev_np1 = u_prev[N1 - 1] * Math.Exp(-((lambda[j - 1] + arr_a[j - 1]) * h) / 2d);
            f[N1 - 1] = 0.5d * (gamma1_jm1 * u_prev[N1 - 1] + gamma2_jm1 * u_prev_np1) +
                        nu_j * u_prev[N1 - 1] / (tau * h);

            var a0 = new double[N1];
            var b0 = new double[N1];
            var c0 = new double[N1];

            // a - below main diagonal (indexed as [1;n-1])
            a0[0] = 0d;
            for (var i = 1; i < N1 - 1; ++i)
            {
                if (is_condensed_meshes && i < N0)
                    a0[i] = -sigma2 / (Get_h2(i, h, beta) + Get_h2(i + 1, h, beta));
                else
                    a0[i] = -sigma2 / (2d * h * h);
            }

            a0[N1 - 1] = -sigma2 / (2d * h * h);

            // main diagonal of matrix (indexed as [0;n-1])
            // b0[0] = sigma2 / (h * h) + r + (1d / tau);
            b0[0] = 1d;
            for (var i = 1; i < N1 - 1; ++i)
            {
                if (is_condensed_meshes && i < N0)
                {
                    var tau0 = j < M0 ? Get_tau2(j, tau, alpha) : tau;
                    b0[i] = sigma2 / (Get_h2(i, h, beta) * Get_h2(i + 1, h, beta)) + r +
                            1d / tau0;
                }
                else
                {
                    b0[i] = sigma2 / (h * h) + r + 1d / tau;
                }
            }

            var val1 = sigma2 / (2d * h * h);
            var val2 = sigma2 * (arr_a[j] + lambda[j]) / (4d * h);
            var val3 = 0.5d * (1d / tau + r);
            var val4 = 1d / (lambda[j] * tau * h);
            b0[N1 - 1] = val1 + val2 + val3 + val4;

            // c - up to main diagonal (indexed as [0;n-2])
            // c0[0] = -sigma2 / (h * h);
            c0[0] = 0d;
            for (var i = 1; i < N1 - 2; ++i)
            {
                if (is_condensed_meshes && i < N0)
                    c0[i] = -sigma2 / (Get_h2(i, h, beta) + Get_h2(i + 1, h, beta));
                else
                    c0[i] = -sigma2 / (2d * h * h);
            }

            c0[N1 - 2] = 0d;

            // Utils.Print(db, "db");
            // Utils.Print(dc, "dc");
            // Utils.Print(dd, "dd");
            var u = SolveByTridiagonalMatrixAlgorithm(N1, a0, b0, c0, f);
            // Utils.Print(u_curr, "u_curr");
            return u;
        }

        private static double Get_h1(int in_i, double in_h, double in_beta)
        {
            if (in_i > N0)
                throw new InvalidOperationException("i > N0!");
            return Math.Pow(in_i * in_h, in_beta);
        }

        private static double Get_tau1(int in_j, double in_tau, double in_alpha)
        {
            if (in_j > M0)
                throw new InvalidOperationException("j > M0!");
            return Math.Pow(in_j * in_tau, in_alpha);
        }

        private static double Get_h2(int in_i, double in_h, double in_beta)
        {
            if (in_i > N0)
                throw new InvalidOperationException("i > N0!");
            // начинаем уменьшать шаг после некоторой точки
            var x = in_i * in_h;
            if (x > 1d)
            {
                var d = in_beta;
                var e = 1d - in_beta;
                return d * x + e;
            }

            return Math.Pow(x, in_beta);
        }

        private static double Get_tau2(int in_j, double in_tau, double in_alpha)
        {
            if (in_j > M0)
                throw new InvalidOperationException("j > M0!");

            var t = in_j * in_tau;
            // начинаем уменьшать шаг после некоторой точки M00
            if (t > 1d)
            {
                var d = in_alpha;
                var e = 1d - in_alpha;
                return d * t + e;
            }

            return Math.Pow(t, in_alpha);
        }

        private static void CheckParameters()
        {
            if (N0 >= N1)
                throw new ArgumentException("N0 >= N1!");
            if (M0 >= M)
                throw new ArgumentException("M0 >= M!");
        }

        [SuppressMessage("ReSharper", "UnusedMember.Local")]
        private static void PrintCondensedMeshesTau()
        {
            var t = false;
            Console.WriteLine("V1" + new string(' ', 9) + "V2");
            for (var j = 0; j < M0 + 10; j++)
            {
                if (j < M0)
                {
                    var v1 = Get_tau1(j, tau, alpha);
                    var v2 = Get_tau2(j, tau, alpha);
                    Console.Write(v1.ToString("F8") + " " + v2.ToString("F8") + (v1.Equals(v2) ? '=' : '\0'));
                    Console.WriteLine();
                }
                else
                {
                    if (!t)
                    {
                        Console.WriteLine("constant");
                        t = true;
                    }

                    Console.Write((j * tau).ToString("F8") + " " + (j * tau).ToString("F8"));
                    Console.WriteLine();
                }
            }
        }

        [SuppressMessage("ReSharper", "UnusedMember.Local")]
        private static void PrintCondensedMeshesH()
        {
            var t = false;
            Console.WriteLine("V1" + new string(' ', 9) + "V2");
            for (var i = 0; i < N0 + 10; i++)
            {
                if (i < N0)
                {
                    var v1 = Get_h1(i, h, beta);
                    var v2 = Get_h2(i, h, beta);
                    Console.Write(v1.ToString("F8") + " " + v2.ToString("F8") + (v1.Equals(v2) ? '=' : '\0'));
                    Console.WriteLine();
                }
                else
                {
                    if (!t)
                    {
                        Console.WriteLine("constant");
                        t = true;
                    }

                    Console.Write((i * h).ToString("F8") + " " + (i * h).ToString("F8"));
                    Console.WriteLine();
                }
            }
        }

        // [SuppressMessage("ReSharper", "UnusedMember.Local")]
        // private static void PrintCondensedMeshes1()
        // {
        //     // for (var i = 0; i < N0 + 10; i++)
        //     //     if (i < N0)
        //     //         Console.Write(Get_h1(i, h, beta) + " ");
        //     //     else
        //     //         Console.Write(i * h + " ");
        //
        //     // Console.WriteLine("\r\n=========================\r\n");
        //     for (var j = 0; j < M0 + 10; j++)
        //     {
        //         if (j < M0)
        //             Console.Write(Get_tau1(j, tau, alpha).ToString("F8") + " ");
        //         else
        //             Console.Write((j * tau).ToString("F8") + " ");
        //     }
        //
        //     Console.WriteLine("\r\n=========================\r\n");
        // }

        // [SuppressMessage("ReSharper", "UnusedMember.Local")]
        // private static void PrintCondensedMeshes2()
        // {
        //     // for (var i = 0; i < N0 + 10; i++)
        //     //     if (i < N0)
        //     //         Console.Write(Get_h1(i, h, beta) + " ");
        //     //     else
        //     //         Console.Write(i * h + " ");
        //
        //     // Console.WriteLine("\r\n=========================\r\n");
        //     for (var j = 0; j < M0 + 10; j++)
        //     {
        //         if (j < M0)
        //             Console.Write(Get_tau2(j, tau, alpha).ToString("G8") + " ");
        //         else
        //             Console.Write((j * tau).ToString("G8") + " ");
        //     }
        // }

        private static double[] SolveByTridiagonalMatrixAlgorithm(int n, IReadOnlyList<double> a,
            IReadOnlyList<double> b, IReadOnlyList<double> c, IReadOnlyList<double> d)
        {
            // https://pro-prof.com/forums/topic/sweep-method-for-solving-systems-of-linear-algebraic-equations
            for (var i = 0; i < b.Count; i++)
            {
                if (Math.Abs(b[i]) < Math.Abs(a[i]) + Math.Abs(c[i]))
                    throw new Exception(
                        $"There is no diagonal dominance! i={i} {Math.Abs(a[i])} {Math.Abs(b[i])} {Math.Abs(c[i])} sum={Math.Abs(a[i] + c[i])} ");
            }
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

        private static void PrintS0ReversedTime(string name, IReadOnlyList<double> arr, double T, double in_tau, double in_alpha, double in_M0, bool in_is_condensed_meshes)
        {
            using var writer = new StreamWriter(name!, false);
            writer.WriteLine("TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA');");
            writer.WriteLine("VARIABLES = S0 t");
            writer.WriteLine("ZONE T='SubZone'");
            writer.WriteLine($"I={arr.Count - 1} K={1} ZONETYPE=Ordered");
            writer.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE)");
            for (var i = 0; i < arr.Count; i++)
            {
                var tau0 = in_tau;
                if (in_is_condensed_meshes && i < in_M0)
                    tau0 = Get_tau2(i, in_tau, in_alpha);
                writer.WriteLine("{0:e16} {1:e16}", arr[i], T - tau0 * i);
            }
        }

        private static void PrintS0WaveDirectTime(string name, IReadOnlyList<double> arr, double in_tau, double in_alpha, double in_M0, bool in_is_condensed_meshes)
        {
            using var writer = new StreamWriter(name!, false);
            writer.WriteLine("TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA');");
            writer.WriteLine("VARIABLES = S0_W t");
            writer.WriteLine("ZONE T='SubZone'");
            writer.WriteLine($"I={arr.Count - 1} K={1} ZONETYPE=Ordered");
            writer.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE)");
            for (var i = 0; i < arr.Count; i++)
            {
                var tau0 = in_tau;
                if (in_is_condensed_meshes && i < in_M0)
                    tau0 = Get_tau2(i, in_tau, in_alpha);
                writer.WriteLine("{0:e16} {1:e16}", arr[i], tau0 * i);
            }
        }

        private static void PrintS0HatDirectTime(string name, IReadOnlyList<double> arr, double in_tau, double in_alpha, double in_M0, bool in_is_condensed_meshes)
        {
            using var writer = new StreamWriter(name!, false);
            writer.WriteLine("TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA');");
            writer.WriteLine("VARIABLES = S0_Deriv t");
            writer.WriteLine("ZONE T='SubZone'");
            writer.WriteLine($"I={arr.Count - 1} K={1} ZONETYPE=Ordered");
            writer.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE)");
            for (var i = 0; i < arr.Count; i++)
            {
                var tau0 = in_tau;
                if (in_is_condensed_meshes && i < in_M0)
                    tau0 = Get_tau2(i, in_tau, in_alpha);
                writer.WriteLine("{0:e16} {1:e16}", arr[i], tau0 * i);
            }
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