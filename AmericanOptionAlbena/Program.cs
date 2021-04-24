#pragma warning disable 219
namespace AmericanOptionAlbena
{
    using System;
    using System.Diagnostics.CodeAnalysis;
    using CoreLib.Utils;

    [SuppressMessage("ReSharper", "ConvertToCompoundAssignment")]
    [SuppressMessage("ReSharper", "SuggestBaseTypeForParameter")]
    public static class Program
    {
        [SuppressMessage("ReSharper", "RedundantAssignment")]
        [SuppressMessage("ReSharper", "TooWideLocalVariableScope")]
        [SuppressMessage("ReSharper", "RedundantIfElseBlock")]
        public static void Main()
        {
            const int l_max_iterations = 1000; // max iteration on l
            const double Tol = 10e-3; // eps to refine eta
            const double eps = 10e-3; // eps to refine mu
            const double lb = 0d; // left bound
            const double rb = 50d; // right bound
            const double T0 = 0d; // the start time
            const double Tn = 1d; // the finish time
            const double T = Tn - T0; // time interval
            const double sigma = 0.1d; // the sigma = volatility
            const double sigma2 = sigma * sigma; // the squared sigma
            const double r = 0.1d; // the risk-free rate
            const double K = 10; // the strike price
            const double q = 0.0d; // the dividend rate
            const double b = -((2d * r) / sigma2);
            const int N = 800; // the number of space intervals
            const double h = (rb - lb) / N; // the space step
            const int N1 = N + 1; // the number of points
            const int M = 1000; // the number of time intervals
            const double tau = T / M; // the time step
            var s0 = new double[M + 1]; // the s0 solution
            var u_curr = new double[N1]; // the current solution
            var u_prev = new double[N1]; // the previous solution
            var s0_hat_deriv = new double[M + 1];
            var s0_dash_deriv = new double[M + 1];
            var s0_dash = new double[M + 1];
            var lambda = new double[M + 1];
            var a = new double[M];

            for (var j = 1; j <= M; j++)
            {
                var mu_j = new double[l_max_iterations];
                var eta_j = new double[l_max_iterations];
                var rho_j = new double[l_max_iterations];
                var l = 0;
                mu_j[l] = s0_hat_deriv[j - 1];
                Console.WriteLine("l = {0} eta_j[l-2] = {1}", l, (l > 1 ? eta_j[l - 2] : 0d));
                Console.WriteLine("l = {0} eta_j[l-1] = {1}", l, (l > 0 ? eta_j[l - 1] : 0d));
                Console.WriteLine("l = {0} eta_j[l] = {1}", l, eta_j[l]);
                Console.WriteLine("l = {0} mu_j[l-2] = {1}", l, 0d);
                Console.WriteLine("l = {0} mu_j[l-1] = {1}", l, 0d);
                Console.WriteLine("l = {0} mu_j[l] = {1}", l, mu_j[l]);
                Console.WriteLine("----------------------------------");
                while (l < l_max_iterations)
                {
                    a[j] = (2d * r - 2d * q - sigma2 + 2d * mu_j[l]) / sigma2;
                    lambda[j] = Math.Sqrt(a[j] * a[j] - 4d * b);
                    u_curr = Solve(l, j, N1, r, q, s0_dash_deriv[j], tau, h, u_prev, sigma2, rho_j, mu_j[l], a[j], lambda);
                    rho_j[l] = s0[j - 1] + tau * mu_j[l];
                    
                    eta_j[l] = K - rho_j[l] - u_curr[0];
                    if (Math.Abs(eta_j[l]) <= Tol)
                    {
                        s0_dash_deriv[j] = s0_hat_deriv[j] = mu_j[l];
                        s0[j] = rho_j[l];
                        break;
                    }

                    if (++l == 1)
                    {
                        mu_j[l] = mu_j[0] + eps;
                        Console.WriteLine("l = {0} eta_j[l-2] = {1}", l, 0d);
                        Console.WriteLine("l = {0} eta_j[l-1] = {1}", l, (l > 0 ? eta_j[l - 1] : 0d));
                        Console.WriteLine("l = {0} eta_j[l] = {1}", l, eta_j[l]);
                        Console.WriteLine("l = {0} mu_j[l-2] = {1}", l, 0d);
                        Console.WriteLine("l = {0} mu_j[l-1] = {1}", l, mu_j[l - 1]);
                        Console.WriteLine("l = {0} mu_j[l] = {1}", l, mu_j[l]);
                        Console.WriteLine("----------------------------------");
                        
                    }
                    else
                    {
                        Console.WriteLine("l = {0} mu_j[l-2] = {1}", l, mu_j[l - 2]);
                        Console.WriteLine("l = {0} mu_j[l-1] = {1}", l, mu_j[l - 1]);
                        Console.WriteLine("l = {0} mu_j[l] = {1}", l, mu_j[l]);
                        Console.WriteLine("l = {0} eta_j[l-2] = {1}", l, eta_j[l - 2]);
                        Console.WriteLine("l = {0} eta_j[l-1] = {1}", l, eta_j[l - 1]);
                        mu_j[l] = (eta_j[l - 1] * mu_j[l - 2] - eta_j[l - 2] * mu_j[l - 1]) / (eta_j[l - 1] - eta_j[l - 2]);
                        Console.WriteLine(
                            "l = {0} (eta_j[l - 1] * mu_j[l - 2] - eta_j[l - 2] * mu_j[l - 1]) = {1:g6}",
                            l,
                            eta_j[l - 1] * mu_j[l - 2] - eta_j[l - 2] * mu_j[l - 1]);
                        Console.WriteLine("l = {0} (eta_j[l - 1] - eta_j[l - 2]) = {1:g6}", l, eta_j[l - 1] - eta_j[l - 2]);
                        Console.WriteLine("l = {0} mu_j[l-2] = {1}", l, mu_j[l - 2]);
                        Console.WriteLine("l = {0} mu_j[l-1] = {1}", l, mu_j[l - 1]);
                        Console.WriteLine("l = {0} calculated mu_j[l] = {1}", l, mu_j[l]);
                        Console.WriteLine("----------------------------------");
                    }
                }

                if (l >= l_max_iterations)
                    throw new InvalidOperationException("The 'l' exceed the maximum number of iterations.");
                
                s0_dash[j] = rho_j[l];
                for (var i = 0; i < u_prev.Length; i++)
                    u_prev[i] = u_curr[i];

                break;
            }

            FileUtils.WriteVectorToFile("S0.dat", s0);

            // // check convexity
            // for (var i = 1; i <= s0.Length - 1; i++)
            // {
            //     if (s0[i - 1] > s0[i])
            //     {
            //         throw new InvalidOperationException("The s0 is not convex!");
            //     }
            // }
            //
            // // print s0 to tecplot file
            // var name = $"{s0}_nx={s0.Length}_tau={tau}.dat";
            // using (var writer = new StreamWriter(name, false))
            // {
            //     writer.WriteLine(
            //         "TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA'\nVARIABLES = S {0}\nZONE T='{1}'",
            //         "t",
            //         "SubZone");
            //     writer.WriteLine($"I={s0.Length} K={1} ZONETYPE=Ordered");
            //     writer.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE)");
            //     for (var i = 0; i < s0.Length; i++)
            //     {
            //         writer.WriteLine("{0:e8} {1:e8}", s0[i], i);
            //     }
            // }
        }

        private static double[] Solve(
            int l,
            int j,
            int N1,
            double r,
            double q,
            double s0_deriv_dash,
            double tau,
            double h,
            double[] u_prev,
            double sigma2,
            double[] rho_j, // to toke l-th approximation of s_dash_derivative (t_j)
            double mu_j_l,
            double a_j,
            double[] lambda)
        {
            // calculate the right part
            var alpha = r - q - ((sigma2) / 2d) + s0_deriv_dash;
            var gamma1 = 1d / tau - alpha / h;
            var gamma2 = alpha / h;

            Console.WriteLine("l = {0} j = {1} s0_dash[j] = {2}", l, j, rho_j[j]);

            var f = new double[N1];
            f[0] = gamma1 * u_prev[0] + gamma2 * u_prev[1] - (sigma2 / h) * rho_j[j];
            for (var i = 1; i < N1 - 1; i++)
                f[i] = gamma1 * u_prev[i] + gamma2 * u_prev[i + 1];
            var nu_j = 2d / (lambda[j - 1] + lambda[j] + 2d * (rho_j[j - 1] - mu_j_l) / sigma2);
            var u_prev_n_plus_1 = 0d; // u_prev[N1]; TODO: How to calculate the value at the index N1
            f[N1 - 1] = 0.5d * (gamma1 * u_prev[N1 - 1] + gamma2 * u_prev_n_plus_1) + (nu_j / (tau * h)) * u_prev[N1 - 1];

            // if (l == 1)
            //     Console.WriteLine();

            double s = 0;
            for (var i = 0; i < f.Length; i++)
                s += f[i];

            // Console.WriteLine(s);

            // Utils.Print(f, "f");

            // b - below main diagonal (indexed as [1;n-1])
            var db = new double[N1];
            db[0] = 0d;
            for (var i = 1; i < N1 - 1; ++i)
                db[i] = -sigma2 / (2d * h * h);

            db[N1 - 1] = -sigma2 / (2d * h * h); // right boundary condition

            // main diagonal of matrix (indexed as [0;n-1])
            var dc = new double[N1];
            dc[0] = sigma2 / (h * h) + r + (1d / tau) + (-sigma2 / (2d * h * h));
            for (var i = 1; i < N1 - 1; ++i)
                dc[i] = sigma2 / (h * h) + r + (1d / tau);

            var val1 = sigma2 / (2d * h * h);
            var val2 = sigma2 * (a_j + lambda[j - 1]) / (4d * h);
            var val3 = 0.5d * (1d / tau + r);
            var val4 = 1d / (lambda[j - 1] * tau * h);
            dc[N1 - 1] = val1 + val2 + val3 + val4; // right boundary condition

            // d - up to main diagonal (indexed as [0;n-2])
            var dd = new double[N1];
            dd[0] = -sigma2 / (h * h);
            for (var i = 1; i < N1 - 2; ++i)
                dd[i] = -sigma2 / (2d * h * h);

            // right boundary condition
            dd[N1 - 2] = 0d;
            // Utils.Print(db, "db");
            // Utils.Print(dc, "dc");
            // Utils.Print(dd, "dd");
            double[] u_curr = SolveByTridiagonalMatrixAlgorithm(N1, db, dc, dd, f);
            // Utils.Print(u_curr, "u_curr");
            return u_curr;
        }

        private static double[] SolveByTridiagonalMatrixAlgorithm(int n, double[] b, double[] c, double[] d, double[] r)
        {
            for (var i = 0; i < c.Length; i++)
                if (Math.Abs(c[i]) < Math.Abs(b[i]) + Math.Abs(d[i]))
                    throw new Exception(
                        $"There is no diagonal dominance! i={i} {Math.Abs(b[i])} {Math.Abs(c[i])} {Math.Abs(d[i])} sum={Math.Abs(b[i] + d[i])} ");
            if (Math.Abs(c[0]) < double.Epsilon)
                throw new InvalidOperationException("c[0] == 0");
            var delta = new double[n];
            var beta = new double[n];
            var lambda = new double[n];
            delta[0] = c[0];
            beta[0] = -d[0] / delta[0];
            lambda[0] = r[0] / delta[0];
            for (var i = 1; i < n - 1; ++i)
            {
                delta[i] = c[i] + b[i] * beta[i - 1];
                beta[i] = -d[i] / delta[i];
                lambda[i] = (r[i] - b[i] * lambda[i - 1]) / delta[i];
            }

            var x = new double[n];
            x[n - 1] = lambda[n - 1];
            for (var i = n - 2; i >= 0; i--)
                x[i] = beta[i] * x[i + 1] + lambda[i];
            return x;
        }
    }
}