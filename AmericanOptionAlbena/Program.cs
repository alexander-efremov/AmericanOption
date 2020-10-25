#pragma warning disable 219
namespace AmericanOptionAlbena
{
    using System;
    using System.IO;
    using CoreLib.Utils;

    public static class Program
    {
        public static void Main()
        {
            const double Tol = 10e-3; // eps to refine eta
            const double eps = 10e-3; // eps to refine mu
            const double lb = 0d; // left bound
            const double rb = 50d; // right bound
            const double T0 = 0d; // the start time
            const double Tn = 1d; // the finish time
            const double T = Tn - T0; // time interval
            const int M = 1000; // the number of time intervals
            const double tau = T / M; // the time step
            const double sigma = 0.1d; // the sigma = volatility
            const double sigma2 = sigma * sigma; // the squared sigma
            const double r = 0.1d; // the risk-free rate
            const double K = 10; // the strike price
            const double q = 0.0d; // the dividend rate
            const double b = -((2d * r) / sigma2);
            var lambda = 1d; // the lambda in infinite element approximation
            var s0eps = 10e-3; // eps to refine s0
            var s0 = new double[M + 1]; // the s0 solution
            var N = 800; // the number of space intervals
            var N1 = N + 1; // the number of points
            var h = (rb - lb) / N; // the space step
            var u_curr = new double[N1]; // the current solution
            var u_prev = new double[N1]; // the previous solution

            // for j = 0
            var s0_wave_j = 0d;
            var s0_wave_j_minus_1 = 0d;
            var s0_deriv_hat = 0d;
            var s0_deriv_dash = 0d;
            var s0_dash_j = 0d;
            var s0_dash_j_minus_1 = 0d;
            var eta_j_l_m2 = 0d;
            var eta_j_l_m1 = 0d;
            var eta_j_l = 0d;
            var lambda_j = 0d;
            var lambda_j_minus_1 = 0d;

            for (var j = 1; j <= M; j++)
            {
                var l = 0;
                var mu_j_l_m2 = 0d;
                var mu_j_l_m1 = 0d;
                var mu_j_l = s0_deriv_hat;

                var iterCount = 0;
                while (true)
                {
                    iterCount++;
                    
                    var a_j = (2d * r - 2d * q - sigma2 + 2d * mu_j_l) / sigma2;
                    lambda_j_minus_1 = lambda_j;
                    lambda_j = Math.Sqrt(a_j * a_j - 4d * b);

                    u_curr = Solve(N1, r, q, sigma, s0_deriv_dash, tau, h, u_prev, sigma2, s0_dash_j, mu_j_l, a_j, s0_dash_j_minus_1, lambda_j, lambda_j_minus_1);

                    var rho_j_l = s0[j - 1] + tau * mu_j_l;
                    eta_j_l_m2 = eta_j_l_m1;
                    eta_j_l_m1 = eta_j_l;
                    eta_j_l = K - rho_j_l - u_curr[0];
                    if (Math.Abs(eta_j_l) <= Tol)
                    {
                        s0_deriv_dash = s0_deriv_hat = mu_j_l;
                        s0_dash_j = rho_j_l;
                        s0[j] = s0_dash_j;
                        break;
                    }
                    // ReSharper disable once RedundantIfElseBlock
                    else
                    {
                        l += 1;
                        if (l == 1)
                        {
                            mu_j_l_m1 = mu_j_l;
                            mu_j_l = mu_j_l_m1 + eps;
                        }
                        else
                        {
                            mu_j_l_m2 = mu_j_l_m1;
                            mu_j_l_m1 = mu_j_l;
                            var mu_val1 = eta_j_l_m1 * mu_j_l_m2 - eta_j_l_m2 * mu_j_l_m1;
                            var mu_val2 = eta_j_l_m1 - eta_j_l_m2;
                            mu_j_l = mu_val1 / mu_val2;
                        }
                    }

                    for (var i = 0; i < u_prev.Length; i++)
                    {
                        u_curr[i] = u_prev[i];
                    }

                    s0_dash_j_minus_1 = s0_dash_j;
                }
            }

            // check convexity
            for (var i = 1; i <= s0.Length - 1; i++)
            {
                if (s0[i - 1] > s0[i])
                {
                    throw new InvalidOperationException("The s0 is not convex!");
                }
            }

            // print s0 to tecplot file
            var name = $"{s0}_nx={s0.Length}_tau={tau}.dat";
            using (var writer = new StreamWriter(name, false))
            {
                writer.WriteLine(
                    "TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA'\nVARIABLES = S {0}\nZONE T='{1}'",
                    "t",
                    "SubZone");
                writer.WriteLine($"I={s0.Length} K={1} ZONETYPE=Ordered");
                writer.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE)");
                for (var i = 0; i < s0.Length; i++)
                {
                    writer.WriteLine("{0:e8} {1:e8}", s0[i], i);
                }
            }
        }

        private static double GetAlpha(double r, double q, double sigma, double s0_deriv_dash)
        {
            var val = (sigma * sigma) / 2d;
            return r - q - val + s0_deriv_dash;
        }

        private static double[] Solve(
            int N1,
            double r,
            double q,
            double sigma,
            double s0_deriv_dash,
            double tau,
            double h,
            double[] u_prev,
            double sigma2,
            double s0_dash_tj,
            double mu_j_l,
            double a_j,
            double s0_dash_j_minus_1,
            double lambda_j,
            double lambda_j_minus_1)
        {
            // calculate right part
            var f = new double[N1];
            var alpha = GetAlpha(r, q, sigma, s0_deriv_dash);
            var gamma1 = 1d / tau - alpha / h;
            var gamma2 = alpha / h;

            f[0] = gamma1 * u_prev[0] + gamma2 * u_prev[1] - (sigma2 * s0_dash_tj) / h;
            for (var i = 1; i < N1 - 1; i++)
            {
                f[i] = gamma1 * u_prev[i] + gamma2 * u_prev[i + 1];
            }

            var nu_j = 2d / (lambda_j_minus_1 + lambda_j + 2d * (s0_dash_j_minus_1 - mu_j_l) / sigma2);
            f[N1 - 1] = 0.5d * (gamma1 * u_prev[N1 - 1] + gamma2 * (u_prev[N1 - 1] + h)) + nu_j / (tau * h) * u_prev[N1 - 1];

            Utils.Print(f, "f");
            
            // b - below main diagonal (indexed as [1;n-1])
            var db = new double[N1];
            db[0] = 0d;
            for (var i = 1; i < N1 - 1; ++i)
            {
                db[i] = -sigma2 / (2d * h * h);
            }

            db[N1 - 1] = -sigma2 / (2d * h * h); // right boundary condition

            // main diagonal of matrix (indexed as [0;n-1])
            var dc = new double[N1];
            dc[0] = sigma2 / (h * h) + r + 1d / tau + (-sigma2 / (2d * h * h));
            for (var i = 1; i < N1 - 1; ++i)
            {
                dc[i] = sigma2 / (h * h) + r + 1d / tau;
            }

            var val1 = sigma2 / (2d * h * h);
            var val2 = sigma2 * (a_j + lambda_j) / (4d * h);
            var val3 = 0.5d * (1d / tau + r);
            var val4 = 1d / (lambda_j * tau * h);
            dc[N1 - 1] = val1 + val2 + val3 + val4; // right boundary condition

            // d - up to main diagonal (indexed as [0;n-2])
            var dd = new double[N1];
            dd[0] = -sigma2 / (h * h);
            for (var i = 0; i <= N1 - 2; ++i)
            {
                dd[i] = -sigma2 / (2d * h * h);
            }

            // right boundary condition
            dd[N1 - 2] = 0d;
            // Utils.Print(db, "db");
            // Utils.Print(dc, "dc");
            // Utils.Print(dd, "dd");
            double[] u_curr = SolveByTridiagonalMatrixAlgorithm(N1, db, dc, dd, f);
            Utils.Print(u_curr, "u_curr");
            return u_curr;
        }

        private static double[] SolveByTridiagonalMatrixAlgorithm(int n, double[] b, double[] c, double[] d, double[] r)
        {
            for (var i = 0; i < c.Length; i++)
            {
                if (Math.Abs(c[i]) < Math.Abs(b[i]) + Math.Abs(d[i]))
                {
                    throw new Exception(
                        $"There is no diagonal dominance! i={i} {Math.Abs(b[i])} {Math.Abs(c[i])} {Math.Abs(d[i])} sum={Math.Abs(b[i] + d[i])} ");
                }
            }

            var delta = new double[n];
            var beta = new double[n];
            var lambda = new double[n];

            if (Math.Abs(c[0]) < double.Epsilon)
            {
                throw new InvalidOperationException("c[0] == 0");
            }

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
            {
                x[i] = beta[i] * x[i + 1] + lambda[i];
            }

            return x;
        }
    }
}