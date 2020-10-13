#pragma warning disable 219
namespace AmericanOptionAlbena
{
    using System;
    using System.Diagnostics.CodeAnalysis;
    using System.IO;

    public static class Program
    {
        [SuppressMessage("Style", "IDE0059:Unnecessary assignment of a value", Justification = "<Pending>")]
        [SuppressMessage("ReSharper", "UnusedVariable")]
        public static void Main()
        {
            double a = 0d; // left bound
            double b = 50d; // right bound
            int N = 10; // the number of space intervals
            int N1 = N + 1; // the number of points
            double h = (b - a) / N; // the space step
            double T0 = 0d; // the start time
            double Tn = 1d; // the finish time
            double T = (Tn - T0); // time interval
            int M = 1000; // the number of time intervals
            int M1 = M + 1; // the number of time points
            double tau = T / M; // the time step
            double sigma = 0.1d; // the sigma = volatility
            double sigma2 = sigma * sigma; // the squared sigma
            double r = 0.1d; // the risk-free rate
            double K = 10; // the strike price
            double K2 = K * K; // the squared strike price
            double q = 0.1d; // the dividend rate
            double lambda = 1d; // the lambda in infinite element approximation
            double s0eps = 10e-3; // eps to refine s0
            double[] s0 = new double[M1]; // the s0 solution
            for (int i = 0; i < s0.Length; i++)
            {
                s0[i] = 0d;
            }
            double[] u1 = new double[N1]; // the current solution
            for (int i = 0; i < u1.Length; i++)
            {
                u1[i] = 0d;
            }

            double[] u0 = new double[N1]; // the previous solution

            // lets start the calculations
            // initial condition of u(x,0) <=> 0 
            for (int i = 0; i < u0.Length; i++)
            {
                u0[i] = 0d;
            }

            // initial condition for s0(0) is 0
            // as far as after 2 changes of variables 
            // condition (S(T) = K, T) in (S,t)-coordinates
            // was transformed into (0,0) in (x,t) coordinates
            s0[0] = 0d;

            // for each time step
            for (int k = 1; k <= M; k++)
            {
                // save previous s0 from k-1
                double s0p = s0[k - 1];

                double s0c;
                do
                {
                    // calculate uc
                    // calculate right part
                    double[] f = new double[N1];
                    double s_t_shtrih = 1d; //TODO!
                    f[0] = 0d; // u(0,t)=K-s_0(t)
                    for (int i = 1; i < N1; i++)
                    {
                        var ga = -(r - q - (sigma2 / 2d) + (s_t_shtrih / K)) * K;
                        double x = i * h;
                        double xl = x - h / 2d;
                        double xr = x + h / 2d;
                        double xlp = xl - ga*tau;
                        double xrp = xr - ga*tau;

                        // will use gammas from 
                        // Vladimir V. Shaidurov*, Alexander V. Vyatkin, and Elena V. Kuchunova
                        // Semi-Lagrangian difference approximations with different stability requirements 2018
                        // (2.25) - gammas for L_inf for du/dt+a(du/dx)
                        double gammmal = (1d / (8d * tau)) * (1d - ((4d * tau) / h) * ga);
                        double gammmac = (1d / (8d * tau)) * (6d + ((4d * tau) / h) * ga - ((4d * tau) / h) * ga);
                        double gammmar = (1d / (8d * tau)) * (1d + ((4d * tau) / h) * ga);

                        double uvl = u0[i - 1];
                        double uvc = u0[i];
                        double uvr = u0[i + 1];
                        f[i] = gammmal * uvl + gammmac * uvc + gammmar * uvr;
                    }

                    // thomas
                    // b - диагональ, лежащая под главной (нумеруется: [1;n-1])
                    var db = new double[N1];
                    for (var i = 1; i <= N1 - 1; ++i)
                    {
                        var x = i * h;
                        db[i] = 0d;
                    }
                    // right boundary cond
                    db[N1 - 1] = 0d;
                    // главная диагональ матрицы A (нумеруется: [0;n-1])
                    double[] dc = new double[N1];
                    for (var i = 1; i <= N1 - 1; ++i)
                    {
                        var x = i * h;
                        dc[i] = 0d;
                    }
                    // right boundary condition
                    dc[N1 - 1] = 0d;
                    // диагональ, лежащая над главной (нумеруется: [0;n-2])
                    double[] dd = new double[N1];
                    for (var i = 0; i <= N1 - 2; ++i)
                    {
                        var x = i * h;
                        dd[i] = 0d;
                    }
                    // right boundary condition
                    dd[N1 - 2] = 0d;

                    u1 = ThomasAlgo(N1, db, dc, dd, f);

                    // calculate s0c = s0_k
                    // from (20)
                    s0c = K - u1[0];
                } while (Math.Abs(s0c - s0p) > s0eps);

                // update current values before get on next time step
                s0[k] = s0c;
                // copy current to prev and zero current
                for (int tmp = 0; tmp < u1.Length; tmp++)
                {
                    u0[tmp] = u1[tmp];
                    u1[tmp] = 0d;
                }
            }

            // check convexity
            for (var i = 1; i < s0.Length; i--)
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

        private static double[] ThomasAlgo(int n, double[] b, double[] c, double[] d, double[] r)
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
