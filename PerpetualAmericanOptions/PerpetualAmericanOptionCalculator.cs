using System;
using CoreLib;

// ReSharper disable CommentTypo

namespace PerpetualAmericanOptions
{
    // from new presentation with FEM
    public class PerpetualAmericanOptionCalculator : AmericanOptionCalculator
    {
        public PerpetualAmericanOptionCalculator(PerpetualParameters parameters) : base(parameters)
        {
        }
        
        private double[] CalculateRightPart(double S0)
        {
            var V = new double[GetN1() + 1]; // add + 1 to calculate right boundary
            for (var i = 0; i <= GetN1(); i++)
            {
                var si = S0 + i * GetH();
                V[i] = GetV(GetSquaredSigma(), GetR(), GetK(), si);
            }

            var rp = new double[GetN1()];
            var si0 = S0 + 0 * GetH();
            var hph0 = S0 + (0d + 1d) * GetH() - si0; // h_{i+1/2}
            var hmh0 = si0 - (S0 + (0d - 1d) * GetH()); // h_{i-1/2}
            var sph0 = si0 + 0.5d * hph0; // s_{i+1/2}
            var smh0 = si0 - 0.5d * hmh0; // s_{i-1/2}
            CheckHCorrectness(hph0, GetTau(), sph0, GetR());
            CheckHCorrectness(hmh0, GetTau(), smh0, GetR());
            var beta20 = 1d / (8d * GetTau()) *
                         (3d - 2d * GetTau() * GetR() * sph0 / hph0) *
                         (1d + 2d * GetTau() * GetR() * sph0 / hph0);
            var beta30 = 1d / (8d * GetTau()) *
                         (1d - 2d * GetTau() * GetR() * sph0 / hph0) *
                         (1d - 2d * GetTau() * GetR() * sph0 / hph0);
            rp[0] = hph0 * (beta20 * V[0] + beta30 * V[1]) + (GetSquaredSigma() * si0 * si0) / 2d;

            for (var i = 1; i <= GetN1() - 1; i++)
            {
                var si = S0 + i * GetH();
                var hmh = si - (S0 + (i - 1) * GetH()); // h_{i-1/2}
                var hph = S0 + (i + 1) * GetH() - si; // h_{i+1/2}
                var smh = si - 0.5d * hmh; // s_{i-1/2}
                var sph = si + 0.5d * hph; // s_{i+1/2}
                CheckHCorrectness(hph, GetTau(), sph, GetR());
                CheckHCorrectness(hmh, GetTau(), sph, GetR());
                var beta1 = 1d / (8d * GetTau()) *
                            (1d + 2d * GetTau() * GetR() * smh / GetH()) *
                            (1d + 2d * GetTau() * GetR() * smh / GetH());
                var beta2 = 1d / (8d * GetTau()) *
                            (3d + 2d * GetTau() * GetR() * smh / GetH()) *
                            (1d - 2d * GetTau() * GetR() * smh / GetH())
                            +
                            1d / (8d * GetTau()) *
                            (3d - 2d * GetTau() * GetR() * sph / GetH()) *
                            (1d + 2d * GetTau() * GetR() * sph / GetH());
                var beta3 = 1d / (8d * GetTau()) *
                            (1d - 2d * GetTau() * GetR() * sph / GetH()) *
                            (1d - 2d * GetTau() * GetR() * sph / GetH());
                var val = ((hph + hmh) / 2d) * (beta1 * V[i - 1] + beta2 * V[i] + beta3 * V[i + 1]);

                rp[i] = GetF(i) + val;
            }

            return rp;
        }

        // it is zero for our case
        private double GetF(int i)
        {
            return 0d;
        }

        private static void CheckHCorrectness(double h, double tau, double sph, double r)
        {
            if (tau / h > 1d / (2d * r * sph))
            {
                throw new ArgumentOutOfRangeException("tau/h");
            }
        }

        public Tuple<double[], double> Solve()
        {
            var tecplotPrinter = new TecplotPrinterSpecial(GetN1(),
                0d,
                GetRightBoundary(),
                GetTau()); 
            var printer = new ThomasArrayPrinter();

            double[] V = null;
            var S0 = GetK();
            var iter = 0;
            while (Math.Abs(GetExactS0() - S0) > GetS0Eps())
            {
                iter++;
                UpdateH(S0);
                Console.WriteLine("S0 = " + S0);
                Console.WriteLine("h = " + GetH());
                double[] b_t = GetB(GetN1(), S0, GetH(), GetSquaredSigma(), GetTau());
                double[] c_t = GetC(GetN1(), S0, GetH(), GetSquaredSigma(), GetTau(), GetR());
                double[] d_t = GetD(GetN1(), S0, GetH(), GetSquaredSigma(), GetTau());
                double[] rp = CalculateRightPart(S0);
                V = ThomasAlgorithmCalculator.Calculate(b_t, c_t, d_t, rp);
                
                //printer.PrintThomasArrays(b_t, c_t, d_t);
                tecplotPrinter.PrintXY(Environment.GetFolderPath(Environment.SpecialFolder.Desktop) + "\\" + "perpetual-rp", 0d, GetH(), rp, S0);

                S0 = GetK() - V[0];
            }

            Console.WriteLine("Iteration count = " + iter);
            tecplotPrinter.PrintXY(Environment.GetFolderPath(Environment.SpecialFolder.Desktop) + "\\" + "v", 0d,
                GetH(), V, S0);
            double[] KS = GetVKS();
            tecplotPrinter.PrintXYSpecial(GetTau(), 0d, GetRightBoundary(),
                Environment.GetFolderPath(Environment.SpecialFolder.Desktop) + "\\" + "ks_v", 0d, GetH(), GetH(), KS, V,
                S0);

            return Tuple.Create(V, S0);
        }

        public double[] GetExactSolution()
        {
            var arr = new double[GetN1()];
            for (var i = 0; i < arr.Length; ++i)
            {
                var si = 0 + i * GetH();
                arr[i] = GetV(GetSquaredSigma(), GetR(), GetK(), si);
            }

            // V(S) does not touch 0
            // adjust value for better visualization
            arr[0] = arr[1];
            return arr;
        }

        public double[] GetExactSolution(double S0)
        {
            var arr = new double[GetN1()];
            for (var i = 0; i < arr.Length; ++i)
            {
                var si = S0 + i * GetH();
                arr[i] = GetV(GetSquaredSigma(), GetR(), GetK(), si);
            }

            return arr;
        }

        public double GetExactS0()
        {
            return GetK() / (1d + GetSquaredSigma() / (2d * GetR()));
        }

        public double[] GetVKS()
        {
            var res = new double[GetN1()];
            for (var i = 0; i < res.Length; i++)
            {
                res[i] = GetK() - i * GetH();
                if (res[i] < 0d) res[i] = 0d;
            }

            return res;
        }

        private static double GetV(double sigma_sq, double r, double K, double si)
        {
            var p1 = sigma_sq / (2d * r);

            var arg = K / (1 + sigma_sq / (2d * r));
            var pow = (2d * r + sigma_sq) / sigma_sq;
            var p2 = Math.Pow(arg, pow);

            var arg2 = si;
            var pow2 = -2d * r / sigma_sq;
            var p3 = Math.Pow(arg2, pow2);

            var v = p1 * p2 * p3;
            return v;
        }

        /**
         * n - число уравнений (строк матрицы)
         * b - диагональ, лежащая под главной (нумеруется: [1;n-1])
         * c - главная диагональ матрицы A (нумеруется: [0;n-1])
         * d - диагональ, лежащая над главной (нумеруется: [0;n-2])
         * f - правая часть (столбец)
         * x - решение, массив x будет содержать ответ
         */
        // pass n = n + 1!
        private static double[] GetB(int n, double S0, double h, double sigma_sq, double tau)
        {
            var b = new double[n];
            for (var i = 1; i <= n - 1; ++i)
            {
                var si = S0 + i * h;
                var hmh = si - (S0 + (i - 1) * h); // h_{i - 1/2}
                if (hmh * hmh > 4d * tau * sigma_sq * si) throw new ArgumentException("hmh is invalid");

                b[i] = hmh / (4d * tau) - sigma_sq * si * si / (2d * hmh);
            }

            // right boundary cond
            b[n - 1] = 0d;

            return b;
        }

        // главная диагональ матрицы A (нумеруется: [0;n-1])
        // pass n = n + 1!
        private static double[] GetC(int n, double S0, double h, double sigma_sq, double tau, double r)
        {
            var c = new double[n];
            for (var i = 1; i <= n - 1; ++i)
            {
                var si = S0 + i * h;
                var hmh = si - (S0 + (i - 1) * h); // h_{i-1/2}
                var hph = S0 + (i + 1) * h - si; // h_{i+1/2}
                if (hmh * hmh > 4d * tau * sigma_sq * si) throw new ArgumentException("hmh is invalid");

                c[i] = sigma_sq * si * si / (2d * hmh) +
                       sigma_sq * si * si / (2d * hph) +
                       (hmh + hph) *
                       (1d / (4d * tau) + r / 2d);
            }

            // left boundary condition
            var si0 = S0 + 0 * h;
            var hph0 = S0 + (0 + 1) * h - si0; // h_{i+1/2}
            c[0] = sigma_sq * si0 * si0 / (2d * hph0) + hph0 / 2d * r + hph0 / (4d * tau);

            // right boundary condition
            c[n - 1] = 0d;

            return c;
        }

        // диагональ, лежащая над главной (нумеруется: [0;n-2])
        // pass n = n + 1!
        private static double[] GetD(int n, double S0, double h, double sigma_sq, double tau)
        {
            var d = new double[n];
            for (var i = 0; i <= n - 2; ++i)
            {
                var si = S0 + i * h;
                var hph = S0 + (i + 1) * h - si; // h_{i+1/2}
                d[i] = hph / (4d * tau) - sigma_sq * si * si / (2d * hph);
            }

            // left boundary condition
            var si0 = S0 + 0 * h;
            var hph0 = S0 + (0 + 1) * h - si0; // h_{i+1/2}
            d[0] = -(sigma_sq * si0 * si0) / (2d * hph0) + hph0 / (4d * tau);

            // right boundary condition
            d[n - 2] = 0d;
            return d;
        }
    }
}