// ReSharper disable All
namespace PerpetualAmericanOptions
{
    using System;

    using CoreLib;

    // from new presentation with FEM
    public class PerpetualAmericanOptionCalculator : AmericanOptionCalculatorBase
    {
        public PerpetualAmericanOptionCalculator(Parameters parameters)
            : base(parameters)
        {
        }

        public double GetExactS0()
        {
            return this.GetK() / (1d + this.GetSquaredSigma() / (2d * this.GetR()));
        }

        public double[] GetExactSolution()
        {
            var arr = new double[this.GetN1()];
            for (var i = 0; i < arr.Length; ++i)
            {
                var si = 0 + i * this.GetH();
                arr[i] = GetV(this.GetSquaredSigma(), this.GetR(), this.GetK(), si);
            }

            // V(S) does not touch 0
            // adjust value for better visualization
            arr[0] = arr[1];
            return arr;
        }

        public double[] GetExactSolution(double S0)
        {
            var arr = new double[this.GetN1()];
            for (var i = 0; i < arr.Length; ++i)
            {
                var si = S0 + i * this.GetH();
                arr[i] = GetV(this.GetSquaredSigma(), this.GetR(), this.GetK(), si);
            }

            return arr;
        }

        public double[] GetVKS()
        {
            var res = new double[this.GetN1()];
            for (var i = 0; i < res.Length; i++)
            {
                res[i] = this.GetK() - i * this.GetH();
                if (res[i] < 0d)
                {
                    res[i] = 0d;
                }
            }

            return res;
        }

        public Tuple<double[], double> Solve()
        {
            var tecplotPrinter = new TecplotPrinterSpecial(0d, this.GetRightBoundary(), this.GetTau());

            double[] V = new double[this.GetN()];
            var S0 = this.GetK();
            var iter = 0;
            while (Math.Abs(this.GetExactS0() - S0) > this.GetS0Eps())
            {
                iter++;
                this.UpdateH(S0);
                Console.WriteLine("S0 = " + S0);
                Console.WriteLine("h = " + this.GetH());
                double[] b_t = GetB(this.GetN1(), S0, this.GetH(), this.GetSquaredSigma(), this.GetTau());
                double[] c_t = GetC(this.GetN1(), S0, this.GetH(), this.GetSquaredSigma(), this.GetTau(), this.GetR());
                double[] d_t = GetD(this.GetN1(), S0, this.GetH(), this.GetSquaredSigma(), this.GetTau());
                double[] rp = this.CalculateRightPart(S0);
                V = this.ThomasAlgorithmCalculator.Calculate(b_t, c_t, d_t, rp);

                // var printer = new ThomasArrayPrinter();
                // printer.PrintThomasArrays(b_t, c_t, d_t);
                tecplotPrinter.PrintXY(this.GetWorkDir() + "perpetual-rp", 0d, this.GetH(), rp, S0);

                S0 = this.GetK() - V[0];
            }

            Console.WriteLine("Iteration count = " + iter);
            tecplotPrinter.PrintXY(this.GetWorkDir() + "v", 0d, this.GetH(), V, S0);
            double[] KS = this.GetVKS();
            TecplotPrinterSpecial.PrintXYSpecial(this.GetTau(), 0d, this.GetRightBoundary(), this.GetWorkDir() + "ks_v", 0d, this.GetH(), this.GetH(), KS, V, S0);

            return Tuple.Create(V, S0);
        }

        private static void CheckHCorrectness(double h, double tau, double sph, double r)
        {
            if (tau / h > 1d / (2d * r * sph))
            {
                throw new InvalidOperationException("tau / h");
            }
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
                if (hmh * hmh > 4d * tau * sigma_sq * si)
                {
                    throw new ArgumentException("hmh is invalid");
                }

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
                if (hmh * hmh > 4d * tau * sigma_sq * si)
                {
                    throw new ArgumentException("hmh is invalid");
                }

                c[i] = sigma_sq * si * si / (2d * hmh) + sigma_sq * si * si / (2d * hph) +
                       (hmh + hph) * (1d / (4d * tau) + r / 2d);
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

        private double[] CalculateRightPart(double S0)
        {
            var V = new double[this.GetN1() + 1]; // add + 1 to calculate right boundary
            for (var i = 0; i <= this.GetN1(); i++)
            {
                var si = S0 + i * this.GetH();
                V[i] = GetV(this.GetSquaredSigma(), this.GetR(), this.GetK(), si);
            }

            var rp = new double[this.GetN1()];
            var si0 = S0 + 0 * this.GetH();
            var hph0 = S0 + (0d + 1d) * this.GetH() - si0; // h_{i+1/2}
            var hmh0 = si0 - (S0 + (0d - 1d) * this.GetH()); // h_{i-1/2}
            var sph0 = si0 + 0.5d * hph0; // s_{i+1/2}
            var smh0 = si0 - 0.5d * hmh0; // s_{i-1/2}
            CheckHCorrectness(hph0, this.GetTau(), sph0, this.GetR());
            CheckHCorrectness(hmh0, this.GetTau(), smh0, this.GetR());
            var beta20 = 1d / (8d * this.GetTau()) * (3d - 2d * this.GetTau() * this.GetR() * sph0 / hph0) * (1d + 2d * this.GetTau() * this.GetR() * sph0 / hph0);
            var beta30 = 1d / (8d * this.GetTau()) * (1d - 2d * this.GetTau() * this.GetR() * sph0 / hph0) * (1d - 2d * this.GetTau() * this.GetR() * sph0 / hph0);
            rp[0] = hph0 * (beta20 * V[0] + beta30 * V[1]) + (this.GetSquaredSigma() * si0 * si0) / 2d;

            for (var i = 1; i <= this.GetN1() - 1; i++)
            {
                var si = S0 + i * this.GetH();
                var hmh = si - (S0 + (i - 1) * this.GetH()); // h_{i-1/2}
                var hph = S0 + (i + 1) * this.GetH() - si; // h_{i+1/2}
                var smh = si - 0.5d * hmh; // s_{i-1/2}
                var sph = si + 0.5d * hph; // s_{i+1/2}
                CheckHCorrectness(hph, this.GetTau(), sph, this.GetR());
                CheckHCorrectness(hmh, this.GetTau(), sph, this.GetR());
                var beta1 = 1d / (8d * this.GetTau()) * (1d + 2d * this.GetTau() * this.GetR() * smh / this.GetH())
                                                      * (1d + 2d * this.GetTau() * this.GetR() * smh / this.GetH());
                var beta2 =
                    1d / (8d * this.GetTau()) * (3d + 2d * this.GetTau() * this.GetR() * smh / this.GetH()) * (1d - 2d * this.GetTau() * this.GetR() * smh / this.GetH())
                    + 1d / (8d * this.GetTau()) * (3d - 2d * this.GetTau() * this.GetR() * sph / this.GetH())
                                                * (1d + 2d * this.GetTau() * this.GetR() * sph / this.GetH());
                var beta3 = 1d / (8d * this.GetTau()) * (1d - 2d * this.GetTau() * this.GetR() * sph / this.GetH())
                                                      * (1d - 2d * this.GetTau() * this.GetR() * sph / this.GetH());
                var val = ((hph + hmh) / 2d) * (beta1 * V[i - 1] + beta2 * V[i] + beta3 * V[i + 1]);

                rp[i] = this.GetF() + val;
            }

            return rp;
        }

        // it is zero for our case
        private double GetF()
        {
            return 0d;
        }
    }
}