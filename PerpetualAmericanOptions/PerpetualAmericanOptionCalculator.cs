using System;
using NUnit.Framework;

// ReSharper disable CommentTypo

namespace PerpetualAmericanOptions
{
    internal class PerpetualAmericanOptionCalculator
    {
        private readonly ThomasAlgorithmCalculator _thomasAlgorithmCalculator;
        private readonly int M;
        private readonly int n;
        private readonly int n_1;
        private readonly double a;
        private readonly double A;
        private readonly double b;
        private readonly double sigma;
        private readonly double sigma_sq;
        private readonly double tau;
        private readonly double r;
        private readonly double h;
        private readonly double h_2;
        private readonly double h_sq;
        
        public int GetN()
        {
            return n;
        }
        
        public int GetN1()
        {
            return n_1;
        }
        
        public int GetM()
        {
            return M;
        }

        public double GetLeftBoundary()
        {
            return a;
        }
        
        public double GetRightBoundary()
        {
            return b;
        }
        
        public double GetA()
        {
            return A;
        }
        
        public double GetSigma()
        {
            return sigma;
        }
        
        public double GetSquaredSigma()
        {
            return sigma_sq;
        }
        
        public double GetTau()
        {
            return tau;
        }
        
        public double GetR()
        {
            return r;
        }
        
        public double GetH()
        {
            return h;
        }
        
        public double GetHalfH()
        {
            return h_2;
        }
        
        public double GetSquaredH()
        {
            return h_sq;
        }

        public PerpetualAmericanOptionCalculator(double b, double a, double sigma, double tau, double A, int n, int M,
            double r)
        {
            this.a = a;
            this.b = b;
            this.sigma = sigma;
            sigma_sq = sigma * sigma;
            this.tau = tau;
            this.A = A;
            this.n = n;
            n_1 = n + 1;
            this.M = M;
            this.r = r;

            h = (this.b - this.a) / this.n;
            h_2 = h * 0.5;
            h_sq = h * h;

            CheckParameters();

            _thomasAlgorithmCalculator = new ThomasAlgorithmCalculator(n_1, h_sq, this.tau, this.sigma);
        }
        
        public double[] Solve()
        {
            var calculatedV = new double[n_1];
            var rp = new double [n_1];

            var b_t = _thomasAlgorithmCalculator.GetB();
            var c_t = _thomasAlgorithmCalculator.GetC();
            var d_t = _thomasAlgorithmCalculator.GetD();
            
            var Vpr = GetExactSolution(tau);
            Utils.Print(Vpr, "V_PR");

            for (var k = 0; k <= M; ++k)
            {
                FillRightPart(rp, Vpr, tau * k);
            }
            //Assert.GreaterOrEqual(tau/h, (r*h)/(4.0*sigma_sq*S[i]));
            //Assert.LessOrEqual(tau/h, (S[i-1/2]/h)/(2.0*tau));
            //print_XY("rp", n_1, h, M, a, b, tau, rp);
            //thomas_algo_verzh_modified(n_1, b_t, c_t, d_t, rp, u);
            //print_XY("u", n_1, h, M, a, b, tau, u);
            //memcpy(u_pr, u, n_1);

            //Vex = GetExactSolution(tau * M);
            // fill_arr_diff(err, u_ex, rp, n_1);
            // print_XY("ex_minus_rp", n_1, h, M, a, b, tau, err);

            //Console.WriteLine("numerical arr sum = %le\n", calc_array_sum(u, n_1, 1, false));
            //Console.WriteLine("exact arr sum = %le\n", calc_array_sum(u_ex, n_1, 1, false));

            //Utils.PrintXY("exact_number", n_1, h, M, a, b, tau, Vex, calculatedV);
            //print_XY("exact", n_1, h, M, a, b, tau, u_ex);
            //Utils.PrintXY("numerical", n_1, h, M, a, b, tau, calculatedV);
            //print_vector("err", err, n_1);
            //Utils.FillArrayDiff(err, Vex, calculatedV, n_1);
            //print_XY("err", n_1, h, M, a, b, tau, err);
            //Console.WriteLine("uniform norm of l1 norm of error %22.014le\n", get_l1_norm(n_1, err));

            return calculatedV;
        }
        
        public double GetL1Error(double[] calculated)
        {
            var exact = GetExactSolution(tau * M);
            double[] err = Utils.GetError(exact, calculated, n_1);
            return Utils.GetL1(h, n_1, err);
        }
        
        public double GetL1Solution(double[] calculatedV)
        {
            return Utils.GetL1(h, n_1, calculatedV);
        }
        
        public double[] GetExactSolution(double time)
        {
            double[] arr = new double[n+1];
            for (var i = 0; i <= n; ++i)
            {
                arr[i] = GetV(A, time, a + i * h);
            }

            return arr;
        }
        
        private double GetV(double inA, double t, double x)
        {
            return inA * Math.Sin(x * (t + 1.0)) * (1.0 - x * x);
        }

        private double get_a(double inA, double t, double x)
        {
            if (x <= 0.0)
            {
                return 0.0;
            }
            return inA * Math.Sin(x * (t + 1.0)) * (1.0 - x * x);
        }

        private double calc_by_betas(int ii, double[] u_pr, double time)
        {
            Assert.IsTrue(h_sq <= 8.0 * sigma * tau);

            var xi = a + ii * h;

            var alpha1 = get_a(A, time, xi - h_2);
            var alpha2 = get_a(A, time, xi + h_2);

            var beta1 = 1.0 / (8.0 * tau) * (1.0 + 2.0 * tau * alpha1 / h) *
                        (1.0 + 2.0 * tau * alpha1 / h);

            var beta2 = 1.0 / (8.0 * tau) * (3.0 - 2.0 * tau * alpha1 / h) *
                        (1.0 + 2.0 * tau * alpha1 / h)
                        + 1.0 / (8.0 * tau) * (3.0 + 2.0 * tau * alpha2 / h) *
                        (1.0 - 2.0 * tau * alpha2 / h);

            var beta3 = 1.0 / (8.0 * tau) * (1.0 - 2.0 * tau * alpha2 / h) *
                        (1.0 - 2.0 * tau * alpha2 / h);

            var val = beta1 * u_pr[ii - 1] + beta2 * u_pr[ii] + beta3 * u_pr[ii + 1];

            return val;
        }

        private void AssertACondition(double A, double time)
        {
            var max = -1.0e-16;
            for (var j = 0; j < n_1; ++j)
            {
                var alpha = get_a(A, time, j * h - h_2);
                if (Math.Abs(alpha) > max) max = Math.Abs(alpha);
            }

            if (tau * max >= h_2)
            {
                Console.WriteLine("!!!!!! tau * Math.Abs(max_alpha) < h_2 FAILED\n");
            }
        }

        private void FillRightPart(double[] rp, double[] u_pr, double time)
        {
            rp[0] = rp[n_1 - 1] = 0.0;

            AssertACondition(A, time);

            for (var i = 1; i < n_1 - 1; ++i)
            {
                rp[i] = calc_by_betas(i, u_pr, time);
            }
        }

        private void CheckParameters()
        {
            Assert.True(h > 0.0);
            Assert.AreEqual(h * h, h_sq);
            Assert.AreEqual(h * 0.5, h_2);
            Assert.True(n > 1);
            Assert.AreEqual(n + 1, n_1);
            Assert.AreEqual(0, a);
            Assert.True(tau > 0.0);
            Assert.True(A > 0.0);
            Assert.True(M >= 1);
        }
    }

    [TestFixture]
    public class PerpetualAmericanOptionTests
    {
        [Test]
        public void PerpetualAmericanOption()
        {
            var lb = 0.0;
            var rb = 1.0;
            var sigma = 1.0;
            var A = 1.0;
            var tau = 4.0e-3;
            var n = 100;
            var M = 10;
            var r = 1.0;

            var calculator = new PerpetualAmericanOptionCalculator(lb, rb, sigma, tau, A, n, M, r);
            PrintParameters(calculator);
            var calculatedV = calculator.Solve();
            var exactV = calculator.GetExactSolution(calculator.GetM() * calculator.GetTau());
            var l1Error = calculator.GetL1Error(calculatedV);
            var l1Solution = calculator.GetL1Solution(calculatedV);
            
            var tecplotPrinter = new TecplotPrinter(calculator.GetN1(),
                calculator.GetH(), 
                calculator.GetLeftBoundary(),
                calculator.GetRightBoundary(),
                calculator.GetTau());
            tecplotPrinter.PrintXY("exact_number", calculator.GetM() * calculator.GetTau(), exactV, calculatedV);
            tecplotPrinter.PrintXY("numerical", calculator.GetM() * calculator.GetTau(), calculatedV);
            
            Console.WriteLine(l1Error);
            Console.WriteLine(l1Solution);
        }

        private void PrintParameters(PerpetualAmericanOptionCalculator calculator)
        {
            Console.WriteLine("A = " + calculator.GetA());
            Console.WriteLine("a = " + calculator.GetLeftBoundary());
            Console.WriteLine("b = " + calculator.GetRightBoundary());
            Console.WriteLine("h = " + calculator.GetH());
            Console.WriteLine("h_2 = " + calculator.GetHalfH());
            Console.WriteLine("h_sq = " + calculator.GetSquaredH());
            Console.WriteLine("r = " + calculator.GetR());
            Console.WriteLine("M = " + calculator.GetM());
            Console.WriteLine("N = " + calculator.GetN());
            Console.WriteLine("N_1 = " + calculator.GetN1());
            Console.WriteLine("tau = " + calculator.GetTau());
            Console.WriteLine("sigma = " + calculator.GetSigma());
            Console.WriteLine("sigma_sq = " + calculator.GetSquaredSigma());
        }
    }
    
    [TestFixture]
    public class ThomasAlgorithmTests
    {
        [Test]
        public void PerpetualAmericanOption()
        {
            var lb = 0.0;
            var rb = 1.0;
            var sigma = 1.0;
            var tau = 4.0e-3;
            var n = 100;
            var h_sq = (rb-lb) / n;
            h_sq = h_sq * h_sq;

            var calculator = new ThomasAlgorithmCalculator(n+1, h_sq, tau, sigma);
            double[] b = calculator.GetB();
            double[] c = calculator.GetC();
            double[] d = calculator.GetD();
            PrintThomasArrays(b, c, d);
        }

        private void PrintThomasArrays(double[] b, double[] c, double[] d)
        {
            Console.WriteLine("B:");
            for (int i = 0; i < b.Length; i++)
            {
                Console.Write(b[i] + " ");
            }

            Console.WriteLine();
            Console.WriteLine("C:");
            for (int i = 0; i < c.Length; i++)
            {
                Console.Write(c[i] + " ");
            }

            Console.WriteLine();
            Console.WriteLine("D:");
            for (int i = 0; i < d.Length; i++)
            {
                Console.Write(d[i] + " ");
            }
        }
    }
}