using System;
using NUnit.Framework;

// ReSharper disable CommentTypo

namespace PerpetualAmericanOptions
{
    internal class PerpetualAmericanOptionCalculator
    {
        private readonly ThomasAlgorithmCalculator _thomasAlgorithmCalculator;
        private readonly int n;
        private readonly int n_1;
        private readonly double a;
        private readonly double b;
        private readonly double sigma;
        private readonly double sigma_sq;
        private readonly double tau;
        private readonly double r;
        private readonly double h;
        private readonly double h_2;
        private readonly double h_sq;
        private readonly double K;
        
        public int GetN()
        {
            return n;
        }
        
        public int GetN1()
        {
            return n_1;
        }

        public double GetLeftBoundary()
        {
            return a;
        }
        
        public double GetRightBoundary()
        {
            return b;
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
        
        public double GetK()
        {
            return K;
        }

        public PerpetualAmericanOptionCalculator(Parameters parameters)
        {
            a = parameters.A;
            b = parameters.B;
            sigma = parameters.Sigma;
            sigma_sq = sigma * sigma;
            tau = parameters.Tau;
            n = parameters.N;
            n_1 = n + 1;
            r = parameters.R;
            K = parameters.K;

            h = (b - a) / n;
            h_2 = h * 0.5d;
            h_sq = h * h;

            CheckParameters();

            _thomasAlgorithmCalculator = new ThomasAlgorithmCalculator(n_1);
        }
        
        public Tuple<double[], double> Solve()
        {
            var tecplotPrinter = new TecplotPrinter(GetN1(),
                GetH(), 
                GetLeftBoundary(),
                GetRightBoundary(),
                GetTau());
            var thomasArrayPrinter = new ThomasArrayPrinter();
            double[] calculatedV = new double[n_1];
            // TODO:TEMPORARILY!
            double S0 = GetExactS0();
            var prevV = GetExactSolution();
            //while (Math.Abs(GetExactS0() - S0) > 10e-5)
            {
                var b_t = GetB(n_1, S0, h, h_sq, sigma_sq, tau, r);
                var c_t = GetC(n_1, S0, h, h_sq, sigma_sq, tau, r);
                var d_t = GetD(n_1, S0, h, h_sq, sigma_sq, tau, r);
                thomasArrayPrinter.PrintThomasArrays(b_t, c_t, d_t);
                
                var rp = CalculateRightPart(prevV, S0);
                Utils.Print(rp, "rp");
                
                calculatedV = _thomasAlgorithmCalculator.Calculate(b_t, c_t, d_t, rp);

                tecplotPrinter.PrintXY(Environment.GetFolderPath(Environment.SpecialFolder.Desktop) + "\\" + "rp", 0d, rp);
                tecplotPrinter.PrintXY(Environment.GetFolderPath(Environment.SpecialFolder.Desktop) + "\\" + "v", 0d, calculatedV);

                var err = Utils.FillArrayDiff(GetExactSolution(), rp);
                tecplotPrinter.PrintXY("ex_minus_rp", 0d, err);
                
                Array.Copy(calculatedV, prevV, n_1);
            }

            return Tuple.Create(calculatedV, S0);
        }
        
        public double GetL1Error(double[] calculated)
        {
            var exact = GetExactSolution();
            var err = Utils.GetError(exact, calculated, n_1);
            return Utils.GetL1(h, n_1, err);
        }
        
        public double GetL1Solution(double[] calculatedV)
        {
            return Utils.GetL1(h, n_1, calculatedV);
        }
        
        public double[] GetExactSolution()
        {
            var arr = new double[n_1];
            
            for (var i = 1; i < n_1; ++i)
            {
                double si = a + i * h;
                arr[i] = GetV(sigma_sq, r, K, si);
            }
            //todo: BOUNDARY CONDITION!
            arr[0] = arr[1];
            return arr;
        }
        
        public double GetExactS0()
        {
            return K / (1d + (sigma_sq / (2d * r)));
        }

        public double[] GetVKS()
        {
            var res = new double[n_1];
            for (int i = 0; i < res.Length; i++)
            {
                res[i] = K - (i * h);
                if (res[i] < 0d)
                {
                    res[i] = 0d;
                }
            }

            return res;
        }
        
        private static double GetV(double sigma_sq, double r, double K, double si)
        {
            var p1 = sigma_sq / (2d * r);
            
            var arg = ( K / (1 + (sigma_sq / (2d * r))) );
            var pow = (2d * r + sigma_sq) / sigma_sq;
            var p2 = Math.Pow(arg, pow);
            
            var arg2 = si;
            var pow2 = (-2d * r) / sigma_sq;
            var p3 = Math.Pow(arg2, pow2);

            var v = p1 * p2 * p3;
            return v;
        }

        private double CalculateBetas(int i, double[] v, double S0)
        {
            var Si = S0 + i * h;
            var s_m_h = Si - h_2;
            var s_p_h = Si + h_2;
            
            var leftArg = tau / h;
            var rightArg = (s_p_h / (2d * r));
            Assert.LessOrEqual(leftArg, rightArg);

            double beta1 = 1d / (8d * tau) * (1d + (2d * tau * r) / (h * s_m_h)) *
                        (1d + (2d * tau * r) / h * s_m_h);

            double beta2 = 1d / (8d * tau) * (3d + (2d * tau * r) / (h * s_m_h)) *
                        (1d - (2d * tau * r) / (h * s_m_h))
                        + 1d / (8d * tau) * (3d - (2d * tau * r) / (h * s_p_h)) *
                        (1d + (2d * tau * r) / (h * s_p_h));

            double beta3 = 1d / (8d * tau) * (1d - (2d * tau * r) / (h * s_p_h)) *
                        (1d - (2d * tau * r) / (h * s_p_h));
            
            double val = beta1 * v[i - 1] + beta2 * v[i] + beta3 * v[i + 1];
            return val;
        }

        private double[] CalculateRightPart(double[] v, double S0)
        {
            var rp = new double[n_1];
            for (var i = 1; i < n_1 - 1; ++i)
            {
                rp[i] = CalculateBetas(i, v, S0);
            }
            // todo: take attention! boundary conditions
            rp[0] = rp[1];
            rp[n_1 - 1] = rp[n_1 - 2];
            return rp;
        }
        
        /**
         * n - число уравнений (строк матрицы)
         * b - диагональ, лежащая под главной (нумеруется: [1;n-1])
         * c - главная диагональ матрицы A (нумеруется: [0;n-1])
         * d - диагональ, лежащая над главной (нумеруется: [0;n-2])
         * f - правая часть (столбец)
         * x - решение, массив x будет содержать ответ
         */
        private static double[] GetB(int n_1, double S0, double h, double h_sq, double sigma_sq, double tau, double r)
        {
            var b = new double[n_1];
            b[0] = b[1] = b[n_1 - 1] = 0d;
            for (var i = 2; i < n_1 - 1; ++i)
            {
                double si = S0 + i * h;
                double si_p_h = si + 0.5 * h;
                var leftArg = tau / h;
                var rightArg = (r * h) / (4d * sigma_sq * si);
                Assert.GreaterOrEqual(leftArg, rightArg, "leftArg >= rightArg GetB");
                Assert.LessOrEqual(h, (2d * si * si_p_h * sigma_sq) / (r * r));
                b[i] = ( ((-sigma_sq) / (2d * h_sq)) + (r / (8d * tau * si)) );
            }

            return b;
        }

        private static double[] GetC(int n_1, double S0, double h, double h_sq, double sigma_sq, double tau, double r)
        {
            var c = new double[n_1];
            c[0] = c[n_1 - 1] = 0d;
            for (var i = 1; i < n_1 - 1; ++i)
            {
                double si = S0 + i * h;
                double si_p_h = si + 0.5 * h;
                double si_sq = si * si;
                var leftArg = tau / h;
                var rightArg = (r * h) / (4d * sigma_sq * si);
                Assert.GreaterOrEqual(leftArg, rightArg, "leftArg >= rightArg GetC");
                Assert.LessOrEqual(h, (2d * si * si_p_h * sigma_sq) / (r * r));
                c[i] = ((3d * r) / (4d * tau * si)) + (r / si_sq) + sigma_sq / h_sq;
            }

            return c;
        }

        private static double[] GetD(int n_1, double S0, double h, double h_sq, double sigma_sq, double tau, double r)
        {
            var d = new double[n_1];
            d[0] = d[n_1 - 2] = d[n_1 - 1] = 0d;
            for (var i = 1; i < n_1 - 2; ++i)
            {
                double si = S0 + i * h;
                double si_p_h = si + 0.5 * h;
                var leftArg = tau / h;
                var rightArg = (r * h) / (4d * sigma_sq * si);
                Assert.GreaterOrEqual(leftArg, rightArg, "leftArg >= rightArg GetC");
                Assert.LessOrEqual(h, (2d * si * si_p_h * sigma_sq) / (r * r));
                d[i] = ( ((-sigma_sq) / (2d * h_sq)) + (r / (8d * tau * si)) );
            }

            return d;
        }

        private void CheckParameters()
        {
            Assert.True(h > 0.0);
            Assert.AreEqual(h * h, h_sq);
            Assert.AreEqual(h * 0.5, h_2);
            Assert.True(n > 1);
            Assert.AreEqual(n + 1, n_1);
            Assert.AreEqual(0, a);
            Assert.True(tau > 0d);
            Assert.True(K > 0d);
            
            // let's check scheme t/h restrictions
            // left for Thomas algo, right for postiive defined M matrix
            // (r * h) / (4*sigmaSq*Si) <= tau/h <= S_i+1/2 / (2 * r)
            // Si = a + i*h S_i+1/2 = a + i*h + h/2
            // then 
            // (r * h) / (4*sigmaSq*(a+i*h)) <= tau/h <= (a+i*h + 0.5h) / (2 * r)
            // min(Si) = h (except i = 0) and min(Si+1/2) = 0.5h
            // then
            // (r * h) / (4*sigmaSq*h) <= tau/h <= 0.5h / (2 * r)
            // 0.25r / sigmaSq <= tau/h <= 0.25h / r

            Assert.GreaterOrEqual(tau / h, (0.25d * r) / sigma_sq);
            
            Assert.GreaterOrEqual((0.25d * h) / r, tau / h);
        }
    }

    [TestFixture]
    public class PerpetualAmericanOptionTests : UnitTestBase
    {
        protected override string SetWorkingDir()
        {
            return Environment.GetFolderPath(Environment.SpecialFolder.Desktop) + "\\";
        }

        [Test]
        public void PerpetualAmericanOption()
        {
            var parameters = GetParameters();
            var calculator = new PerpetualAmericanOptionCalculator(parameters);
            var exactS0 = calculator.GetExactS0();
            Console.WriteLine("Exact S0 = " + exactS0);
            Console.WriteLine();

            PrintParameters(calculator);
            Console.WriteLine();
            
            var answer = calculator.Solve();
            var exactV = calculator.GetExactSolution();
            var l1Error = calculator.GetL1Error(answer.Item1);
            var l1Solution = calculator.GetL1Solution(answer.Item1);
            
            var tecplotPrinter = new TecplotPrinter(calculator.GetN1(),
                calculator.GetH(), 
                calculator.GetLeftBoundary(),
                calculator.GetRightBoundary(),
                calculator.GetTau());
            tecplotPrinter.PrintXY("exact_number", 0d, exactV, answer.Item1);
            tecplotPrinter.PrintXY("numerical", 0d, answer.Item1);
            
            Utils.Print(exactV, "V_exact");
            Utils.Print(answer.Item1, "V_num");
            Console.WriteLine("S0 = {0}", answer.Item2);
            Console.WriteLine("L1 of error = " + l1Error);
            Console.WriteLine("L1 of solution = " + l1Solution);
        }

        [Test]
        public void PerpetualAmericanOptionDrawVS0()
        {
            var parameters = GetParameters();
            var calculator = new PerpetualAmericanOptionCalculator(parameters);
            var exactS0 = calculator.GetVKS();
            var printer = new TecplotPrinter(calculator.GetN1(),
                calculator.GetH(), 
                calculator.GetLeftBoundary(),
                calculator.GetRightBoundary(),
                calculator.GetTau());
            printer.PrintXY(WorkingDirPath + "VKS", 0d, exactS0);
        }

        [Test]
        public void PerpetualAmericanOptionDrawExactSolution()
        {
            var parameters = GetParameters();
            var calculator = new PerpetualAmericanOptionCalculator(parameters);
            var V = calculator.GetExactSolution();
            var VS0 = calculator.GetVKS();
            var printer = new TecplotPrinter(calculator.GetN1(),
                calculator.GetH(), 
                calculator.GetLeftBoundary(),
                calculator.GetRightBoundary(),
                calculator.GetTau());
            printer.PrintXY(WorkingDirPath + "exact", 0d, VS0, V);
        }

        private Parameters GetParameters()
        {
            double a = 0d;
            double b = 1d;
            double sigma = 1d;
            double tau = 2e-4;
            int n = 100;
            double r = 0.08d;
            double K = 0.5d;
            return new Parameters(a, b, n, r, tau , sigma, K);
        }

        private static void PrintParameters(PerpetualAmericanOptionCalculator calculator)
        {
            Console.WriteLine("a = " + calculator.GetLeftBoundary());
            Console.WriteLine("b = " + calculator.GetRightBoundary());
            Console.WriteLine("h = " + calculator.GetH());
            Console.WriteLine("h_2 = " + calculator.GetHalfH());
            Console.WriteLine("h_sq = " + calculator.GetSquaredH());
            Console.WriteLine("r = " + calculator.GetR());
            Console.WriteLine("N = " + calculator.GetN());
            Console.WriteLine("N_1 = " + calculator.GetN1());
            Console.WriteLine("tau = " + calculator.GetTau());
            Console.WriteLine("sigma = " + calculator.GetSigma());
            Console.WriteLine("sigma_sq = " + calculator.GetSquaredSigma());
            Console.WriteLine("K = " + calculator.GetK());
        }
    }
}