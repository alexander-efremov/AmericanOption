using System;
using NUnit.Framework;

namespace PerpetualAmericanOptions
{
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
            var l1Error = GetL1Error(calculator, answer.Item1);
            var l1Solution = GetL1Solution(calculator, answer.Item1);
            
            var tecplotPrinter = new TecplotPrinter(calculator.GetN1(),
                calculator.GetH(), 
                0d,
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
                0d,
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
                0d,
                calculator.GetRightBoundary(),
                calculator.GetTau());
            printer.PrintXY(WorkingDirPath + "exact", 0d, VS0, V);
        }

        private Parameters GetParameters()
        {
            double a = 0d;
            double b = 1d;
            double sigma = 1d;
            double tau = 1e-3;
            int n = 400;
            double r = 0.08d;
            double K = 0.5d;
            return new Parameters(a, b, n, r, tau , sigma, K);
        }
        
        internal double GetL1Error(PerpetualAmericanOptionCalculator cal, double[] calculated)
        {
            var exact = cal.GetExactSolution();
            var err = Utils.GetError(exact, calculated, exact.Length);
            return Utils.GetL1(cal.GetH(), err);
        }
        
        internal double GetL1Solution(PerpetualAmericanOptionCalculator cal, double[] calculatedV)
        {
            return Utils.GetL1(cal.GetH(), calculatedV);
        }

        private static void PrintParameters(PerpetualAmericanOptionCalculator calculator)
        {
            //Console.WriteLine("a = " + calculator.GetLeftBoundary());
            Console.WriteLine("b = " + calculator.GetRightBoundary());
            Console.WriteLine("h = " + calculator.GetH());
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