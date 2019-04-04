using System;
using CoreLib;
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

        private PerpetualParameters GetParameters()
        {
            var a = 0d;
            var b = 1d;
            var sigma = 1d;
            var tau = 1e-3;
            var n = 400;
            var r = 0.08d;
            var K = 0.5d;
            var S0Eps = 10e-4;
            return new PerpetualParameters(a, b, n, r, tau, sigma, K, S0Eps);
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

            Tuple<double[], double> answer = calculator.Solve();
            double[] exactV = calculator.GetExactSolution(calculator.GetExactS0());
            var l1Error = GetL1Error(calculator, calculator.GetExactSolution(calculator.GetExactS0()), answer.Item1);
            var l1Solution = GetL1Solution(calculator, answer.Item1);

            Utils.Print(exactV, "V_exact");
            Utils.Print(answer.Item1, "V_num");
            Console.WriteLine("S0 = {0}", answer.Item2);
            Console.WriteLine("S0 - exactS0 = {0}", Math.Abs(answer.Item2 - calculator.GetExactS0()));
            Console.WriteLine("L1 of error = " + l1Error);
            Console.WriteLine("L1 of solution = " + l1Solution);
        }

        [Test]
        public void PerpetualAmericanOptionDrawExactSolution()
        {
            var parameters = GetParameters();
            var calculator = new PerpetualAmericanOptionCalculator(parameters);
            double[] V = calculator.GetExactSolution();
            double[] VS0 = calculator.GetVKS();
            var h = (calculator.GetRightBoundary() - calculator.GetExactS0()) / calculator.GetN();
            var printer = new TecplotPrinter(calculator.GetN1(),
                0d,
                calculator.GetRightBoundary(),
                calculator.GetTau());
            printer.PrintXY(WorkingDirPath + "exact", 0d, h, VS0, V);
        }

        [Test]
        public void PerpetualAmericanOptionDrawExactSolutionFromS0()
        {
            var parameters = GetParameters();
            var calculator = new PerpetualAmericanOptionCalculator(parameters);
            double[] V = calculator.GetExactSolution(calculator.GetExactS0());
            double[] VS0 = calculator.GetVKS();
            var h = calculator.GetH();
            var printer = new TecplotPrinterSpecial(calculator.GetN1(),
                0d,
                calculator.GetRightBoundary(),
                calculator.GetTau());
            printer.PrintXYSpecial(parameters.Tau, parameters.A, parameters.B, WorkingDirPath + "exact-S0", 0d, h, h,
                VS0, V, calculator.GetExactS0());
        }

        [Test]
        public void PerpetualAmericanOptionDrawVS0()
        {
            var parameters = GetParameters();
            var calculator = new PerpetualAmericanOptionCalculator(parameters);
            var h = (calculator.GetRightBoundary() - calculator.GetExactS0()) / calculator.GetN();
            double[] exactS0 = calculator.GetVKS();
            var printer = new TecplotPrinter(calculator.GetN1(),
                0d,
                calculator.GetRightBoundary(),
                calculator.GetTau());
            printer.PrintXY(WorkingDirPath + "VKS", 0d, h, exactS0);
        }
    }
}