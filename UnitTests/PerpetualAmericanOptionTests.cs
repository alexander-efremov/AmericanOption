namespace UnitTests
{
    using System;
    using System.IO;
    using CoreLib;
    using CoreLib.Utils;
    using NUnit.Framework;
    using PerpetualAmericanOptions;

    [TestFixture]
    public class PerpetualAmericanOptionTests : UnitTestBase
    {
        [Test]
        public void PerpetualAmericanOption()
        {
            var parameters = this.GetParameters();
            var calculator = new PerpetualAmericanOptionCalculator(parameters);
            var exactS0 = calculator.GetExactS0();
            Console.WriteLine("Exact S0 = " + exactS0);
            Console.WriteLine();

            this.PrintParameters(calculator);
            Console.WriteLine();

            (double[] item1, var item2) = calculator.Solve();
            double[] exactV = calculator.GetExactSolution(calculator.GetExactS0());
            var l1Error = GetL1Error(calculator, calculator.GetExactSolution(calculator.GetExactS0()), item1);
            var l1Solution = GetL1Solution(calculator, item1);

            Utils.Print(exactV, "V_exact");
            Utils.Print(item1, "V_num");
            Console.WriteLine("S0 = {0}", item2);
            Console.WriteLine("S0 - exactS0 = {0}", Math.Abs(item2 - calculator.GetExactS0()));
            Console.WriteLine("L1 of error = " + l1Error);
            Console.WriteLine("L1 of solution = " + l1Solution);
        }

        [Test]
        public void PerpetualAmericanOptionDrawExactSolution()
        {
            var parameters = this.GetParameters();
            var calculator = new PerpetualAmericanOptionCalculator(parameters);
            double[] V = calculator.GetExactSolution();
            double[] VS0 = calculator.GetVKS();
            var h = (calculator.GetRightBoundary() - calculator.GetExactS0()) / calculator.GetN();
            var printer = new TecplotPrinter(
                0d,
                calculator.GetRightBoundary(),
                calculator.GetTau());
            printer.PrintXY(this.GetWorkingDir() + "exact", 0d, h, VS0, V, 0d);
        }

        [Test]
        public void PerpetualAmericanOptionDrawExactSolutionFromS0()
        {
            var parameters = this.GetParameters();
            var calculator = new PerpetualAmericanOptionCalculator(parameters);
            double[] V = calculator.GetExactSolution(calculator.GetExactS0());
            double[] VS0 = calculator.GetVKS();
            var h = calculator.GetH();

            // var printer = new TecplotPrinterSpecial(
            //    calculator.GetN1(),
            //    0d,
            //    calculator.GetRightBoundary(),
            //    calculator.GetTau());
            TecplotPrinterSpecial.PrintXYSpecial(
                parameters.Tau,
                parameters.A,
                parameters.B,
                this.GetWorkingDir() + "exact-S0",
                0d,
                h,
                h,
                VS0,
                V,
                calculator.GetExactS0());
        }

        [Test]
        public void PerpetualAmericanOptionDrawVS0()
        {
            var parameters = this.GetParameters();
            var calculator = new PerpetualAmericanOptionCalculator(parameters);
            var h = (calculator.GetRightBoundary() - calculator.GetExactS0()) / calculator.GetN();
            double[] exactS0 = calculator.GetVKS();
            var printer = new TecplotPrinter(
                0d,
                calculator.GetRightBoundary(),
                calculator.GetTau());
            printer.PrintXY(this.GetWorkingDir() + "VKS", 0d, h, exactS0);
        }

        private PerpetualParameters GetParameters()
        {
            const double a = 0d;
            const double b = 2d;
            const double sigmaSq = 0.1d;
            const double tau = 1e-3;
            const int n = 400;
            const double r = 0.08d;
            const double K = 0.5d;
            const double S0Eps = 10e-4;
            return new PerpetualParameters(a, b, n, r, tau, sigmaSq, K, S0Eps, this.GetWorkingDir());
        }

        private string GetWorkingDir()
        {
            return Environment.GetFolderPath(Environment.SpecialFolder.Desktop) + Path.DirectorySeparatorChar;
        }
    }
}