using System;
using System.IO;
using System.Security.Cryptography;
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
            var exactV = calculator.GetExactSolution(calculator.GetExactS0());
            var l1Error = GetL1Error(calculator, calculator.GetExactSolution(calculator.GetExactS0()), answer.Item1);
            var l1Solution = GetL1Solution(calculator, answer.Item1);
            
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
            double h = (calculator.GetRightBoundary() - calculator.GetExactS0()) / calculator.GetN();
            var exactS0 = calculator.GetVKS();
            var printer = new TecplotPrinter(calculator.GetN1(),
                0d,
                calculator.GetRightBoundary(),
                calculator.GetTau());
            printer.PrintXY(WorkingDirPath + "VKS", 0d, h, exactS0);
        }

        [Test]
        public void PerpetualAmericanOptionDrawExactSolution()
        {
            var parameters = GetParameters();
            var calculator = new PerpetualAmericanOptionCalculator(parameters);
            var V = calculator.GetExactSolution();
            var VS0 = calculator.GetVKS();
            double h = (calculator.GetRightBoundary() - calculator.GetExactS0()) / calculator.GetN();
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
            var V = calculator.GetExactSolution(calculator.GetExactS0());
            var VS0 = calculator.GetVKS();
            double h = calculator.GetH();
             
            PrintXY(parameters.Tau, parameters.A, parameters.B, WorkingDirPath + "exact-S0", 0d, h, h, VS0, V, calculator.GetExactS0());
        }
        
        internal void PrintXY(double tau, double a, double b, string filename, double t, double h1, double h2, double[] KS, double[] V, double S0 = 0)
        {
            var name = string.Format("{0}_hx={1}_t={2}_tau={3}_a={4}_c={5}.dat", filename, h1, t, tau, a, b);
            using (var writer = new StreamWriter(name, false))
            {
                writer.WriteLine("TITLE = 'DEM DATA'\nVARIABLES = 'x' {0}", "u");
                writer.WriteLine("ZONE T='ONE'");
                writer.WriteLine("I={0} K={1} ZONETYPE=Ordered", KS.Length, 1);
                writer.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE)");
                for (var i = 0; i < KS.Length; i++)
                {
                    writer.WriteLine("{0:e8}  {1:e8}", a + i * h1, KS[i]);
                }

                writer.WriteLine("\nZONE T='TWO'");
                writer.WriteLine("I={0} K={1} ZONETYPE=Ordered", V.Length, 1);
                writer.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE)");
                for (var i = 0; i < V.Length; i++)
                {
                    writer.WriteLine("{0:e8}  {1:e8}", S0 + i * h2, V[i]);
                }
            }
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
        
        internal double GetL1Error(PerpetualAmericanOptionCalculator cal, double[] exact, double[] calculated)
        {
            var err = Utils.GetError(exact, calculated, exact.Length);
            return Utils.GetL1(cal.GetH(), err);
        }
        
        internal double GetL1Solution(PerpetualAmericanOptionCalculator cal, double[] calculatedV)
        {
            return Utils.GetL1(cal.GetH(), calculatedV);
        }

        private static void PrintParameters(PerpetualAmericanOptionCalculator calculator)
        {
            Console.WriteLine("b = " + calculator.GetRightBoundary());
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