using System;
using NUnit.Framework;
using TemporalAmericanOption;

namespace PerpetualAmericanOptions
{
    public class TemporalAmericanOptionTests : UnitTestBase
    {
        [Test]
        public void TemporalAmericanOption()
        {
            var parameters = GetParameters();
            var calculator = new TemporalAmericanOptionCalculator(parameters);

            PrintParameters(calculator);
            Console.WriteLine();

            var S0Arr = calculator.Solve();

            Console.WriteLine("Calculated S0");
            for (int i = S0Arr.Length - 1; i >= 0; i--)
            {
                Console.WriteLine("k = " + (i + 1) + " -> " + S0Arr[i]);
            }

            Console.WriteLine();
        }
        
        private TemporalParameters GetParameters()
        {
            var a = 0d;
            var b = 1d;
            var sigma = 1d;
            var tau = 1e-3;
            var n = 400;
            var r = 0.08d;
            var K = 0.5d;
            var M = 50;
            var S0eps = 1e-5;
            return new TemporalParameters(a, b, n, r, tau, sigma, K, S0eps, M);
        }
        
        private static void PrintParameters(TemporalAmericanOptionCalculator calculator)
        {
            Console.WriteLine("b = " + calculator.GetRightBoundary());
            Console.WriteLine("r = " + calculator.GetR());
            Console.WriteLine("N = " + calculator.GetN());
            Console.WriteLine("N_1 = " + calculator.GetN1());
            Console.WriteLine("tau = " + calculator.GetTau());
            Console.WriteLine("sigma = " + calculator.GetSigma());
            Console.WriteLine("sigma_sq = " + calculator.GetSquaredSigma());
            Console.WriteLine("K = " + calculator.GetK());
            Console.WriteLine("M = " + calculator.GetM());
            Console.WriteLine("S0 Eps = " + calculator.GetS0Eps());
        }

        protected override string SetWorkingDir()
        {
            return Environment.GetFolderPath(Environment.SpecialFolder.Desktop) + "\\";
        }
    }
}