namespace Runner
{
    using System;
    using System.IO;
    using AmericanOption;
    using AmericanOptions;

    internal class Program
    {
        public static void Main(string[] args)
        {
            var parameters = GetParameters(true, GetWorkingDir() + "AO/");
            var calculator = new AmericanOptionCalculator(parameters, true, false);
            PrintParameters(calculator);
            double[] S0t = calculator.Solve();
        }
        
        private static string GetWorkingDir()
        {
            return Environment.GetFolderPath(Environment.SpecialFolder.Desktop) + Path.DirectorySeparatorChar;
        }
        
        private static void PrintParameters(AmericanOptionCalculator calculator)
        {
            Console.WriteLine("a = " + calculator.GetLeftBoundary());
            Console.WriteLine("b = " + calculator.GetRightBoundary());
            Console.WriteLine("r = " + calculator.GetR());
            Console.WriteLine("N = " + calculator.GetN());
            Console.WriteLine("N_1 = " + calculator.GetN1());
            Console.WriteLine("tau = " + calculator.GetTau());
            Console.WriteLine("sigma_sq = " + calculator.GetSquaredSigma());
            Console.WriteLine("K = " + calculator.GetK());
            Console.WriteLine("M = " + calculator.GetM());
            Console.WriteLine("T = " + calculator.GetT());
            Console.WriteLine("h = " + calculator.GetH());
            Console.WriteLine("S0 Eps = " + calculator.GetS0Eps());
            Console.WriteLine();
        }
        
        private static AmericanOptionParameters GetParameters(bool saveSolutions, string workPath)
        {
            const double a = 0d;
            const double b = 1d;
            const double alpha = 1d;
            const double beta = 1d;
            const double r = 0.2d;
            const int n = 12000;
            const double tau = 1e-3d;
            const double sigmaSq = 0.2d;
            const double S0Eps = 1e-5d;
            const double K = 0.1d;
            const int M = 3650;
            return new AmericanOptionParameters(alpha, beta, a, b, n, r, tau, sigmaSq, K, S0Eps, /*h,*/ M, workPath, saveSolutions, 100d);
        }
    }
}