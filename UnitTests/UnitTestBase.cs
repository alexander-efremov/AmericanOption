namespace UnitTests
{
    using System;
    using System.Collections.Generic;
    using CoreLib;
    using NUnit.Framework;

    public abstract class UnitTestBase
    {
        internal static double GetL1Error(AmericanOptionCalculatorBase cal, double[] exact, double[] calculated)
        {
            IEnumerable<double> err = Utils.GetError(exact, calculated);
            return Utils.GetL1(cal.GetH(), err);
        }

        internal static double GetL1Solution(AmericanOptionCalculatorBase cal, IEnumerable<double> calculatedV)
        {
            return Utils.GetL1(cal.GetH(), calculatedV);
        }

        protected virtual void PrintParameters(AmericanOptionCalculatorBase calculator)
        {
            Console.WriteLine("b = " + calculator.GetRightBoundary());
            Console.WriteLine("r = " + calculator.GetR());
            Console.WriteLine("N = " + calculator.GetN());
            Console.WriteLine("N_1 = " + calculator.GetN1());
            Console.WriteLine("tau = " + calculator.GetTau());
            Console.WriteLine("sigma_sq = " + calculator.GetSquaredSigma());
            Console.WriteLine("K = " + calculator.GetK());
        }

        [SetUp]
        protected void SetUp()
        {
        }
    }
}