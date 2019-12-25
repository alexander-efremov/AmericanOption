namespace PerpetualAmericanOptions
{
    using System;
    using CoreLib;
    using NUnit.Framework;

    public abstract class UnitTestBase
    {
        protected abstract string GetWorkingDir();

        [SetUp]
        protected void SetUp()
        {
        }

        internal double GetL1Error(AmericanOptionCalculator cal, double[] exact, double[] calculated)
        {
            double[] err = Utils.GetError(exact, calculated);
            return Utils.GetL1(cal.GetH(), err);
        }

        internal double GetL1Solution(AmericanOptionCalculator cal, double[] calculatedV)
        {
            return Utils.GetL1(cal.GetH(), calculatedV);
        }

        protected virtual void PrintParameters(AmericanOptionCalculator calculator)
        {
            Console.WriteLine("b = " + calculator.GetRightBoundary());
            Console.WriteLine("r = " + calculator.GetR());
            Console.WriteLine("N = " + calculator.GetN());
            Console.WriteLine("N_1 = " + calculator.GetN1());
            Console.WriteLine("tau = " + calculator.GetTau());
            Console.WriteLine("sigma_sq = " + calculator.GetSquaredSigma());
            Console.WriteLine("K = " + calculator.GetK());
        }
    }
}