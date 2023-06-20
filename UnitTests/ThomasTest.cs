using System;
using System.Collections.Generic;
using CoreLib;
using NUnit.Framework;

namespace UnitTests
{
    public class ThomasTest
    {
        [Test]
        public void Test()
        {
            //https://pro-prof.com/forums/topic/sweep-method-for-solving-systems-of-linear-algebraic-equations
            int n = 3;
            // b - below main diagonal (indexed as [1;n-1])
            var a = new double[n];
            a[0] = 0d;
            a[1] = 5d;
            a[2] = 1d;

            // main diagonal of matrix (indexed as [0;n-1])
            var b = new double[n];
            b[0] = 2d;
            b[1] = 4d;
            b[2] = -3d;

            // d - up to main diagonal (indexed as [0;n-2])
            var c = new double[n];
            c[0] = -1d;
            c[1] = 2d;
            
            double[] f = {3, 6, 2};
            var u = SolveByTridiagonalMatrixAlgorithm(n, a, b, c, f);
            Utils.Print(u, "u_curr");
            Assert.AreEqual(u[0], 1.49d, 0.01d);
            Assert.AreEqual(u[1], -0.02d, 0.01d);
            Assert.AreEqual(u[2], -0.67d, 0.01d);
        }

        private static double[] SolveByTridiagonalMatrixAlgorithm(
            int n,
            IReadOnlyList<double> a,
            IReadOnlyList<double> b,
            IReadOnlyList<double> c,
            IReadOnlyList<double> d)
        {
            // for (var i = 0; i < c.Length; i++)
            //     if (Math.Abs(c[i]) < Math.Abs(b[i]) + Math.Abs(d[i]))
            //         throw new Exception(
            //             $"There is no diagonal dominance! i={i} {Math.Abs(b[i])} {Math.Abs(c[i])} {Math.Abs(d[i])} sum={Math.Abs(b[i] + d[i])} ");
            if (Math.Abs(b[0]) < double.Epsilon)
                throw new InvalidOperationException("c[0] == 0");
            var y = new double[n];
            var alpha = new double[n];
            var beta = new double[n];
            y[0] = b[0];
            alpha[0] = -c[0] / y[0];
            beta[0] = d[0] / y[0];
            Console.WriteLine(y[0]);
            Console.WriteLine(alpha[0]);
            Console.WriteLine(beta[0]);
            for (var i = 1; i < n - 1; ++i)
            {
                y[i] = b[i] + a[i] * alpha[i - 1];
                alpha[i] = -c[i] / y[i];
                beta[i] = (d[i] - a[i] * beta[i - 1]) / y[i];
                Console.WriteLine(y[i]);
                Console.WriteLine(alpha[i]);
                Console.WriteLine(beta[i]);
            }

            y[n - 1] = b[n - 1] + a[n - 1] * alpha[n - 2];
            beta[n - 1] = (d[n - 1] - a[n - 1] * beta[n - 2]) / y[n - 1];
            Console.WriteLine(y[n-1]);
            Console.WriteLine(alpha[n-1]);
            Console.WriteLine(beta[n-1]);

            Console.WriteLine();
            var x = new double[n];
            x[n - 1] = beta[n - 1];
            Console.WriteLine(x[n-1]);
            for (var i = n - 2; i >= 0; --i)
            {
                x[i] = alpha[i] * x[i + 1] + beta[i];
                Console.WriteLine(x[i]);
            }
            return x;
        }
    }
}