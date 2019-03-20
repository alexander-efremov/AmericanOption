// ReSharper disable CommentTypo

using System;

namespace CoreLib
{
    public class ThomasAlgorithmCalculator
    {
        private readonly int n;

        /// <summary>
        /// </summary>
        /// <param name="n">Pass n + 1</param>
        internal ThomasAlgorithmCalculator(int n)
        {
            this.n = n;
        }

        public double[] Calculate(double[] b, double[] c, double[] d, double[] r)
        {
            var delta = new double[n];
            var beta = new double[n];
            var lambda = new double[n];

            if (Math.Abs(c[0]) < double.Epsilon) throw new InvalidOperationException("c[0] == 0");

            delta[0] = c[0];
            beta[0] = -d[0] / delta[0];
            lambda[0] = r[0] / delta[0];

            for (var i = 1; i < n - 1; ++i)
            {
                delta[i] = c[i] + b[i] * beta[i - 1];
                beta[i] = -d[i] / delta[i];
                lambda[i] = (r[i] - b[i] * lambda[i - 1]) / delta[i];
            }

            var x = new double[n];
            x[n - 1] = lambda[n - 1];
            for (var i = n - 2; i >= 0; i--) x[i] = beta[i] * x[i + 1] + lambda[i];

            return x;
        }
    }
}