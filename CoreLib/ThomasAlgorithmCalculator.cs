using System;
using System.Collections.Generic;

namespace CoreLib
{
    public class ThomasAlgorithmCalculator
    {
        private readonly int _n;

        /// <summary>
        ///     Initializes a new instance of the <see cref="ThomasAlgorithmCalculator" /> class.
        /// </summary>
        /// <param name="n">
        ///     Pass n + 1
        /// </param>
        internal ThomasAlgorithmCalculator(int n)
        {
            _n = n;
        }

        public double[] Calculate(double[] b, double[] c, double[] d, double[] r)
        {
            CheckDiagonalDominance(b, c, d);
            var delta = new double[_n];
            var beta = new double[_n];
            var lambda = new double[_n];

            if (Math.Abs(c[0]) < double.Epsilon) throw new InvalidOperationException("c[0] == 0");

            delta[0] = c[0];
            beta[0] = -d[0] / delta[0];
            lambda[0] = r[0] / delta[0];

            for (var i = 1; i < _n - 1; ++i)
            {
                delta[i] = c[i] + b[i] * beta[i - 1];
                beta[i] = -d[i] / delta[i];
                lambda[i] = (r[i] - b[i] * lambda[i - 1]) / delta[i];
            }

            var x = new double[_n];
            x[_n - 1] = lambda[_n - 1];
            for (var i = _n - 2; i >= 0; i--) x[i] = beta[i] * x[i + 1] + lambda[i];

            return x;
        }

        private void CheckDiagonalDominance(IReadOnlyList<double> lower, IReadOnlyList<double> central, IReadOnlyList<double> upper)
        {
            for (var i = 0; i < central.Count; i++)
                if (Math.Abs(central[i]) < Math.Abs(lower[i]) + Math.Abs(upper[i]))
                    throw new Exception(
                        $"There is no diagonal dominance! i={i} {Math.Abs(lower[i])} {Math.Abs(central[i])} {Math.Abs(upper[i])} sum={Math.Abs(lower[i] + upper[i])} ");
        }
    }
}