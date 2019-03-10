// ReSharper disable CommentTypo

namespace PerpetualAmericanOptions
{
    internal class ThomasAlgorithmCalculator
    {
        private readonly int n;

        /// <summary>
        /// </summary>
        /// <param name="n">Pass n + 1</param>
        internal ThomasAlgorithmCalculator(int n)
        {
            this.n = n;
        }
        
        internal double[] Calculate(double[] b, double[] c, double[] d, double[] r) {
            double[] delta = new double[n];
            double[] beta = new double[n];
            double[] lambda = new double[n];

            delta[1] = c[1];
            beta[1] = -d[1] / delta[1];
            lambda[1] = r[1] / delta[1];

            for (int i = 2; i < n - 1; ++i) {
                delta[i] = c[i] + b[i] * beta[i - 1];
                beta[i] = -d[i] / delta[i];
                lambda[i] = (r[i] - b[i] * lambda[i - 1]) / delta[i];
            }

            var x = new double[n];

            x[n - 2] = lambda[n - 2];

            for (int i = n - 3; i >= 1; i--)
            {
                x[i] = beta[i] * x[i + 1] + lambda[i];
            }

            return x;
        }
    }
}