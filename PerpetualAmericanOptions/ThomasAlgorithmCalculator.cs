namespace PerpetualAmericanOptions
{
    internal class ThomasAlgorithmCalculator
    {
        private readonly int n_1;
        private readonly double h_sq;
        private readonly double tau;
        private readonly double sigma;

        internal ThomasAlgorithmCalculator(int n_1, double h_sq, double tau, double sigma)
        {
            this.n_1 = n_1;
            this.h_sq = h_sq;
            this.tau = tau;
            this.sigma = sigma;
        }

        /**
	 * n - число уравнений (строк матрицы)
	 * b - диагональ, лежащая под главной (нумеруется: [1;n-1])
	 * c - главная диагональ матрицы A (нумеруется: [0;n-1])
	 * d - диагональ, лежащая над главной (нумеруется: [0;n-2])
	 * f - правая часть (столбец)
	 * x - решение, массив x будет содержать ответ
	 */
        internal double[] GetB()
        {
            var b = new double[n_1];
            b[0] = b[1] = b[n_1 - 1] = 0.0;
            for (var i = 2; i < n_1 - 1; ++i)
            {
                b[i] = 1.0 / (8.0 * tau) - sigma / h_sq;
            }

            return b;
        }

        internal double[] GetC()
        {
            var c = new double[n_1];
            c[0] = c[n_1 - 1] = 0.0;
            for (var i = 1; i < n_1 - 1; ++i)
            {
                c[i] = 3.0 / (4.0 * tau) + 2.0 * sigma / h_sq;
            }

            return c;
        }

        internal double[] GetD()
        {
            var d = new double[n_1];
            d[0] = d[n_1 - 2] = d[n_1 - 1] = 0.0;
            for (var i = 1; i < n_1 - 2; ++i)
            {
                d[i] = 1.0 / (8.0 * tau) - sigma / h_sq;
            }

            return d;
        }
    }
}