namespace PerpetualAmericanOptions
{
//    internal class JacobiCalculator
//    {
//        private readonly int n;
//        private readonly int iterCount;
//        private readonly double eps;
//
//        public JacobiCalculator(int n, int iterCount, double eps)
//        {
//            this.n = n;
//            this.iterCount = iterCount;
//            this.eps = eps;
//        }
//
//        internal void Calculate(double[] rp, double[] u_pr)
//        {
//            int iter = 0;
//            double maxErr;
//            double[] u = new double[u_pr.Length];
//            do {
//                double a_ii = (((2. * sigma1) / h1_sq) + ((2. * sigma2) / h2_sq) + (1. / tau));
//                double i_coef = -sigma1 / h1_sq;
//                double j_coef = -sigma2 / h2_sq;
//                for (int i = 1; i < n; ++i) {
//                        u[i] = (1d / a_ii) * (rp[i]
//                                              -
//                                              (i_coef * u_pr[i] + // up
//                                               i_coef * u_pr[i] + // bottom
//                                               j_coef * u_pr[i] + // left
//                                               j_coef * u_pr[i]) // right
//                            );
//                }
//
//                maxErr = double.MinValue;
//                for (int i = 0; i < n; ++i) {
//                        double val = Math.Abs(u[i] - u_pr[i]);
//                        if (val > maxErr) {
//                            maxErr = val;
//                        }
//                }
//
//                for (int i = 0; i < n; i++)
//                {
//                    u_pr[i] = u[i];
//                }
//                ++iter;
//            } while (maxErr > eps && iter < iterCount);
//        }
//    }
}