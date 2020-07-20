namespace PerpetualAmericanOptions
{
    using CoreLib;

    public class PerpetualParameters : Parameters
    {
        public PerpetualParameters(
            double alpha,
            double beta,
            double a,
            double b,
            int n,
            double r,
            double tau,
            double sigmaSq,
            double k,
            double S0Eps,
            double h,
            string workDir)
            : base(alpha, beta, a, b, n, r, tau, sigmaSq, k, S0Eps, h, workDir)
        {
        }
    }
}