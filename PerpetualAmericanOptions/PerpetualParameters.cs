namespace PerpetualAmericanOptions
{
    using CoreLib;

    public class PerpetualParameters : Parameters
    {
        public PerpetualParameters(
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
            : base(a, b, n, r, tau, sigmaSq, k, S0Eps, h, workDir)
        {
        }
    }
}