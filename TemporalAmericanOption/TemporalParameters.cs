using CoreLib;

namespace TemporalAmericanOption
{
    public class TemporalParameters : Parameters
    {
        public TemporalParameters(double a, double b, int n, double r, double tau, double sigma, double k, double S0Eps, int M) :
            base(a, b, n, r, tau, sigma, k, S0Eps)
        {
            this.M = M;
        }

        public int M { get; }
    }
}