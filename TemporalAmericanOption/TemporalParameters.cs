using CoreLib;

namespace TemporalAmericanOption
{
    public class TemporalParameters : Parameters
    {
        public TemporalParameters(double a, double b, int n, double r, double tau, double sigma, double k, int M) :
            base(a, b, n, r, tau, sigma, k)
        {
            this.M = M;
        }

        public double M { get; }
    }
}