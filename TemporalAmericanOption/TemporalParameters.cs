using CoreLib;

namespace TemporalAmericanOption
{
    public class TemporalParameters : Parameters
    {
        public TemporalParameters(double a, double b, int n, double r, double tau, double sigma_sq, double k, double S0Eps, int M, double T) :
            base(a, b, n, r, tau, sigma_sq, k, S0Eps)
        {
            this.M = M;
            this.T = T;
        }

        public int M { get; }
        
        public double T { get; }
    }
}