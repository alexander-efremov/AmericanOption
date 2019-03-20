using CoreLib;

namespace PerpetualAmericanOptions
{
    

    public class PerpetualParameters : Parameters
    {
        public PerpetualParameters(double a, double b, int n, double r, double tau, double sigma, double k) : base(a, b,
            n, r, tau, sigma, k)
        {
        }
    }
}