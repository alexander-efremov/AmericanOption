using CoreLib;

namespace AmericanOption
{
    public class AmericanOptionParameters : Parameters
    {
        public AmericanOptionParameters(double alpha, double beta, double a, double b, int n, double r, double tau, double sigma_sq, double k, double S0Eps, //double h,
            int m, string workDir, bool saveVSolutions, double smoothness) : base(alpha, beta, a, b, n, r, tau, sigma_sq, k, S0Eps, /*h,*/ workDir)
        {
            M = m;
            Smoothness = smoothness;
            SaveVSolutions = saveVSolutions;
        }

        public int M { get; }

        public bool SaveVSolutions { get; private set; }

        public double Smoothness { get; private set; }

        public double T => M * Tau;
    }
}