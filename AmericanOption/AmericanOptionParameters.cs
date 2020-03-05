namespace AmericanOption
{
    using CoreLib;

    public class AmericanOptionParameters : Parameters
    {
        public AmericanOptionParameters(double a, double b, int n, double r, double tau, double sigma_sq, double k,
            double S0Eps, int M, string workDir, bool saveVSolutions, double smoothness)
            : base(a, b, n, r, tau, sigma_sq, k, S0Eps, workDir)
        {
            this.M = M;
            this.Smoothness = smoothness;
            this.SaveVSolutions = saveVSolutions;
        }
        
        public double Smoothness { get; private set; }

        public bool SaveVSolutions { get; private set; }

        public int M { get; }

        public double T
        {
            get { return M * Tau; }
        }
    }
}