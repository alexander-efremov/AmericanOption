namespace AmericanOption
{
    using CoreLib;

    public class AmericanOptionParameters : Parameters
    {
        public AmericanOptionParameters(double a, double b, int n, double r, double tau, double sigma_sq, double k,
            double S0Eps, int M, string workDir, bool saveVSolutions)
            : base(a, b, n, r, tau, sigma_sq, k, S0Eps, workDir)
        {
            this.M = M;
            this.SaveVSolutions = saveVSolutions;
        }

        public bool SaveVSolutions { get; private set; }

        public int M { get; }

        public double T
        {
            get { return M * Tau; }
        }
    }
}