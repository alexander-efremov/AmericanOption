namespace TemporalAmericanOption
{
    using CoreLib;

    public class TemporalParameters : Parameters
    {
        public TemporalParameters(double a, double b, int n, double r, double tau, double sigma_sq, double k,
            double S0Eps, int M, double T, string workDir) :
            base(a, b, n, r, tau, sigma_sq, k, S0Eps, workDir)
        {
            this.M = M;
            this.T = T;
        }


        public bool SaveVSolutions { get; set; }

        public int M { get; }

        public double T { get; }
    }
}