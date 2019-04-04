namespace CoreLib
{
    public abstract class Parameters
    {
        public double A { get; }
        public double B { get; }
        public int N { get; }
        public double R { get; }
        public double Tau { get; }
        public double Sigma { get; }
        public double K { get; }
        public double S0Eps { get; }

        public Parameters(double a, double b, int n, double r, double tau, double sigma, double k, double S0Eps)
        {
            A = a;
            B = b;
            N = n;
            R = r;
            Tau = tau;
            Sigma = sigma;
            K = k;
            this.S0Eps = S0Eps;
        }
    }
}