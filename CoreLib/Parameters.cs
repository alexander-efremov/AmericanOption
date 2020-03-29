namespace CoreLib
{
    public abstract class Parameters
    {
        protected Parameters(
            double a,
            double b,
            int n,
            double r,
            double tau,
            double sigmaSq,
            double k,
            double s0eps,
            double h,
            string workDir)
        {
            this.A = a;
            this.B = b;
            this.N = n;
            this.R = r;
            this.Tau = tau;
            this.SigmaSq = sigmaSq;
            this.K = k;
            this.S0Eps = s0eps;
            this.WorkDir = workDir;
        }

        public double A { get; }

        public double B { get; }

        public double K { get; }

        public int N { get; }

        public double R { get; }

        public double S0Eps { get; }

        public double SigmaSq { get; }

        public double Tau { get; }

        public string WorkDir { get; private set; }
    }
}