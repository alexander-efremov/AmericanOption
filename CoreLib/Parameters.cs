namespace CoreLib
{
    using System;

    public abstract class Parameters
    {
        protected Parameters(
            double alpha,
            double beta,
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
            if (alpha <= 0d)
            {
                throw new ArgumentException("Alpha must be greater than zero!");
            }
            
            this.Alpha = alpha;
            this.Beta = beta;
            this.A = a;
            this.A = a;
            this.B = b;
            this.N = n;
            this.N1 = n+1;
            this.R = r;
            this.Tau = tau;
            this.SigmaSq = sigmaSq;
            this.K = k;
            this.S0Eps = s0eps;
            this.WorkDir = workDir;
        }

        public double Alpha { get; }

        public double Beta { get; }
        public double A { get; }

        public double B { get; }

        public double K { get; }

        public int N { get; }
        
        public int N1 { get; }

        public double R { get; }

        public double S0Eps { get; }

        public double SigmaSq { get; }

        public double Tau { get; }

        public string WorkDir { get; private set; }
    }
}