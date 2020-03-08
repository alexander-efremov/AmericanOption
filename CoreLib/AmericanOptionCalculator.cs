namespace CoreLib
{
    using System;

    public abstract class AmericanOptionCalculatorBase
    {
        private readonly double a;
        private readonly double b;
        private readonly double K;
        private readonly int n;
        private readonly int n1;
        private readonly double r;
        private readonly double S0Eps;
        private readonly double sigmaSq;
        private readonly double tau;
        private readonly string workDir;
        private double h;

        protected AmericanOptionCalculatorBase(Parameters parameters)
        {
            this.a = parameters.A;
            this.b = parameters.B;
            this.sigmaSq = parameters.SigmaSq;
            this.tau = parameters.Tau;
            this.n = parameters.N;
            this.n1 = this.n + 1;
            this.r = parameters.R;
            this.K = parameters.K;
            this.h = this.b / this.n;
            this.S0Eps = parameters.S0Eps;
            this.workDir = parameters.WorkDir;

            this.CheckParameters();

            this.ThomasAlgorithmCalculator = new ThomasAlgorithmCalculator(this.n1);
        }

        protected ThomasAlgorithmCalculator ThomasAlgorithmCalculator { get; }

        public double Geta()
        {
            return this.a;
        }

        public double GetH()
        {
            return this.h;
        }

        public double GetK()
        {
            return this.K;
        }

        public double GetLeftBoundary()
        {
            return this.a;
        }

        public int GetN()
        {
            return this.n;
        }

        public int GetN1()
        {
            return this.n1;
        }

        public double GetR()
        {
            return this.r;
        }

        public double GetRightBoundary()
        {
            return this.b;
        }

        public double GetS0Eps()
        {
            return this.S0Eps;
        }

        public double GetSquaredSigma()
        {
            return this.sigmaSq;
        }

        public double GetTau()
        {
            return this.tau;
        }

        public void UpdateH(double S0)
        {
            if (S0 < 0d)
            {
                throw new ArgumentException("S0");
            }

            this.h = (this.b - S0) / this.n;

            if (this.h <= 0d)
            {
                throw new ArgumentException("h");
            }
        }

        protected double Getb()
        {
            return this.b;
        }

        protected string GetWorkDir()
        {
            return this.workDir;
        }

        private void CheckParameters()
        {
            if (this.n < 2)
            {
                throw new ArgumentException("n");
            }

            if (this.n1 != this.n + 1)
            {
                throw new ArgumentException("n_1");
            }

            if (this.tau <= 0d)
            {
                throw new ArgumentException("tau");
            }

            if (this.K <= 0d)
            {
                throw new ArgumentException("K");
            }
        }
    }
}