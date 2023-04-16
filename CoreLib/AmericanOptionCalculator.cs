namespace CoreLib
{
    using System;

    public abstract class AmericanOptionCalculatorBase
    {
        private readonly double alpha;
        private readonly double a;
        private readonly double beta;
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
            alpha = parameters.Alpha;
            beta = parameters.Beta;
            a = parameters.A;
            b = parameters.B;
            sigmaSq = parameters.SigmaSq;
            tau = parameters.Tau;
            n = parameters.N;
            n1 = n + 1;
            r = parameters.R;
            K = parameters.K;
            h = parameters.B / parameters.N;
            S0Eps = parameters.S0Eps;
            workDir = parameters.WorkDir;

            CheckParameters();

            ThomasAlgorithmCalculator = new ThomasAlgorithmCalculator(n1);
        }

        protected ThomasAlgorithmCalculator ThomasAlgorithmCalculator { get; }

        public double Geta()
        {
            return a;
        }

        public double GetH()
        {
            return h;
        }

        public double GetK()
        {
            return K;
        }

        public double GetLeftBoundary()
        {
            return a;
        }

        public int GetN()
        {
            return n;
        }

        public int GetN1()
        {
            return n1;
        }

        public double GetR()
        {
            return r;
        }

        public double GetRightBoundary()
        {
            return b;
        }

        public double GetS0Eps()
        {
            return S0Eps;
        }

        public double GetSquaredSigma()
        {
            return sigmaSq;
        }

        public double GetTau()
        {
            return tau;
        }
        
        // public double GetAlpha(int i, int m)
        // {
        //     // alpha = alpha / m
        //     //return alpha;
        //     return alpha / (i  + m + 1);
        // }
        //
        // public double GetBeta(int m)
        // {
        //     // beta = beta_s - m*h^2
        //     return beta  ;
        //     //return beta - m * (h * h);
        // }

        protected void UpdateH1(double S0)
        {
            if (S0 < 0d)
            {
                throw new ArgumentException("S0");
            }

            h = (b - S0) / n;

            if (h <= 0d)
            {
                throw new ArgumentException("h");
            }
        }

        protected double Getb()
        {
            return b;
        }

        protected string GetWorkDir()
        {
            return workDir;
        }

        private void CheckParameters()
        {
            if (n < 2)
            {
                throw new ArgumentException("n");
            }

            if (n1 != n + 1)
            {
                throw new ArgumentException("n_1");
            }

            if (tau <= 0d)
            {
                throw new ArgumentException("tau");
            }

            if (K <= 0d)
            {
                throw new ArgumentException("K");
            }
        }
        
        public TecplotPrinterSpecial GetTecplotPrinter()
        {
            return new TecplotPrinterSpecial(0d, GetRightBoundary());
        }
    }
}