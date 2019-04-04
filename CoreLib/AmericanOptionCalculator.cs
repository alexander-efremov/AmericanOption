using System;

namespace CoreLib
{
    public abstract class AmericanOptionCalculator
    {
        private readonly ThomasAlgorithmCalculator _thomasAlgorithmCalculator;
        private readonly double a;
        private readonly double b;
        private readonly double K;
        private readonly int n;
        private readonly int n_1;
        private readonly double r;
        private readonly double sigma;
        private readonly double sigma_sq;
        private readonly double tau;
        private readonly double S0Eps;
        private double h;
        
        public double GetS0Eps()
        {
            return S0Eps;
        }

        public int GetN()
        {
            return n;
        }

        public int GetN1()
        {
            return n_1;
        }

        public double GetRightBoundary()
        {
            return b;
        }
        
        public double GetLeftBoundary()
        {
            return a;
        }

        public double GetSigma()
        {
            return sigma;
        }

        public double GetSquaredSigma()
        {
            return sigma_sq;
        }

        public double GetTau()
        {
            return tau;
        }

        public double GetR()
        {
            return r;
        }

        public double GetH()
        {
            return h;
        }

        public double GetK()
        {
            return K;
        }

        public ThomasAlgorithmCalculator ThomasAlgorithmCalculator
        {
            get { return _thomasAlgorithmCalculator; }
        }
        
        public AmericanOptionCalculator(Parameters parameters)
        {
            a = parameters.A;
            b = parameters.B;
            sigma = parameters.Sigma;
            sigma_sq = sigma * sigma;
            tau = parameters.Tau;
            n = parameters.N;
            n_1 = n + 1;
            r = parameters.R;
            K = parameters.K;
            h = b / n;
            S0Eps = parameters.S0Eps;

            CheckParameters();

            _thomasAlgorithmCalculator = new ThomasAlgorithmCalculator(n_1);
        }
        
        protected void UpdateH(double S0)
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

        private void CheckParameters()
        {
            if (n < 2) throw new ArgumentException("n");
            if (n_1 != n + 1) throw new ArgumentException("n_1");
            if (tau <= 0d) throw new ArgumentException("tau");
            if (K <= 0d) throw new ArgumentException("K");
        }
    }
}