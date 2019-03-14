using System;
using System.Configuration;
using NUnit.Framework;

// ReSharper disable CommentTypo

namespace PerpetualAmericanOptions
{
    // from new presentation with FEM
    internal class PerpetualAmericanOptionCalculator1
    {
        private readonly ThomasAlgorithmCalculator _thomasAlgorithmCalculator;
        private readonly int n;
        private readonly int n_1;
        private readonly double b;
        private readonly double sigma;
        private readonly double sigma_sq;
        private readonly double tau;
        private readonly double r;
        private double h;
        private double h_sq;
        private readonly double K;
        
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
        
        public double GetSquaredH()
        {
            return h_sq;
        }
        
        public double GetK()
        {
            return K;
        }

        public PerpetualAmericanOptionCalculator1(Parameters parameters)
        {
            b = parameters.B;
            sigma = parameters.Sigma;
            sigma_sq = sigma * sigma;
            tau = parameters.Tau;
            n = parameters.N;
            n_1 = n + 1;
            r = parameters.R;
            K = parameters.K;
            h = b / n;
            h_sq = h * h;

            CheckParameters();

            _thomasAlgorithmCalculator = new ThomasAlgorithmCalculator(n_1);
        }

        private double[] CalculateRightPart(double S0)
        {
            double[] V = new double[n_1 + 1]; // add + 1 to calculate right boundary
            for (int i = 0; i <= n_1; i++)
            {
                var si = S0 + i * h;
                V[i] = GetV(sigma_sq, r, K, si);
            }

//            return V;
            var rp = new double[n_1];
            var si0 = S0 + 0 * h;
            var hph0 = (S0 + (0 + 1) * h) - si0; // h_{i+1/2}
            var sph0 = si0 + 0.5d * hph0; // s_{i+1/2}
            check_h(hph0, tau, sph0, r);
            var beta20 = (1d / (8d * tau)) * 
                        (3d - (2d * tau * r * sph0) / h) *
                        (1d + (2d * tau * r * sph0) / h);
            var beta30 = (1d / (8d * tau)) * 
                        (1d - (2d * tau * r * sph0) / h) *
                        (1d - (2d * tau * r * sph0) / h);
            var val0 = beta20 * V[0] + beta30 * V[1];
            val0 = val0 * (1d / tau);
            rp[0] = val0 + (sigma_sq*si0*si0)/2d;
            for (int i = 1; i <= n_1 - 1; i++)
            {
                var si = S0 + i * h;
                var hmh = si - (S0 + (i - 1) * h); // h_{i-1/2}
                var hph = (S0 + (i + 1) * h) - si; // h_{i+1/2}
                var smh = si - 0.5d * hmh; // s_{i-1/2}
                var sph = si + 0.5d * hph; // s_{i+1/2}
                check_h(hph, tau, sph, r);
                check_h(hmh, tau, sph, r);
                var beta1 = (1d / (8d * tau)) * 
                            (1d + (2d * tau * r * smh) / h) * 
                            (1d + (2d * tau * r * smh) / h);
                var beta2 = (1d / (8d * tau)) * 
                               (3d + (2d * tau * r * smh) / h) *
                               (1d - (2d * tau * r * smh) / h)
                               + 
                               (1d / (8d * tau)) * 
                               (3d - (2d * tau * r * sph) / h) *
                               (1d + (2d * tau * r * sph) / h);
                var beta3 = (1d / (8d * tau)) * 
                               (1d - (2d * tau * r * sph) / h) *
                               (1d - (2d * tau * r * sph) / h);
                var val = beta1 * V[i - 1] + beta2 * V[i] + beta3 * V[i + 1];
                val = val * (1d / tau);

                rp[i] = GetF(i) + val;
            }

            return rp;
        }

        private double GetF(int i)
        {
            return 0;
        }

        private static void check_h(double h, double tau, double sph, double r)
        {
            Assert.LessOrEqual(tau/h, 1d/(2d*r*sph));
        }

        public Tuple<double[], double> Solve()
        {
            var tecplotPrinter = new TecplotPrinter(GetN1(),
                GetH(), 
                0d,
                GetRightBoundary(),
                GetTau());
            var printer = new ThomasArrayPrinter();
            double[] calculatedV = new double[n_1];
            
//            double S0 = 0.5d * ((b - 0d) / n);
            var S0 = GetExactS0();
            //while (Math.Abs(GetExactS0() - S0) > 10e-5)
            {
                CalculateH(S0);
                CheckParameters2();
                var b_t = GetB(n_1, S0, h, sigma_sq, tau);
                var c_t = GetC(n_1, S0, h, sigma_sq, tau, r);
                var d_t = GetD(n_1, S0, h, sigma_sq, tau);
                printer.PrintThomasArrays(b_t, c_t, d_t);
                var rp = CalculateRightPart(S0);
                tecplotPrinter.PrintXY(Environment.GetFolderPath(Environment.SpecialFolder.Desktop) + "\\" + "rp", 0d, rp, S0);
                calculatedV = _thomasAlgorithmCalculator.Calculate(b_t, c_t, d_t, rp);

                // for better vis
                calculatedV[calculatedV.Length - 1] = calculatedV[calculatedV.Length - 2]; 
                tecplotPrinter.PrintXY(Environment.GetFolderPath(Environment.SpecialFolder.Desktop) + "\\" + "v", 0d,  calculatedV, S0);
                tecplotPrinter.PrintXY(Environment.GetFolderPath(Environment.SpecialFolder.Desktop) + "\\" + "rp_v", 0d, rp, calculatedV, S0);
//                var vs0 = GetVKS();
//                tecplotPrinter.PrintXY(Environment.GetFolderPath(Environment.SpecialFolder.Desktop) + "\\" + "v", 0d, vs0, calculatedV);
//
//                var err = Utils.FillArrayDiff(GetExactSolution(GetExactS0()), rp);
//                tecplotPrinter.PrintXY("ex_minus_rp", 0d, err);
            }

            return Tuple.Create(calculatedV, S0);
        }

        public double[] GetExactSolution()
        {
            var arr = new double[n_1];
            for (var i = 0; i < arr.Length; ++i)
            {
                double si = 0 + i * h;
                arr[i] = GetV(sigma_sq, r, K, si);
            }
            // V(S) does not touch 0
            // adjust value for better visualization
            arr[0] = arr[1];
            return arr;
        }
        
        public double GetExactS0()
        {
            return K / (1d + (sigma_sq / (2d * r)));
        }

        public double[] GetVKS()
        {
            var res = new double[n_1];
            for (int i = 0; i < res.Length; i++)
            {
                res[i] = K - (i * h);
                if (res[i] < 0d)
                {
                    res[i] = 0d;
                }
            }

            return res;
        }
        
        private void CalculateH(double S0)
        {
            h = (b - S0) / n;
            h_sq = h * h;
        }
        
        private static double GetV(double sigma_sq, double r, double K, double si)
        {
            var p1 = sigma_sq / (2d * r);
            
            var arg = ( K / (1 + (sigma_sq / (2d * r))) );
            var pow = (2d * r + sigma_sq) / sigma_sq;
            var p2 = Math.Pow(arg, pow);
            
            var arg2 = si;
            var pow2 = (-2d * r) / sigma_sq;
            var p3 = Math.Pow(arg2, pow2);

            var v = p1 * p2 * p3;
            return v;
        }
        
        /**
         * n - число уравнений (строк матрицы)
         * b - диагональ, лежащая под главной (нумеруется: [1;n-1])
         * c - главная диагональ матрицы A (нумеруется: [0;n-1])
         * d - диагональ, лежащая над главной (нумеруется: [0;n-2])
         * f - правая часть (столбец)
         * x - решение, массив x будет содержать ответ
         */
        // pass n = n + 1!
        private static double[] GetB(int n, double S0, double h, double sigma_sq, double tau)
        {
            var b = new double[n];
            for (var i = 1; i <= n - 1; ++i)
            {
                var si = S0 + i * h;
                var hmh = si - (S0 + (i - 1) * h); // h_{i - 1/2}
                if (hmh * hmh > 4d * tau * sigma_sq * si)
                {
                    throw new ArgumentException("hmh is invalid");
                }

                b[i] = (hmh / (4d * tau)) - (sigma_sq * si * si) / (2d * hmh);
            }
            
            // right boundary cond
            //b[n - 1] = 0d;

            return b;
        }

        // главная диагональ матрицы A (нумеруется: [0;n-1])
        // pass n = n + 1!
        private static double[] GetC(int n, double S0, double h, double sigma_sq, double tau, double r)
        {
            var c = new double[n];
            for (var i = 1; i <= n - 1; ++i)
            {
                var si = S0 + i * h;
                var hmh = si - (S0 + (i - 1) * h); // h_{i-1/2}
                var hph = (S0 + (i + 1) * h) - si; // h_{i+1/2}
                if (hmh * hmh > 4d * tau * sigma_sq * si)
                {
                    throw new ArgumentException("hmh is invalid");
                }

                c[i] = ((sigma_sq * si * si) / (2d * hmh)) +
                       ((sigma_sq * si * si) / (2d * hph)) +
                       (hmh + hph) *
                       ((1d / (4d * tau)) + r / 2d); 
            }
            
            // left boundary condition
            var si0 = S0 + 0 * h;
            var hph0 = (S0 + (0 + 1) * h) - si0; // h_{i+1/2}
            c[0] = (sigma_sq * si0 * si0) / (2d * hph0) + (hph0 / 2d) * r + hph0/(4d*tau);

            // right boundary condition
            // c[n - 1] = 0d;

            return c;
        }

        // диагональ, лежащая над главной (нумеруется: [0;n-2])
        // pass n = n + 1!
        private static double[] GetD(int n, double S0, double h, double sigma_sq, double tau)
                 {
            var d = new double[n];
            for (var i = 0; i <= n - 2; ++i)
            {
                var si = S0 + i * h;
                var hph = (S0 + (i + 1) * h) - si; // h_{i+1/2}
                d[i] = (hph / (4d * tau)) - (sigma_sq * si * si) / (2d * hph);
            }

            // left boundary condition
            var si0 = S0 + 0 * h;
            var hph0 = (S0 + (0 + 1) * h) - si0; // h_{i+1/2}
            d[0] = -(sigma_sq * si0 * si0) / (2d * hph0) + hph0 / (4d*tau);
            
            // right boundary condition
            d[n - 2] = 0d;
            return d;
        }

        private void CheckParameters()
        {
            Assert.True(h > 0.0);
            Assert.AreEqual(h * h, h_sq);
            Assert.True(n > 1);
            Assert.AreEqual(n + 1, n_1);
            Assert.True(tau > 0d);
            Assert.True(K > 0d); 
        }
        
        private void CheckParameters2()
        {
            Assert.True(h > 0.0);
            Assert.AreEqual(h * h, h_sq);
            
//            // let's check scheme t/h restrictions
//            // left for Thomas algo, right for postiive defined M matrix
//            // (r * h) / (4*sigmaSq*Si) <= tau/h <= S_i+1/2 / (2 * r)
//            // Si = a + i*h S_i+1/2 = a + i*h + h/2
//            // then 
//            // (r * h) / (4*sigmaSq*(a+i*h)) <= tau/h <= (a+i*h + 0.5h) / (2 * r)
//            // min(Si) = h (except i = 0) and min(Si+1/2) = 0.5h
//            // then
//            // (r * h) / (4*sigmaSq*h) <= tau/h <= 0.5h / (2 * r)
//            // 0.25r / sigmaSq <= tau/h <= 0.25h / r
//
//            Assert.GreaterOrEqual(tau / h, (0.25d * r) / sigma_sq);
//            Assert.GreaterOrEqual((0.25d * h) / r, tau / h);
        }
    } 
}