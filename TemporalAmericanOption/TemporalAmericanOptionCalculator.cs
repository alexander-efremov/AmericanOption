using System;
using System.IO;
using CoreLib;
using TemporalAmericanOption;

// ReSharper disable CommentTypo

namespace PerpetualAmericanOptions
{
    // from new presentation with FEM
    public class TemporalAmericanOptionCalculator : AmericanOptionCalculator
    {
        private readonly int M;
        public TemporalAmericanOptionCalculator(TemporalParameters parameters) : base(parameters)
        {
            M = parameters.M;
        }

        public int GetM()
        {
            return M;
        }
        
        private double[] GetVT()
        {
            var arr = new double[GetN1()];
            for (var i = 0; i < arr.Length; ++i)
            {
                var S = GetLeftBoundary() + i * GetH();
                if (S < GetK())
                {
                    arr[i] = GetK() - S;
                }
                else
                {
                    arr[i] = 0d;
                }
            }

            return arr;
        }

        public double[] Solve()
        {
            var tecplotPrinter = new TecplotPrinterSpecial(GetN1(),
                0d,
                GetRightBoundary(),
                GetTau());
            var printer = new ThomasArrayPrinter();
            
            var VK1 = GetVT();
            var St = new double[GetM()];
            St[St.Length - 1] = GetK();

            Console.WriteLine("Time step = " + GetM() + string.Format(" h = {0} S(T) = {1}", GetH(), St[St.Length - 1]));
            Console.WriteLine("--------------------------------------------------");
            File.Delete(Environment.GetFolderPath(Environment.SpecialFolder.Desktop) + "\\" + "stat.txt");
            
            for (int k = GetM() - 1; k >= 1; --k)
            {
                Console.WriteLine("Time step = " + k);
                
                int iter = 0;
                double S0Old;
                double S0New = St[k];
                do
                {
                    iter++;
                    
                    S0Old = S0New;
                    UpdateH(S0Old);
                    var rp = CalculateRightPart(S0Old, VK1, GetH(), GetTau());
                    var b_t = GetB(GetN1(), S0Old, GetH(), GetSquaredSigma(), GetTau());
                    var c_t = GetC(GetN1(), S0Old, GetH(), GetSquaredSigma(), GetTau(), GetR());
                    var d_t = GetD(GetN1(), S0Old, GetH(), GetSquaredSigma(), GetTau());
                    var VK = ThomasAlgorithmCalculator.Calculate(b_t, c_t, d_t, rp);
                   
                    if (GetK() - VK[0] <= 0d) throw new Exception("GetK() - VK[0] <= 0d");
                    
                    S0New = GetK() - VK[0];
                    for (int i = 0; i < VK.Length; i++)
                    {
                        VK1[i] = VK[i];
                    }

                    using (var streamWriter =
                        File.AppendText(
                            Environment.GetFolderPath(Environment.SpecialFolder.Desktop) + "\\" + "stat.txt"))
                    {
                        streamWriter.WriteLine(
                            new string(' ', 2) + " Time step = {0} Iteration = " + iter +
                            " h = {1} S0 = {2} Abs(S0New-S0Old)={3} S0Eps={4} Cnd={5}", GetH(), S0Old,
                            Math.Abs(S0New - S0Old), GetS0Eps(), Math.Abs(S0New - S0Old) > GetS0Eps(), k);
                    }
                    Console.WriteLine(new string(' ', 2) + "Iteration = " + iter + " h = {0} S0_old = {1} S0_new = {2}", GetH(), S0Old, S0New);
                    
                } while (Math.Abs(S0New - S0Old) > GetS0Eps());

                St[k - 1] = S0New;
                Console.WriteLine("--------------------------------------------------");
            }

            return St;
        }

        private double[] CalculateRightPart(double S0, double[] VK1, double h, double tau)
        {
            var rp = new double[GetN1()];
            
            for (var i = 1; i < rp.Length - 1; ++i)
            {
                var si = S0 + i * h;
                var hmh = si - (S0 + (i - 1) * h); // h_{i-1/2}
                var hph = S0 + (i + 1) * h - si; // h_{i+1/2}
                var smh = si - 0.5d * hmh; // s_{i-1/2}
                var sph = si + 0.5d * hph; // s_{i+1/2}
                CheckHCorrectness(hph, tau, sph, GetR());
                CheckHCorrectness(hmh, tau, sph, GetR());
                var beta1 = 1d / (8d * tau) *
                            (1d + 2d * tau * GetR() * smh / hmh) *
                            (1d + 2d * tau * GetR() * smh / hmh);
                
                var beta2 = 1d / (8d * tau) *
                            (3d - 2d * tau * GetR() * smh / hmh) *
                            (1d + 2d * tau * GetR() * smh / hmh)
                            +
                            1d / (8d * tau) *
                            (3d + 2d * tau * GetR() * sph / hph) *
                            (1d - 2d * tau * GetR() * sph / hph);
                
                var beta3 = 1d / (8d * tau) *
                            (1d - 2d * tau * GetR() * sph / hph) *
                            (1d - 2d * tau * GetR() * sph / hph);
                
                var f = GetF(GetSquaredSigma(), GetR(), GetK(), i, S0);
                rp[i] = ((hmh + hph) / 2d) * f + ((hph + hmh) / 2d)*(beta1 * VK1[i - 1] + beta2 * VK1[i] + beta3 * VK1[i + 1]);
            }
            
            var si0 = S0 + 0 * h;
            var hph0 = S0 + (0 + 1) * h - si0; // h_{i+1/2}
            var sph0 = si0 + 0.5d * hph0; // s_{i+1/2}
            CheckHCorrectness(hph0, tau, sph0, GetR());
            var beta20 = 1d / (8d * tau) *
                         (3d + 2d * tau * GetR() * sph0 / hph0) *
                         (1d - 2d * tau * GetR() * sph0 / hph0);
            var beta30 = 1d / (8d * tau) *
                         (1d - (2d * tau * GetR() * sph0) / hph0) *
                         (1d - (2d * tau * GetR() * sph0) / hph0);
            var f0 = GetF(GetSquaredSigma(), GetR(), GetK(), 0, S0);
            
            rp[0] = hph0 * f0 + ((GetSquaredSigma() * si0 * si0) / 2d) + (hph0) * (beta20 * VK1[0] + beta30 * VK1[1]);
            
            rp[rp.Length - 1] = 0d;
            
            return rp;
        }

        private double GetF(double sigma_sq, double r, double K, int i, double S0)
        {
            return 0;
//            var si = S0 + i * GetH();
//            var p1 = sigma_sq / (2d * r);
//
//            var arg = K / (1 + sigma_sq / (2d * r));
//            var pow = (2d * r + sigma_sq) / sigma_sq;
//            var p2 = Math.Pow(arg, pow);
//
//            var arg2 = si;
//            var pow2 = -2d * r / sigma_sq;
//            var p3 = Math.Pow(arg2, pow2);
//
//            var v = p1 * p2 * p3;
//            return v;
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
        private double[] GetB(int n, double S0, double h, double sigma_sq, double tau)
        {
            var b = new double[n];
            for (var i = 1; i <= n - 1; ++i)
            {
                var si = S0 + i * h;
                var hmh = (S0 + i * h) - (S0 + (i - 1) * h); // h_{i - 1/2}
                if (hmh * hmh > 4d * tau * sigma_sq * si) throw new ArgumentException("hmh is invalid");

                b[i] = (hmh / (4d * tau)) - (sigma_sq * si * si) / (2d * hmh);
            }

            // right boundary cond
            b[n - 1] = 0d;

            return b;
        }

        // главная диагональ матрицы A (нумеруется: [0;n-1])
        // pass n = n + 1!
        private double[] GetC(int n, double S0, double h, double sigma_sq, double tau, double r)
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
                
                c[i] = (sigma_sq * si * si) / (2d * hmh) +
                       (sigma_sq * si * si) / (2d * hph) +
                       (hmh + hph) * (1d / (4d * tau) + r / 2d);
            }

            // left boundary condition
            var si0 = S0 + 0 * h;
            var hph0 = S0 + (0 + 1) * h - si0; // h_{i+1/2}
            var hmh0 = si0 - (S0 + (0d - 1d) * h); // h_{i-1/2}
            var sph0 = si0 + 0.5d * hph0; // s_{i+1/2}
            var smh0 = si0 - 0.5d * hmh0; // s_{i-1/2}
            CheckHCorrectness(hph0, tau, sph0, r);
            CheckHCorrectness(hmh0, tau, smh0, r);
            
            c[0] = (sigma_sq * si0 * si0) / (2d * hph0) + (hph0 * r) / 2d + hph0 / (4d * tau);

            // right boundary condition
            c[n - 1] = 0d;

            return c;
        }

        // диагональ, лежащая над главной (нумеруется: [0;n-2])
        // pass n = n + 1!
        private double[] GetD(int n, double S0, double h, double sigma_sq, double tau)
        {
            var d = new double[n];
            for (var i = 0; i <= n - 2; ++i)
            {
                var si = S0 + i * h;
                var hph = S0 + (i + 1) * h - si; // h_{i+1/2}
               
                d[i] = hph / (4d * tau) - (sigma_sq * si * si) / (2d * hph);
            }

            // left boundary condition
            var si0 = S0 + 0 * h;
            var hph0 = S0 + (0 + 1) * h - si0; // h_{i+1/2}
           
            d[0] = hph0 / (4d * tau) - (sigma_sq * si0 * si0) / (2d * hph0);

            // right boundary condition
            d[n - 2] = 0d;
            return d;
        }
        
        private static void CheckHCorrectness(double h, double tau, double sph, double r)
        {
            if (tau / h > 1d / (2d * r * sph))
            {
                throw new ArgumentOutOfRangeException("tau/h");
            }
        }
    }
}