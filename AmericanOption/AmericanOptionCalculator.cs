using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.IO;
using System.Linq;
using CoreLib;

namespace AmericanOption
{
    // from new presentation with FEM
    [SuppressMessage("ReSharper", "CommentTypo", Justification = "OK.")]
    public class AmericanOptionCalculator : AmericanOptionCalculatorBase
    {
        private readonly bool allowOutputConsole;

        private readonly bool allowOutputFile;

        private readonly int MValue;

        private readonly string outputPath;

        private readonly string outputPathRp;

        private readonly string outputPathStat;

        private readonly bool saveSolutions;

        private readonly double smoothness;

        private readonly List<SolutionData> solutions = new List<SolutionData>();

        private readonly double TValue;

        public AmericanOptionCalculator(AmericanOptionParameters parameters, bool allowOutputFile, bool allowOutputConsole)
            : base(parameters)
        {
            this.allowOutputFile = allowOutputFile;
            this.allowOutputConsole = allowOutputConsole;
            outputPath = parameters.WorkDir;
            if (string.IsNullOrEmpty(outputPath))
                outputPath = GetWorkDir() + "AO/";

            MValue = parameters.M;
            TValue = parameters.T;
            smoothness = parameters.Smoothness;
            saveSolutions = parameters.SaveVSolutions;

            outputPathStat = Path.Combine(outputPath, "stat");
            if (!Directory.Exists(outputPathStat))
                Directory.CreateDirectory(outputPathStat);

            outputPathRp = Path.Combine(outputPath, "rp");
            if (!Directory.Exists(outputPathRp))
                Directory.CreateDirectory(outputPathRp);
        }

        /// <summary>
        /// Returns exact solutions V(S(t)) in reversed order (from the M to the 0 time step)
        /// </summary>
        /// <param name="S0t">S0 in direct order from 0 to M time layer</param>
        /// <returns></returns>
        public List<SolutionData> GetExactSolutions(double[] S0t)
        {
            if (Math.Abs(GetM() * GetTau() - GetT()) > double.Epsilon)
                throw new Exception("GetExactSolutions");

            var list = new List<SolutionData>
            {
                new SolutionData(S0t, GetM())
                {
                    Solution = GetVonT() // V(S, T) = (K - S)+
                }
            };

            for (var k = GetM() - 1; k >= 0; k--)
            {
                Point[] solution = GetExactSolution(GetTau(), k, S0t[k]);
                var solutionData = new SolutionData(S0t, k)
                {
                    Solution = solution
                };
                list.Add(solutionData);
            }

            return list;
        }

        /// <summary>
        /// Returns exact solutions V(S(t)) in reversed order (from the M to the 0 time step)
        /// </summary>
        /// <returns></returns>
        public List<SolutionData> GetExactSolutions2(double[] S0t)
        {
            if (Math.Abs(GetM() * GetTau() - GetT()) > double.Epsilon)
                throw new Exception("GetExactSolutions");

            int startPosOfS01;
            Point[] solution2 = GetExactSolution2(GetM(), GetM(), 0d, out startPosOfS01);
            var item = new SolutionData(0d, GetM())
            {
                Solution = solution2,
                StartPosOfS0 = startPosOfS01
            };
            var list = new List<SolutionData>
            {
                item
            };
            
            for (var k = GetM() - 1; k >= 0; k--)
            {
                int startPosOfS02;
                Point[] solution = GetExactSolution2(GetTau(), k, S0t[k], out  startPosOfS02);
                var solutionData = new SolutionData(S0t, k)
                {
                    Solution = solution,
                    StartPosOfS0 = startPosOfS02
                };
                list.Add(solutionData);
            }

            return list;
        }

        public int GetM()
        {
            return MValue;
        }

        /// <summary>
        /// Returns numeric solutions V(S(t)) in reversed order (from the M to the 0 time step)
        /// </summary>
        /// <returns></returns>
        public List<SolutionData> GetNumericSolutions()
        {
            // solutions were filled in reversed order from the M to 0
            // then here we reverse it
            var list = new List<SolutionData>();

            foreach (var data in solutions)
            {
                Point[] array = data.Solution.ToArray();
                var solutionData = new SolutionData(data.S0, data.k);
                solutionData.Solution = array;
                list.Add(solutionData);
            }
            
            return list;
        }

        public double GetT()
        {
            return TValue;
        }

        /// <summary>
        /// Solves the american option problem
        /// </summary>
        /// <returns>Array of S0 values in REVERSED order (S0t[0] refers to the M time step, S0t[0] refers to the 0 time step)</returns>
        public double[] Solve()
        {
            // var printer = new ThomasArrayPrinter();
            double[] VNext = GetVST(); // V(S, T) = (K - S)+

            var printer = GetTecplotPrinter();
            // tecplotPrinter.PrintXY(Path.Combine(_outputPath, "VST"), 0d, GetH(), V);
            
            var S0t = new double[GetM() + 1];
            S0t[GetM()] = GetK();
            
            if (saveSolutions)
                SaveNumericSolution(VNext, GetM(), S0t[GetM()]);

            PrintHeader(S0t);
            CreateFileWithStatistics();
            PrintVST(printer, VNext, S0t);

            for (var m = GetM() - 1; m >= 0; --m)
            {
                var tau = GetTau();
                if (allowOutputConsole)
                    Console.WriteLine("Time step = " + m);

                double[] VCurrent;
                var iter = 0;
                var S0Next = S0t[m + 1];
                var S0Current = S0Next;
                double S0Diff;
                if (allowOutputConsole)
                {
                    Console.WriteLine("m = " + m);
                    Console.WriteLine("S0Current = " + S0Current);
                }
                
                //var h = this.GetH();
                var h = (Getb() - GetK()) / GetN();
                
                do
                {
                    // calculate new V(S, t_k)
                    VCurrent = CalculateV(S0Current, VNext, tau, m, h, out double[] rp);

                    if (VCurrent[0] <= 0d)
                        throw new Exception("VCurrent[0] <= 0d");

                    // calculate new S0(t)_i
                    var S0CurrentPrev = S0Current;
                    S0Current = GetK() - VCurrent[0];
                    CheckS0(GetK() - VCurrent[0], S0CurrentPrev);
                    PrintData(S0Current, printer, rp, ++iter, m, S0CurrentPrev, VCurrent);
                    S0Diff = Math.Abs(GetK() - VCurrent[0] - S0Current);
                } while (S0Diff > GetS0Eps());

                // Console.WriteLine(m + " " + S0Current);

                // TODO: why we get some negative values instead of zero?
                RemoveNegativeValues(VCurrent);
                
                if (saveSolutions)
                    SaveNumericSolution(VCurrent, m, S0Current);

                S0t[m] = S0Current;
                for (var i = 0; i < VCurrent.Length; i++)
                    VNext[i] = VCurrent[i];

                if (allowOutputFile)
                    printer.PrintXY(Path.Combine(outputPath, "V" + m), GetTau(), m, GetH(), VCurrent, S0Current);

                if (allowOutputConsole)
                    Console.WriteLine("--------------------------------------------------");
            }

            return S0t;
        }

        private static void CheckHCorrectness(double h, double tau, double sph, double r)
        {
            if (tau / h > 1d / (2d * r * sph))
            {
                // ReSharper disable once NotResolvedInText
                //throw new ArgumentOutOfRangeException("tau/h");
            }
        }

        private static void RemoveNegativeValues(IList<double> VCurrent)
        {
            for (var i = 0; i < VCurrent.Count; i++)
                if (VCurrent[i] < 0d)
                    VCurrent[i] = 0d;
        }

        private double[] CalculateRightPart(double S0, IReadOnlyList<double> Vk1, double h, double tau, int m)
        {
            var rp = new double[GetN1()];
            var r = GetR();
            var sigma = GetSquaredSigma();

            // var hph0 = Math.Pow(this.GetAlpha(0, m) * (0+1) * h, this.GetBeta(m)); // h_{i+1/2}
            var hph0 = h;
            var hmh0 = hph0; // h_{i-1/2}
            var smh0 = S0 - 0.5d * hmh0; // s_{i-1/2}
            var sph0 = S0 + 0.5d * hph0; // s_{i+1/2}
            CheckHCorrectness(hmh0, tau, sph0, r);
            CheckHCorrectness(hph0, tau, sph0, r);
            // var betam1 = 1d / (8d * tau) * (1d + (2d * tau * r * smh0) / hmh0) * (1d + (2d * tau * r * smh0) / hmh0);
            var beta0 =  1d / (8d * tau) * (3d - 2d * tau * r * smh0 / hmh0) * (1d + 2d * tau * r * smh0 / hmh0) 
                + 1d / (8d * tau) * (3d + 2d * tau * r * sph0 / hph0) * (1d - 2d * tau * r * sph0 / hph0);
            var betap1 = 1d / (8d * tau) * (1d - 2d * tau * r * sph0 / hph0) * (1d - 2d * tau * r * sph0 / hph0);
            var f0 = GetF(sigma, r, GetK(), h, S0, GetTau(), m, GetT(), S0, 0);

            var x0 = -hph0 / 2d * f0 + sigma * S0 * S0 / 2d
                                     + hph0 / 2d * (beta0 * Vk1[0] + betap1 * Vk1[1]);
            rp[0] = x0;

            for (var i = 1; i < rp.Length - 1; ++i)
            {
                //Console.WriteLine("i={0} beta={1} alpha={2} h={3}", i, GetBeta(m), GetAlpha(m), Math.Pow(GetAlpha(m) * i * h, GetBeta(m)));
                var xi = S0 + i * h;
                // var hmh = xi - (S0 + Math.Pow(this.GetAlpha(i,m) * (i - 1) * h, this.GetBeta(m))); // h_{i-1/2}
                // var hph = S0 + Math.Pow(this.GetAlpha(i,m) * (i + 1) * h, this.GetBeta(m)) - xi; // h_{i+1/2}
                var hmh = h;
                var hph = h;
                var simh = xi - 0.5d * hmh; // s_{i-1/2}
                var siph = xi + 0.5d * hph; // s_{i+1/2}
                CheckHCorrectness(hph, tau, siph, r);
                CheckHCorrectness(hmh, tau, siph, r);
                var betaIm1 = 1d / (8d * tau) * (1d + 2d * tau * r * simh / hmh) * (1d + 2d * tau * r * simh / hmh);

                var betaI = 1d / (8d * tau) * (3d - 2d * tau * r * simh / hmh) * (1d + 2d * tau * r * simh / hmh)
                            + 1d / (8d * tau) * (3d + 2d * tau * r * siph / hph) * (1d - 2d * tau * r * siph / hph);

                var betaIp1 = 1d / (8d * tau) * (1d - 2d * tau * r * siph / hph) * (1d - 2d * tau * r * siph / hph);

                var f = GetF(sigma, r, GetK(), h, xi, tau, m, GetT(), S0, i);
                rp[i] = (hmh + hph) / 2d * f + (hmh + hph) / 2d * (betaIm1 * Vk1[i - 1] + betaI * Vk1[i] + betaIp1 * Vk1[i + 1]);
            }

            rp[rp.Length - 1] = 0d;

            // Console.WriteLine();
            // Console.WriteLine("Print rp on " + k);
            // var sb = new StringBuilder();
            // foreach (var d in rp.Take(5))
            // {
            //     sb.Append(d.ToString("e8", CultureInfo.InvariantCulture) + " ");
            // }
            // Console.WriteLine(k + " " + S0.ToString("e8", CultureInfo.InvariantCulture) + " " + sb);
            
            return rp;
        }

        private double[] CalculateV(double s0Old, IReadOnlyList<double> Vk1, double tau, int m, double h, out double[] rp)
        {
            rp = CalculateRightPart(s0Old, Vk1, h, tau, m);
            double[] b_t = GetB(GetN1(), s0Old, h, GetSquaredSigma(), tau);
            double[] c_t = GetC(GetN1(), s0Old, h, GetSquaredSigma(), tau, GetR());
            double[] d_t = GetD(GetN1(), s0Old, h, GetSquaredSigma(), tau);
            double[] Vk = ThomasAlgorithmCalculator.Calculate(b_t, c_t, d_t, rp);

            // PrintThomasArraysToConsole(b_t, c_t, d_t);
            return Vk;
        }

        private void CheckS0(double S0Prev, double S0Current)
        {
            if (S0Prev <= 0d)
                throw new Exception("S0New <= 0d");

            if (S0Prev >= GetK())
                throw new Exception($"S0New >= K: S0New = {S0Prev} K = {GetK()}");

            if (S0Prev >= S0Current)
                throw new Exception($"S0CurrentPrev: {S0Prev.ToString("e8", CultureInfo.InvariantCulture)} >= S0Current: {S0Current.ToString("e8", CultureInfo.InvariantCulture)}");
        }

        private void CreateFileWithStatistics()
        {
            if (!allowOutputFile)
                return;

            if (File.Exists(Path.Combine(outputPathStat, "stat.txt")))
                File.Delete(Path.Combine(outputPathStat, "stat.txt"));
        }

        private double Get_a0(double T,  double K, double smooth)
        {
            var calculateA0 = K / (smooth * T);
            return calculateA0;
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
        private double[] GetB(int n, double S0, double h, double sigmaSq, double tau)
        {
            var b = new double[n];
            for (var i = 1; i <= n - 1; ++i)
            {
                var si = S0 + i * h;
                var hmh = S0 + i * h - (S0 + (i - 1) * h); // h_{i - 1/2}
                if (!((hmh * hmh) < (4d * tau * sigmaSq * si * si)))
                    throw new ArgumentException("hmh is invalid: tau: " + tau + " sigmaSq: " + sigmaSq + " si: " + si + " si^2: " + si*si);

                b[i] = hmh / (4d * tau) - (sigmaSq * si * si) / (2d * hmh);
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
                var hph = S0 + (i + 1) * h - si; // h_{i+1/2}
                if (hmh * hmh > 4d * tau * sigma_sq * si * si)
                    throw new ArgumentException("hmh is invalid");

                c[i] = (sigma_sq * si * si) / (2d * hmh) + (sigma_sq * si * si) / (2d * hph) + (hmh + hph) * (1d / (4d * tau) + r / 2d);
            }

            // left boundary condition
            var si0 = S0 + 0 * h;
            var hph0 = S0 + (0 + 1) * h - si0; // h_{i+1/2}
            var hmh0 = si0 - (S0 + (0d - 1d) * h); // h_{i-1/2}
            var sph0 = si0 + 0.5d * hph0; // s_{i+1/2}
            var smh0 = si0 - 0.5d * hmh0; // s_{i-1/2}
            CheckHCorrectness(hph0, tau, sph0, r);
            CheckHCorrectness(hmh0, tau, smh0, r);

            c[0] = sigma_sq * si0 * si0 / (2d * hph0) + hph0 * r / 2d + hph0 / (4d * tau);

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

                d[i] = hph / (4d * tau) - sigma_sq * si * si / (2d * hph);
            }

            // left boundary condition
            var si0 = S0 + 0 * h;
            var hph0 = S0 + (0 + 1) * h - si0; // h_{i+1/2}

            d[0] = hph0 / (4d * tau) - sigma_sq * si0 * si0 / (2d * hph0);

            // right boundary condition
            d[n - 2] = 0d;
            return d;
        }

        private Point[] GetExactSolution(double tau, int m, double S0)
        {
            double t = tau * m;
            // this.UpdateH(S0);
            
            var V = new Point[GetN1()];
            var r = GetR();
            var K = GetK();
            var T = GetT();
            var sigma = GetSquaredSigma();

            var list = new List<Point>();
            var superH = GetSuperH();
            int j = 0;
            while (0d + j * superH < S0)
            {
                var S = 0d + j * superH;
                if (S > S0) break;
                list.Add(new Point(S, GetK()-S));
                j++;
            }
            
            for (int i = 0; i < GetN1() ; i++)
            {
                var S = S0 + i * GetH();
                // Vi
                var VS = (T - t)
                         * (sigma / (2d * r))
                         * Math.Pow(K / (1d + sigma / (2d * r)), (2d * r + sigma) / sigma)
                         * Math.Pow(S - Get_a0(T, K, GetSmoothness()) * t, -2d * r / sigma);

                V[i] = new Point(S, VS);
            }

            foreach (var point in V)
                list.Add(point);

            // for (int i = 0; i < this.GetN1() ; i++)
            // {
            //     // the first version (based on t*(exact_soltion_perp_am_option)) that is not working
            //     // because it is not time-dependent
            //     // var si = S0 + i * this.GetH();
            //     // var p1 = this.GetSquaredSigma() / (2d * this.GetR());
            //     //
            //     // var arg = this.GetK() / (1 + this.GetSquaredSigma() / (2d * this.GetR()));
            //     // var pow = (2d * this.GetR() + this.GetSquaredSigma()) / this.GetSquaredSigma();
            //     // var p2 = Math.Pow(arg, pow);
            //     // var arg2 = si;
            //     // var pow2 = -2d * this.GetR() / this.GetSquaredSigma();
            //     // var p3 = Math.Pow(arg2, pow2);
            //     //
            //     // var v = p1 * p2 * p3;
            //     // V[i] = tk * v;
            //     
            //     // the second version
            //     // the second version (check validation-equation-am_option.nb and validation-equation-am_option.docx)
            //     // missed
            //     var S = S0 + i * this.GetH();
            //
            //     // Vi
            //     var VS = (T - t)
            //              * (sigma / (2d * r))
            //              * Math.Pow(K / (1d + sigma / (2d * r)), (2d * r + sigma) / sigma)
            //              * Math.Pow(S - this.Get_a0(T, K, this.GetSmoothness()) * t, -2d * r / sigma);
            //              // * Math.Pow(S - S0 + (K - S0) * ((2d * r) / (sigma_sq)), (-2d * r) / sigma_sq);
            //     
            //     // the third version 
            //     // var S = S0 + i * this.GetH();
            //     //
            //     // // Vi
            //     // var Vi = (T - t)
            //     //          * (sigma_sq / (2d * r))
            //     //          * Math.Pow(K / (1d + (sigma_sq / (2d * r))), (2d * r + sigma_sq) / sigma_sq)
            //     //          * Math.Pow(S - S0 + (K - S0) * ((2d * r) / sigma_sq), (-2d * r) / sigma_sq);
            //
            //     V[j + i] = new Point(S, VS);
            // }

            return list.ToArray();
        }

        private Point[] GetExactSolution2(double tau, int m, double s0, out int startPosOfS0)
        {
            double t = tau * m;
            // this.UpdateH(s0);

            var h = GetH();
            var V = new List<Point>();
            var r = GetR();
            var K = GetK();
            var T = GetT();
            var sigma = GetSquaredSigma();
            var i = 0;
            while (i * GetH() < s0)
            {
                V.Add(new Point(
                    Geta() + i * GetH(), 
                    GetK() - i * GetH()));
                i++;
            }

            startPosOfS0 = i;
            i = 0;
            while (s0 + i * h < GetK())
            {
                var S = s0 + i * h;
                var VS = (T - t)
                        * (sigma / (2d * r))
                        * Math.Pow(K / (1d + sigma / (2d * r)), (2d * r + sigma) / sigma)
                        * Math.Pow(S - Get_a0(T, K, smoothness) * t, -2d * r / sigma);
                V.Add(new Point(S, VS));
                
                i++;
            }
            
            while (s0 + i * h < Getb())
            {
                V.Add(new Point());
                i++;
            }

            return V.ToArray();
        }

        [SuppressMessage("ReSharper", "UnusedParameter.Local")]
        private double GetF(double sigma, double r, double K, double h, double S,
             double tau, int m, double T, double S0, int i)
        {
            return 0d;
            
            //// var t = tau * k; // => k = (M-m) -> t=(M-m)*tau
            
            //// dynamic t mesh
            //var t = Math.Pow(this.GetAlpha(i,m) * m * tau, this.GetBeta(m)); // => k = (M-m) -> t=(alpha*(M-m)*tau)^beta

            //#region MyRegion

            //// return 0d;
            //// the first version (based on t*(exact_soltion_perp_am_option)) that is not working
            //// because it is not time-dependent
            //// var si = S0 + i * h;
            //// var p1 = sigma_sq / (2d * r);
            ////
            //// var arg = K / (1 + sigma_sq / (2d * r));
            //// var pow = (2d * r + sigma_sq) / sigma_sq;
            //// var p2 = Math.Pow(arg, pow);
            ////
            //// var arg2 = si;
            //// var pow2 = -2d * r / sigma_sq;
            //// var p3 = Math.Pow(arg2, pow2);
            ////
            //// var v = p1 * p2 * p3;
            //// return v;
            
            //// the second version (check validation-equation-am_option.nb and validation-equation-am_option.docx)
            
            //var a0 = this.Get_a0(T, K, this.GetSmoothness());
            
            //// dV/dt
            //var dvdtNumerator = K * Math.Pow(K / (sigma / (2d * r) + 1d), 2d * r / sigma)
            //                      * Math.Pow(S - a0 * t, -2d * r / sigma - 1d)
            //                      * (2d * a0 * r * t - 2d * a0 * r * T - a0 * sigma * t + sigma * S);
            //var dvdtDenominator = 2d * r + sigma;
            //var dvdt = -1d * (dvdtNumerator / dvdtDenominator);
            //var p1 = dvdt;
            
            //// d^2V/dS^2
            //var d2VdS2Numerator = 2d * K * r * (t - T)
            //                      * Math.Pow(K / (sigma / (2d * r) + 1), 2d * r / sigma)
            //                      * Math.Pow(S - a0 * t, -2d * r / sigma - 2d);
            //var d2VdS2Denominator = sigma;
            //var d2VdS2 = -1d * (d2VdS2Numerator / d2VdS2Denominator);
            //// (sigma^2/2)*S^2*d^2V/dS^2
            //var p2 = sigma / 2d * S * S * d2VdS2;
            
            //// dV/dS
            //var dVdS = (t - T)
            //           * Math.Pow(K / (sigma / (2d * r) + 1d), (2d * r + sigma) / sigma)
            //           * Math.Pow(S - a0 * t, -2d * r / sigma - 1d);
            //// r*S*dV/dS
            //var p3 = r * S * dVdS;
            
            //// V
            //var V = (T - t)
            //        * (sigma / (2d * r))
            //        * Math.Pow(K / (1d + sigma / (2d * r)), (2d * r + sigma) / sigma)
            //        * Math.Pow(S - a0 * t, -2d * r / sigma);
            //// r*V
            //var p4 = r * V;
            
            //return p1 + p2 + p3 - p4;
            
            //// the third version (check third-try-validation-equation-am_option.nb)
            ////
            //// // dV/dt
            //// var dVdtNumerator = -sigma
            ////                     * Math.Pow((K * r) / (2d * r + sigma), ((2d * r) / sigma) + 1d)
            ////                     * Math.Pow((r * (K - S)) / sigma, (-2d * r) / sigma);
            //// if (double.IsInfinity(Math.Pow((r * (K - S)) / sigma, (-2d * r) / sigma)))
            //// {
            ////     if (Math.Abs(K - S) > double.Epsilon)
            ////     {
            ////         throw new Exception();
            ////     }
            ////
            ////     dVdtNumerator = 0d;
            //// }
            //// var dVdtDenominator = r;
            //// var dVdt = dVdtNumerator / dVdtDenominator;
            //// var p1 = dVdt;
            ////
            //// // d^2V/dS^2
            //// var d2VdS2Numerator = -2d * K * (t - T)
            ////                       * Math.Pow((K * r) / (2d * r + sigma), (2d * r) / sigma)
            ////                       * Math.Pow((r * (K - S)) / sigma, 1d - ((2d * r) / sigma));
            //// if (double.IsInfinity(Math.Pow((r * (K - S)) / sigma, 1d - ((2d * r) / sigma))))
            //// {
            ////     if (Math.Abs(K - S) > double.Epsilon)
            ////     {
            ////         throw new Exception();
            ////     }
            ////
            ////     d2VdS2Numerator = 0d;
            //// }
            //// var d2VdS2Denominator = Math.Pow(K - S, 3d);
            //// if (Math.Abs(K - S) < double.Epsilon)
            //// {
            ////     d2VdS2Denominator = 0d;
            //// }
            ////
            //// double d2VdS2;
            //// if (Math.Abs(d2VdS2Denominator) < double.Epsilon)
            //// {
            ////     d2VdS2 = 0d;
            //// }
            //// else
            //// {
            ////     d2VdS2 = d2VdS2Numerator / d2VdS2Denominator;
            //// }
            //// // (sigma^2/2)*S^2*d^2V/dS^2
            //// var p2 = (sigma / 2d) * S * S * d2VdS2;
            ////
            //// // dV/dS
            //// var dVdSNumerator = -2d * r
            ////                         * (t - T)
            ////                         * Math.Pow((K * r) / (2d * r + sigma), ((2d * r) / sigma) + 1d)
            ////                         * Math.Pow((r * (K - S)) / sigma, ((-2d * r) / sigma) - 1d);
            //// if (double.IsInfinity(Math.Pow((r * (K - S)) / (sigma), ((-2d * r) / sigma) - 1d)))
            //// {
            ////     if (Math.Abs(K - S) > double.Epsilon)
            ////     {
            ////         throw new Exception();
            ////     }
            ////
            ////     dVdSNumerator = 0d;
            //// }
            ////
            //// var dVdSDenominator = sigma;
            //// // r*S*dV/dS
            //// var dVdS = dVdSNumerator / dVdSDenominator;
            //// var p3 = r * S * dVdS;
            ////
            //// // V
            //// var V = (T - t)
            ////         * (sigma / (2d * r))
            ////         * Math.Pow(K / (1d + (sigma / (2d * r))), (2d * r + sigma) / sigma)
            ////         * Math.Pow(((2d*r)/(sigma))*(K-S), (-2d * r) / sigma);
            //// if (Math.Abs((K - S)) < double.Epsilon)
            //// {
            ////     if (Math.Abs(K - S) > double.Epsilon)
            ////     {
            ////         throw new Exception();
            ////     }
            ////
            ////     V = 0d;
            //// }
            //// // r*V
            //// var p4 = r * V;
            ////
            //// return p1 + p2 + p3 - p4;

            //#endregion
            
            //// // the third updated version (check third-try-validation-equation-am_option.nb)
            ////
            //// // dV/dt
            //// var dVdtNumerator = -K * sigma
            ////                        * Math.Pow(K / ((sigma / (2d * r)) + 1d), (2d * r) / sigma)
            ////                        * Math.Pow(( (2d * r * (K - S0)) / sigma) - S0 + S, (-2d * r) / sigma);
            //// var dVdtDenominator = 2d * r + sigma;
            //// var dVdt = dVdtNumerator / dVdtDenominator;
            //// var p1 = dVdt;
            ////
            //// // d^2V/dS^2
            //// var d2VdS2 = ( ((-2d * r) / sigma) - 1d)
            ////              * (t-T) // = -(T-t)
            ////              * Math.Pow(K / ((sigma / (2d * r)) + 1d), (2d * r + sigma) / sigma)
            ////              * Math.Pow( ((2d * r * (K - S0)) / sigma) - S0 + S, ((-2d * r) / sigma) - 2d);
            //// // (sigma^2/2)*S^2*d^2V/dS^2
            //// var p2 = (sigma / 2d) * S * S * d2VdS2;
            ////
            //// // dV/dS
            //// var dVdS = (t - T)
            ////            * Math.Pow(K / ((sigma / (2d * r)) + 1d), (2d * r + sigma) / sigma)
            ////            * Math.Pow(((2d * r * (K - S0)) / sigma) - S0 + S, ((-2d * r) / sigma) - 1d);
            //// // r*S*dV/dS
            //// var p3 = r * S * dVdS;
            ////
            //// // V
            //// var V = (T - t)
            ////         * (sigma / (2d * r))
            ////         * Math.Pow(K / (1d + (sigma / (2d * r))), (2d * r + sigma) / sigma)
            ////         * Math.Pow((S - S0 + (K - S0) * ((2d * r) / sigma)), (-2d * r) / sigma);
            //// // r*V
            //// var p4 = r * V;
            
            //// return p1 + p2 + p3 - p4;
            //// return 0;
        }

        private double GetSmoothness()
        {
            return smoothness;
        }

        private double[] GetVST()
        {
            // this.UpdateH(this.GetK());
            var arr = new double[GetN1()];
            for (var i = 0; i < arr.Length; ++i)
            {
                var S = GetK() + i * GetH();
                if (S <= GetK())
                    arr[i] = GetK() - S;
                else
                    arr[i] = 0d;
            }

            return arr;
        }

        private Point[] GetVonT()
        {
            var superH = GetSuperH();
            
            var list = new List<Point>();

            var s = 0d;
            int i = 0;
            while (s < GetK())
            {
                s = 0d + i * superH;
                if (s > GetK())
                    break;

                list.Add(new Point(s, GetK()-s));
                i++;
            }

            for (int j = 0; j < GetN1(); j++)
            {
                s = GetK() + j * superH;
                if (s > Getb())
                    break;

                list.Add(new Point(s, 0));
            }

            return list.ToArray();
        }

        private double GetSuperH()
        {
            return GetH() / 1d;
        }

        private void PrintData(double s0Current, TecplotPrinter tecplotPrinter, double[] rp, int iter, int m, double S0Next, IReadOnlyList<double> Vk)
        {
            PrintRpToTecplot(tecplotPrinter, rp, s0Current);
            PrintStatistics(iter, m, GetH(), S0Next, s0Current);
            PrintValuesToConsole(iter, GetH(), S0Next, Vk, s0Current);
        }

        private void PrintHeader(IReadOnlyList<double> St)
        {
            if (!allowOutputConsole)
                return;

            Console.WriteLine($"Number of time steps = {GetM()} h = {GetH()} S(T) = {St[St.Count - 1]}");
            Console.WriteLine("--------------------------------------------------");
        }

        private void PrintRpToTecplot(TecplotPrinter tecplotPrinter, double[] rp, double S0New)
        {
            if (allowOutputFile)
                tecplotPrinter.PrintXY(Path.Combine(outputPathRp, "temporal-rp"), 0d, 0, GetH(), rp, S0New);
        }

        private void PrintStatistics(int iter, int m, double h_old, double S0Old, double S0New)
        {
            if (!allowOutputFile)
                return;

            using (var streamWriter = File.AppendText(Path.Combine(outputPathStat, "stat.txt")))
            {
                streamWriter.WriteLine(
                    new string(' ', 2) + " Time step = {0} Iteration = " + iter + " h = {1} S0 = {2} Abs(S0New-S0Old)={3} S0Eps={4} Cnd={5}",
                    m,
                    h_old,
                    S0Old,
                    Math.Abs(S0New - S0Old),
                    GetS0Eps(),
                    Math.Abs(S0New - S0Old) > GetS0Eps());
            }
        }

        private void PrintValuesToConsole(int iter, double h_old, double S0Old, IReadOnlyList<double> Vk, double S0New)
        {
            if (allowOutputConsole)
                Console.WriteLine(new string(' ', 2) + "Iteration = " + iter + " h = {0} S0_old = {1} Vk[0] = {2} S0_new = {3}", h_old, S0Old, Vk[0], S0New);
        }

        private void PrintVST(TecplotPrinter tecplotPrinter, double[] Vk1, IReadOnlyList<double> St)
        {
            if (allowOutputFile)
                tecplotPrinter.PrintXY(outputPath + "V" + GetM(), GetTau(), GetM(), GetH(), Vk1, St[St.Count - 1]);
        }

        private void SaveNumericSolution(IReadOnlyList<double> V, int m, double S0)
        {
            var temp = new Point[V.Count];
            for (var i = 0; i < V.Count; i++)
                temp[i] = new Point(S0+i*GetH(),V[i]);

            solutions.Add(new SolutionData(S0, m)
            {
                Solution = temp
            });
        }
    }
}