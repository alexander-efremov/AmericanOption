namespace AmericanOptions
{
    using System;
    using System.Collections.Generic;
    using System.Diagnostics.CodeAnalysis;
    using System.Globalization;
    using System.IO;
    using System.Linq;
    using System.Text;
    using CoreLib;

    using AmericanOption;

    // from new presentation with FEM
    [SuppressMessage("ReSharper", "CommentTypo", Justification = "OK.")]
    public class AmericanOptionCalculator : AmericanOptionCalculatorBase
    {
        private readonly bool allowOutputFile;

        private readonly bool allowOutputConsole;

        private readonly string outputPath;

        private readonly string outputPathStat;

        private readonly string outputPathRp;

        private readonly List<SolutionData> solutions = new List<SolutionData>();

        private readonly bool saveSolutions;

        private readonly int MValue;

        private readonly double TValue;
        
        public AmericanOptionCalculator(AmericanOptionParameters parameters, bool allowOutputFile, bool allowOutputConsole)
            : base(parameters)
        {
            this.allowOutputFile = allowOutputFile;
            this.allowOutputConsole = allowOutputConsole;
            this.outputPath = parameters.WorkDir;
            if (string.IsNullOrEmpty(this.outputPath))
            {
                this.outputPath = this.GetWorkDir() + "AO/";
            }

            this.MValue = parameters.M;
            this.TValue = parameters.T;
            this.saveSolutions = parameters.SaveVSolutions;

            this.outputPathStat = Path.Combine(this.outputPath, "stat");
            if (!Directory.Exists(this.outputPathStat))
            {
                Directory.CreateDirectory(this.outputPathStat);
            }

            this.outputPathRp = Path.Combine(this.outputPath, "rp");
            if (!Directory.Exists(this.outputPathRp))
            {
                Directory.CreateDirectory(this.outputPathRp);
            }
        }

        public int GetM()
        {
            return this.MValue;
        }

        public double GetT()
        {
            return this.TValue;
        }

        /// <summary>
        /// Returns numeric solutions V(S(t)) in reversed order (from the M to the 0 time step)
        /// </summary>
        /// <returns></returns>
        public List<SolutionData> GetNumericSolutions()
        {
            // solutions were filled in reversed order from the M to 0
            // then here we reverse it
            List<SolutionData> list = this.solutions.ToList();
            return list;
        }

        /// <summary>
        /// Solves the american option problem
        /// </summary>
        /// <returns>Array of S0 values in REVERSED order (S0t[0] refers to the M time step, S0t[0] refers to the 0 time step)</returns>
        public double[] Solve()
        {
            // var printer = new ThomasArrayPrinter();
            double[] VNext = this.GetVST(); // V(S, T) = (K - S)+

            var printer = new TecplotPrinterSpecial(0d, this.GetRightBoundary(), this.GetTau());
            // tecplotPrinter.PrintXY(Path.Combine(_outputPath, "VST"), 0d, GetH(), Vk1);
            
            var S0t = new double[this.GetM() + 1];
            S0t[this.GetM()] = this.GetK();
            
            this.SaveNumericSolutions(VNext, this.GetM(), S0t[this.GetM()]);

            this.PrintHeader(S0t);
            this.CreateFileWithStatistics();
            this.PrintVST(printer, VNext, S0t);

            for (var k = this.GetM() - 1; k >= 0; --k)
            {
                if (this.allowOutputConsole)
                {
                    Console.WriteLine("Time step = " + k);
                }

                double[] VCurrent;
                var iter = 0;
                var S0Next = S0t[k + 1];
                double S0Current = S0Next;
                double S0Diff;
                do
                {
                    this.UpdateH(S0Current);
                    
                    // calculate new V(S, t_k)
                    VCurrent = this.CalculateV(S0Current, VNext, this.GetTau(), k, out double[] rp);

                    // calculate new S0(t)_i
                    var S0CurrentPrev = S0Current;
                    S0Diff = Math.Abs((this.GetK() - VCurrent[0]) - S0Current);
                    S0Current = this.GetK() - VCurrent[0];
                    this.CheckS0(S0Current, printer, rp, ++iter, k,  S0CurrentPrev, VCurrent);
                } while (S0Diff > this.GetS0Eps());

                // TODO: why we get some negative values instead of zero?
                RemoveNegativeValues(VCurrent);
                
                this.SaveNumericSolutions(VCurrent, k, S0Current);
                S0t[k] = S0Current;
                for (var i = 0; i < VCurrent.Length; i++)
                {
                    VNext[i] = VCurrent[i];
                }

                if (this.allowOutputFile)
                {
                    printer.PrintXY(Path.Combine(this.outputPath, "V" + k), this.GetTau() * k, this.GetH(), VCurrent, S0Current);
                }

                if (this.allowOutputConsole)
                {
                    Console.WriteLine("--------------------------------------------------");
                }
            }

            return S0t;
        }

        private static void RemoveNegativeValues(double[] VCurrent)
        {
            for (var i = 0; i < VCurrent.Length; i++)
            {
                if (VCurrent[i] < 0d)
                {
                    VCurrent[i] = 0d;
                }
            }
        }

        public TecplotPrinterSpecial GetTecplotPrinter()
        {
            return new TecplotPrinterSpecial(0d, this.GetRightBoundary(), this.GetTau());
        }
        
        /// <summary>
        /// Returns exact solutions V(S(t)) in reversed order (from the M to the 0 time step)
        /// </summary>
        /// <param name="S0t">S0 in direct order from 0 to M time layer</param>
        /// <returns></returns>
        public List<SolutionData> GetExactSolutions(double[] S0t)
        {
            if (Math.Abs(this.GetM() * this.GetTau() - this.GetT()) > double.Epsilon)
            {
                throw new Exception("GetExactSolutions");
            }

            var list = new List<SolutionData>
            {
                new SolutionData(S0t, this.GetM())
                {
                    Solution = this.GetVST() // V(S, T) = (K - S)+
                }
            };
            
            for (var k = this.GetM() - 1; k >= 0; k--)
            {
                double[] solution = this.GetExactSolution(k * this.GetTau(), S0t[k]);
                var solutionData = new SolutionData(S0t, k)
                {
                    Solution = solution
                };
                list.Add(solutionData);
            }

            return list;
        }

        private double[] GetExactSolution(double t, double S0)
        {
            var V = new double[this.GetN1()];
            var r = this.GetR();
            var K = this.GetK();
            var T = this.GetT();
            var sigma_sq = this.GetSquaredSigma();
            UpdateH(S0);
            for (var i = 0; i < V.Length; i++)
            {
                // the first version (based on t*(exact_soltion_perp_am_option)) that is not working
                // because it is not time-dependent
                // var si = S0 + i * this.GetH();
                // var p1 = this.GetSquaredSigma() / (2d * this.GetR());
                //
                // var arg = this.GetK() / (1 + this.GetSquaredSigma() / (2d * this.GetR()));
                // var pow = (2d * this.GetR() + this.GetSquaredSigma()) / this.GetSquaredSigma();
                // var p2 = Math.Pow(arg, pow);
                // var arg2 = si;
                // var pow2 = -2d * this.GetR() / this.GetSquaredSigma();
                // var p3 = Math.Pow(arg2, pow2);
                //
                // var v = p1 * p2 * p3;
                // V[i] = tk * v;
                
                // the second version
                // the second version (check validation-equation-am_option.nb and validation-equation-am_option.docx)
                // missed
                var S = S0 + i * this.GetH();

                // Vi
                var Vi = (T - t)
                         * (sigma_sq / (2d * r))
                         * Math.Pow(K / (1d + (sigma_sq / (2d * r))), (2d * r + sigma_sq) / sigma_sq)
                         * Math.Pow(S - (this.Get_a0(T, K) * t), (-2d * r) / sigma_sq);
                         // * Math.Pow(S - S0 + (K - S0) * ((2d * r) / (sigma_sq)), (-2d * r) / sigma_sq);
                
                // the third version 
                // var S = S0 + i * this.GetH();
                //
                // // Vi
                // var Vi = (T - t)
                //          * (sigma_sq / (2d * r))
                //          * Math.Pow(K / (1d + (sigma_sq / (2d * r))), (2d * r + sigma_sq) / sigma_sq)
                //          * Math.Pow(S - S0 + (K - S0) * ((2d * r) / sigma_sq), (-2d * r) / sigma_sq);

                V[i] = Vi;
            }
            
            // TODO: why we get some negative values instead of zero?
            RemoveNegativeValues(V);

            return V;
        }

        private void PrintVST(TecplotPrinter tecplotPrinter, double[] Vk1, IReadOnlyList<double> St)
        {
            if (this.allowOutputFile)
            {
                tecplotPrinter.PrintXY(this.outputPath + "V" + this.GetM(), this.GetTau() * this.GetM(), this.GetH(), Vk1, St[St.Count - 1]);
            }
        }

        private double[] GetVST()
        {
            this.UpdateH(this.GetK());
            var arr = new double[this.GetN1()];
            for (var i = 0; i < arr.Length; ++i)
            {
                var S = this.GetK() + i * this.GetH();
                if (S <= this.GetK())
                {
                    arr[i] = this.GetK() - S;
                }
                else
                {
                    arr[i] = 0d;
                }
            }

            return arr;
        }

        private void PrintHeader(IReadOnlyList<double> St)
        {
            if (!this.allowOutputConsole)
            {
                return;
            }

            Console.WriteLine($"Number of time steps = {this.GetM()} h = {this.GetH()} S(T) = {St[St.Count - 1]}");
            Console.WriteLine("--------------------------------------------------");
        }

        private void CreateFileWithStatistics()
        {
            if (!this.allowOutputFile)
            {
                return;
            }

            if (File.Exists(Path.Combine(this.outputPathStat, "stat.txt")))
            {
                File.Delete(Path.Combine(this.outputPathStat, "stat.txt"));
            }
        }

        private void SaveNumericSolutions(IReadOnlyList<double> Vk1, int k, double s0)
        {
            if (!this.saveSolutions)
            {
                return;
            }

            var r = new double[Vk1.Count];
            for (var i = 0; i < Vk1.Count; i++)
            {
                r[i] = Vk1[i];
            }

            this.solutions.Add(new SolutionData(s0, k)
            {
                Solution = r
            });
        }

        private double[] CalculateV(double s0Old, IReadOnlyList<double> Vk1, double tau, int k, out double[] rp)
        {
            rp = this.CalculateRightPart(s0Old, Vk1, this.GetH(), tau, k);

            double[] b_t = this.GetB(this.GetN1(), s0Old, this.GetH(), this.GetSquaredSigma(), tau);
            double[] c_t = this.GetC(this.GetN1(), s0Old, this.GetH(), this.GetSquaredSigma(), tau, this.GetR());
            double[] d_t = this.GetD(this.GetN1(), s0Old, this.GetH(), this.GetSquaredSigma(), tau);
            double[] Vk = this.ThomasAlgorithmCalculator.Calculate(b_t, c_t, d_t, rp);

            // PrintThomasArraysToConsole(b_t, c_t, d_t);
            return Vk;
        }

        private void PrintStatistics(int iter, int k, double h_old, double S0Old, double S0New)
        {
            if (!this.allowOutputFile)
            {
                return;
            }

            using (var streamWriter = File.AppendText(Path.Combine(this.outputPathStat, "stat.txt")))
            {
                streamWriter.WriteLine(
                    new string(' ', 2) + " Time step = {0} Iteration = " + iter + " h = {1} S0 = {2} Abs(S0New-S0Old)={3} S0Eps={4} Cnd={5}",
                    k,
                    h_old,
                    S0Old,
                    Math.Abs(S0New - S0Old),
                    this.GetS0Eps(),
                    Math.Abs(S0New - S0Old) > this.GetS0Eps());
            }
        }

        private void PrintValuesToConsole(int iter, double h_old, double S0Old, IReadOnlyList<double> Vk, double S0New)
        {
            if (this.allowOutputConsole)
            {
                Console.WriteLine(new string(' ', 2) + "Iteration = " + iter + " h = {0} S0_old = {1} Vk[0] = {2} S0_new = {3}", h_old, S0Old, Vk[0], S0New);
            }
        }

        private void PrintRpToTecplot(TecplotPrinter tecplotPrinter, double[] rp, double S0New)
        {
            if (this.allowOutputFile)
            {
                tecplotPrinter.PrintXY(Path.Combine(this.outputPathRp, "temporal-rp"), 0d, this.GetH(), rp, S0New);
            }
        }

        private double[] CalculateRightPart(double S0, IReadOnlyList<double> Vk1, double h, double tau, int k)
        {
            var rp = new double[this.GetN1()];
            var r = this.GetR();
            var sigma = this.GetSquaredSigma();

            var hmh0 = S0 - (S0 + h); // h_{i-1/2}
            var hph0 = (S0 + h) - S0; // h_{i+1/2}
            var smh0 = S0 - 0.5d * hmh0; // s_{i-1/2}
            var sph0 = S0 + 0.5d * hph0; // s_{i+1/2}
            CheckHCorrectness(hmh0, tau, sph0, r);
            CheckHCorrectness(hph0, tau, sph0, r);
            var betam1 = 1d / (8d * tau) * (1d + (2d * tau * r * smh0) / hmh0) * (1d + (2d * tau * r * smh0) / hmh0);
            var beta0 =  1d / (8d * tau) * (3d - 2d * tau * r * smh0 / hmh0) * (1d + 2d * tau * r * smh0 / hmh0) 
                + 1d / (8d * tau) * (3d + (2d * tau * r * sph0) / hph0) * (1d - (2d * tau * r * sph0) / hph0);
            var betap1 = 1d / (8d * tau) * (1d - (2d * tau * r * sph0) / hph0) * (1d - (2d * tau * r * sph0) / hph0);
            var f0 = this.GetF(sigma, r, this.GetK(), h, S0, this.GetTau(), k, this.GetT(), S0);

            rp[0] = ((-hph0 / 2d) * f0) + ((sigma * S0 * S0) / 2d)
                                        + (hph0 / 2d) * (betam1 * 0d /*Vk1[-1]*/ + beta0 * Vk1[0] + betap1 * Vk1[1]);

            for (var i = 1; i < rp.Length - 1; ++i)
            {
                var si = S0 + i * h;
                var hmh = si - (S0 + (i - 1) * h); // h_{i-1/2}
                var hph = (S0 + (i + 1) * h) - si; // h_{i+1/2}
                var simh = si - 0.5d * hmh; // s_{i-1/2}
                var siph = si + 0.5d * hph; // s_{i+1/2}
                CheckHCorrectness(hph, tau, siph, r);
                CheckHCorrectness(hmh, tau, siph, r);
                var betaIm1 = 1d / (8d * tau) * (1d + 2d * tau * r * simh / hmh) * (1d + 2d * tau * r * simh / hmh);

                var betaI = 1d / (8d * tau) * (3d - 2d * tau * r * simh / hmh) * (1d + 2d * tau * r * simh / hmh)
                            + 1d / (8d * tau) * (3d + 2d * tau * r * siph / hph) * (1d - 2d * tau * r * siph / hph);

                var betaIp1 = 1d / (8d * tau) * (1d - 2d * tau * r * siph / hph) * (1d - 2d * tau * r * siph / hph);

                var f = this.GetF(sigma, r, this.GetK(), h, si, tau, k, this.GetT(), S0);
                rp[i] = ((hmh + hph) / 2d) * f + ((hmh + hph) / 2d) * (betaIm1 * Vk1[i - 1] + betaI * Vk1[i] + betaIp1 * Vk1[i + 1]);
            }

            rp[rp.Length - 1] = 0d;

            Console.WriteLine();
            Console.WriteLine("Print rp on " + k);
            var sb = new StringBuilder();
            foreach (var d in rp.Take(5))
            {
                sb.Append(d.ToString("e8", CultureInfo.InvariantCulture) + " ");
            }

            Console.WriteLine(k + " " + S0.ToString("e8", CultureInfo.InvariantCulture) + " " + sb);
            return rp;
        }

        [SuppressMessage("ReSharper", "UnusedParameter.Local")]
        private double GetF(double sigma, double r, double K, double h, double S,
             double tau, int k, double T, double S0)
        {
            #region MyRegion

            // return 0d;
            // the first version (based on t*(exact_soltion_perp_am_option)) that is not working
            // because it is not time-dependent
            // var si = S0 + i * h;
            // var p1 = sigma_sq / (2d * r);
            //
            // var arg = K / (1 + sigma_sq / (2d * r));
            // var pow = (2d * r + sigma_sq) / sigma_sq;
            // var p2 = Math.Pow(arg, pow);
            //
            // var arg2 = si;
            // var pow2 = -2d * r / sigma_sq;
            // var p3 = Math.Pow(arg2, pow2);
            //
            // var v = p1 * p2 * p3;
            // return v;
            
            // the second version (check validation-equation-am_option.nb and validation-equation-am_option.docx)
            var t = tau * k;
            var a0 = this.Get_a0(T, K);
            
            // dV/dt
            var dvdtNumerator = K * Math.Pow(K / ((sigma / (2d * r)) + 1d), (2d * r) / sigma)
                                  * Math.Pow(S - a0 * t, ((-2d * r) / sigma) - 1d)
                                  * (2d * a0 * r * t - 2d * a0 * r * T - a0 * sigma * t + sigma * S);
            var dvdtDenominator = 2d * r + sigma;
            var dvdt = -1d * (dvdtNumerator / dvdtDenominator);
            var p1 = dvdt;
            
            // d^2V/dS^2
            var d2VdS2Numerator = 2d * K * r * (t - T)
                                  * Math.Pow(K / ((sigma / (2d * r)) + 1), (2d * r) / sigma)
                                  * Math.Pow(S - a0 * t, ((-2d * r) / sigma) - 2d);
            var d2VdS2Denominator = sigma;
            var d2VdS2 = (-1d * (d2VdS2Numerator / d2VdS2Denominator));
            // (sigma^2/2)*S^2*d^2V/dS^2
            var p2 = (sigma / 2d) * S * S * d2VdS2;
            
            // dV/dS
            var dVdS = (t - T)
                       * Math.Pow(K / ((sigma / (2d * r)) + 1d), (2d * r + sigma) / sigma)
                       * Math.Pow(S - a0 * t, ((-2d * r) / sigma) - 1d);
            // r*S*dV/dS
            var p3 = r * S * dVdS;
            
            // V
            var V = (T - t)
                    * (sigma / (2d * r))
                    * Math.Pow(K / (1d + (sigma / (2d * r))), (2d * r + sigma) / sigma)
                    * Math.Pow(S - a0 * t, (-2d * r) / sigma);
            // r*V
            var p4 = r * V;
            
            return p1 + p2 + p3 - p4;
            
            // the third version (check third-try-validation-equation-am_option.nb)
            // var t = tau * k;
            //
            // // dV/dt
            // var dVdtNumerator = -sigma
            //                     * Math.Pow((K * r) / (2d * r + sigma), ((2d * r) / sigma) + 1d)
            //                     * Math.Pow((r * (K - S)) / sigma, (-2d * r) / sigma);
            // if (double.IsInfinity(Math.Pow((r * (K - S)) / sigma, (-2d * r) / sigma)))
            // {
            //     if (Math.Abs(K - S) > double.Epsilon)
            //     {
            //         throw new Exception();
            //     }
            //
            //     dVdtNumerator = 0d;
            // }
            // var dVdtDenominator = r;
            // var dVdt = dVdtNumerator / dVdtDenominator;
            // var p1 = dVdt;
            //
            // // d^2V/dS^2
            // var d2VdS2Numerator = -2d * K * (t - T)
            //                       * Math.Pow((K * r) / (2d * r + sigma), (2d * r) / sigma)
            //                       * Math.Pow((r * (K - S)) / sigma, 1d - ((2d * r) / sigma));
            // if (double.IsInfinity(Math.Pow((r * (K - S)) / sigma, 1d - ((2d * r) / sigma))))
            // {
            //     if (Math.Abs(K - S) > double.Epsilon)
            //     {
            //         throw new Exception();
            //     }
            //
            //     d2VdS2Numerator = 0d;
            // }
            // var d2VdS2Denominator = Math.Pow(K - S, 3d);
            // if (Math.Abs(K - S) < double.Epsilon)
            // {
            //     d2VdS2Denominator = 0d;
            // }
            //
            // double d2VdS2;
            // if (Math.Abs(d2VdS2Denominator) < double.Epsilon)
            // {
            //     d2VdS2 = 0d;
            // }
            // else
            // {
            //     d2VdS2 = d2VdS2Numerator / d2VdS2Denominator;
            // }
            // // (sigma^2/2)*S^2*d^2V/dS^2
            // var p2 = (sigma / 2d) * S * S * d2VdS2;
            //
            // // dV/dS
            // var dVdSNumerator = -2d * r
            //                         * (t - T)
            //                         * Math.Pow((K * r) / (2d * r + sigma), ((2d * r) / sigma) + 1d)
            //                         * Math.Pow((r * (K - S)) / sigma, ((-2d * r) / sigma) - 1d);
            // if (double.IsInfinity(Math.Pow((r * (K - S)) / (sigma), ((-2d * r) / sigma) - 1d)))
            // {
            //     if (Math.Abs(K - S) > double.Epsilon)
            //     {
            //         throw new Exception();
            //     }
            //
            //     dVdSNumerator = 0d;
            // }
            //
            // var dVdSDenominator = sigma;
            // // r*S*dV/dS
            // var dVdS = dVdSNumerator / dVdSDenominator;
            // var p3 = r * S * dVdS;
            //
            // // V
            // var V = (T - t)
            //         * (sigma / (2d * r))
            //         * Math.Pow(K / (1d + (sigma / (2d * r))), (2d * r + sigma) / sigma)
            //         * Math.Pow(((2d*r)/(sigma))*(K-S), (-2d * r) / sigma);
            // if (Math.Abs((K - S)) < double.Epsilon)
            // {
            //     if (Math.Abs(K - S) > double.Epsilon)
            //     {
            //         throw new Exception();
            //     }
            //
            //     V = 0d;
            // }
            // // r*V
            // var p4 = r * V;
            //
            // return p1 + p2 + p3 - p4;

            #endregion
            
            // // the third updated version (check third-try-validation-equation-am_option.nb)
            // var t = tau * k;
            //
            // // dV/dt
            // var dVdtNumerator = -K * sigma
            //                        * Math.Pow(K / ((sigma / (2d * r)) + 1d), (2d * r) / sigma)
            //                        * Math.Pow(( (2d * r * (K - S0)) / sigma) - S0 + S, (-2d * r) / sigma);
            // var dVdtDenominator = 2d * r + sigma;
            // var dVdt = dVdtNumerator / dVdtDenominator;
            // var p1 = dVdt;
            //
            // // d^2V/dS^2
            // var d2VdS2 = ( ((-2d * r) / sigma) - 1d)
            //              * (t-T) // = -(T-t)
            //              * Math.Pow(K / ((sigma / (2d * r)) + 1d), (2d * r + sigma) / sigma)
            //              * Math.Pow( ((2d * r * (K - S0)) / sigma) - S0 + S, ((-2d * r) / sigma) - 2d);
            // // (sigma^2/2)*S^2*d^2V/dS^2
            // var p2 = (sigma / 2d) * S * S * d2VdS2;
            //
            // // dV/dS
            // var dVdS = (t - T)
            //            * Math.Pow(K / ((sigma / (2d * r)) + 1d), (2d * r + sigma) / sigma)
            //            * Math.Pow(((2d * r * (K - S0)) / sigma) - S0 + S, ((-2d * r) / sigma) - 1d);
            // // r*S*dV/dS
            // var p3 = r * S * dVdS;
            //
            // // V
            // var V = (T - t)
            //         * (sigma / (2d * r))
            //         * Math.Pow(K / (1d + (sigma / (2d * r))), (2d * r + sigma) / sigma)
            //         * Math.Pow((S - S0 + (K - S0) * ((2d * r) / sigma)), (-2d * r) / sigma);
            // // r*V
            // var p4 = r * V;
            
            

            // return p1 + p2 + p3 - p4;
            // return 0;
        }

        private double Get_a0(double T,  double K)
        {
            var calculateA0 = K / (2d * T);
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
                var hmh = (S0 + i * h) - (S0 + (i - 1) * h); // h_{i - 1/2}
                if (hmh * hmh > 4d * tau * sigmaSq * si * si)
                {
                    throw new ArgumentException("hmh is invalid");
                }

                b[i] = (hmh / (4d * tau)) - (sigmaSq * si * si) / (2d * hmh);
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
                if (hmh * hmh > 4d * tau * sigma_sq * si * si)
                {
                    throw new ArgumentException("hmh is invalid");
                }

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
                // ReSharper disable once NotResolvedInText
                throw new ArgumentOutOfRangeException("tau/h");
            }
        }
        
        private void CheckS0(double s0Current, TecplotPrinterSpecial tecplotPrinter, double[] rp, int iter, int k, double S0CurrentPrev, double[] Vk)
        {
            try
            {
                if (s0Current <= 0d)
                {
                    throw new Exception("S0New <= 0d");
                }

                if (s0Current >= this.GetK())
                {
                    throw new Exception($"S0New >= K: S0New = {s0Current} K = {GetK()}");
                }
            }
            finally
            {
                this.PrintRpToTecplot(tecplotPrinter, rp, s0Current);
                this.PrintStatistics(iter, k, this.GetH(), S0CurrentPrev, s0Current);
                this.PrintValuesToConsole(iter, this.GetH(), S0CurrentPrev, Vk, s0Current);
            }
        }
    }

    public class SolutionData
    {
        public SolutionData(IReadOnlyList<double> S0arr, int k)
        {
            this.k = k;
            S0 = S0arr[k];
        }
        
        public SolutionData(double S0, int k)
        {
            this.k = k;
            this.S0 = S0;
        }

        public double[] Solution { get; set; }
        public double S0 { get; private set; }
        public double k { get; private set; }
    }
}