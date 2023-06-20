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
        private readonly bool _allowOutputConsole;

        private readonly bool _allowOutputFile;

        private readonly int _mValue;

        private readonly string _outputPath;

        private readonly string _outputPathRp;

        private readonly string _outputPathStat;

        private readonly bool _saveSolutions;

        private readonly double _smoothness;

        private readonly List<SolutionData> _solutions = new List<SolutionData>();

        private readonly double _value;

        public AmericanOptionCalculator(AmericanOptionParameters parameters, bool allowOutputFile, bool allowOutputConsole)
            : base(parameters)
        {
            _allowOutputFile = allowOutputFile;
            _allowOutputConsole = allowOutputConsole;
            _outputPath = parameters.WorkDir;
            if (string.IsNullOrEmpty(_outputPath))
                _outputPath = GetWorkDir() + "AO/";

            _mValue = parameters.M;
            _value = parameters.T;
            _smoothness = parameters.Smoothness;
            _saveSolutions = parameters.SaveVSolutions;

            _outputPathStat = Path.Combine(_outputPath, "stat");
            if (!Directory.Exists(_outputPathStat))
                Directory.CreateDirectory(_outputPathStat);

            _outputPathRp = Path.Combine(_outputPath, "rp");
            if (!Directory.Exists(_outputPathRp))
                Directory.CreateDirectory(_outputPathRp);
        }

        /// <summary>
        /// Returns exact solutions V(S(t)) in reversed order (from the M to the 0 time step)
        /// </summary>
        /// <param name="s0T">S0 in direct order from 0 to M time layer</param>
        /// <returns></returns>
        public List<SolutionData> GetExactSolutions(double[] s0T)
        {
            if (Math.Abs(GetM() * GetTau() - GetT()) > double.Epsilon)
                throw new Exception("GetExactSolutions");

            var list = new List<SolutionData>
            {
                new SolutionData(s0T, GetM())
                {
                    Solution = GetVonT() // V(S, T) = (K - S)+
                }
            };

            for (var k = GetM() - 1; k >= 0; k--)
            {
                Point[] solution = GetExactSolution(GetTau(), k, s0T[k]);
                var solutionData = new SolutionData(s0T, k)
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
        public List<SolutionData> GetExactSolutions2(double[] s0T)
        {
            if (Math.Abs(GetM() * GetTau() - GetT()) > double.Epsilon)
                throw new Exception("GetExactSolutions");

            var solution2 = GetExactSolution2(GetM(), GetM(), 0d, out var startPosOfS01);
            var item = new SolutionData(0d, GetM()) {Solution = solution2, StartPosOfS0 = startPosOfS01};
            var list = new List<SolutionData> {item};
            for (var k = GetM() - 1; k >= 0; k--)
            {
                var solution = GetExactSolution2(GetTau(), k, s0T[k], out var startPosOfS02);
                var solutionData = new SolutionData(s0T, k)
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
            return _mValue;
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

            foreach (var data in _solutions)
            {
                Point[] array = data.Solution.ToArray();
                var solutionData = new SolutionData(data.S0, data.K)
                {
                    Solution = array
                };
                list.Add(solutionData);
            }
            
            return list;
        }

        public double GetT()
        {
            return _value;
        }

        /// <summary>
        /// Solves the american option problem
        /// </summary>
        /// <returns>Array of S0 values in REVERSED order (S0t[0] refers to the M time step, S0t[0] refers to the 0 time step)</returns>
        public double[] Solve()
        {
            // var printer = new ThomasArrayPrinter();
            double[] vNext = GetVst(); // V(S, T) = (K - S)+

            var printer = GetTecplotPrinter();
            // tecplotPrinter.PrintXY(Path.Combine(_outputPath, "VST"), 0d, GetH(), V);
            
            var s0T = new double[GetM() + 1];
            s0T[GetM()] = GetK();
            
            if (_saveSolutions)
                SaveNumericSolution(vNext, GetM(), s0T[GetM()]);

            PrintHeader(s0T);
            CreateFileWithStatistics();
            PrintVst(printer, vNext, s0T);

            for (var m = GetM() - 1; m >= 0; --m)
            {
                var tau = GetTau();
                if (_allowOutputConsole)
                    Console.WriteLine("Time step = " + m);

                double[] vCurrent;
                var iter = 0;
                var s0Next = s0T[m + 1];
                var s0Current = s0Next;
                double s0Diff;
                if (_allowOutputConsole)
                {
                    Console.WriteLine("m = " + m);
                    Console.WriteLine("S0Current = " + s0Current);
                }
                
                //var h = GetH();
                var h = (Getb() - GetK()) / GetN();
                
                do
                {
                    // calculate new V(S, t_k)
                    vCurrent = CalculateV(s0Current, vNext, tau, m, h, out double[] rp);

                    if (vCurrent[0] <= 0d)
                        throw new Exception("VCurrent[0] <= 0d");

                    // calculate new S0(t)_i
                    var s0CurrentPrev = s0Current;
                    s0Current = GetK() - vCurrent[0];
                    CheckS0(GetK() - vCurrent[0], s0CurrentPrev);
                    PrintData(s0Current, printer, rp, ++iter, m, s0CurrentPrev, vCurrent);
                    s0Diff = Math.Abs(GetK() - vCurrent[0] - s0Current);
                } while (s0Diff > GetS0Eps());

                // Console.WriteLine(m + " " + S0Current);

                // TODO: why we get some negative values instead of zero?
                RemoveNegativeValues(vCurrent);
                
                if (_saveSolutions)
                    SaveNumericSolution(vCurrent, m, s0Current);

                s0T[m] = s0Current;
                for (var i = 0; i < vCurrent.Length; i++)
                    vNext[i] = vCurrent[i];

                if (_allowOutputFile)
                    printer.PrintXY(Path.Combine(_outputPath, "V" + m), GetTau(), m, GetH(), vCurrent, s0Current);

                if (_allowOutputConsole)
                    Console.WriteLine("--------------------------------------------------");
            }

            return s0T;
        }

        private static void CheckHCorrectness(double h, double tau, double sph, double r)
        {
            if (tau / h > 1d / (2d * r * sph))
            {
                // ReSharper disable once NotResolvedInText
                //throw new ArgumentOutOfRangeException("tau/h");
            }
        }

        private static void RemoveNegativeValues(IList<double> vCurrent)
        {
            for (var i = 0; i < vCurrent.Count; i++)
                if (vCurrent[i] < 0d)
                    vCurrent[i] = 0d;
        }

        private double[] CalculateRightPart(double s0, IReadOnlyList<double> vk1, double h, double tau, int m)
        {
            var rp = new double[GetN1()];
            var r = GetR();
            var sigma = GetSquaredSigma();

            // var hph0 = Math.Pow(GetAlpha(0, m) * (0+1) * h, GetBeta(m)); // h_{i+1/2}
            var hph0 = h;
            var hmh0 = hph0; // h_{i-1/2}
            var smh0 = s0 - 0.5d * hmh0; // s_{i-1/2}
            var sph0 = s0 + 0.5d * hph0; // s_{i+1/2}
            CheckHCorrectness(hmh0, tau, sph0, r);
            CheckHCorrectness(hph0, tau, sph0, r);
            // var betam1 = 1d / (8d * tau) * (1d + (2d * tau * r * smh0) / hmh0) * (1d + (2d * tau * r * smh0) / hmh0);
            var beta0 =  1d / (8d * tau) * (3d - 2d * tau * r * smh0 / hmh0) * (1d + 2d * tau * r * smh0 / hmh0) 
                + 1d / (8d * tau) * (3d + 2d * tau * r * sph0 / hph0) * (1d - 2d * tau * r * sph0 / hph0);
            var betap1 = 1d / (8d * tau) * (1d - 2d * tau * r * sph0 / hph0) * (1d - 2d * tau * r * sph0 / hph0);
            var f0 = GetF(sigma, r, GetK(), h, s0, GetTau(), m, GetT(), s0, 0);

            var x0 = -hph0 / 2d * f0 + sigma * s0 * s0 / 2d
                                     + hph0 / 2d * (beta0 * vk1[0] + betap1 * vk1[1]);
            rp[0] = x0;

            for (var i = 1; i < rp.Length - 1; ++i)
            {
                //Console.WriteLine("i={0} beta={1} alpha={2} h={3}", i, GetBeta(m), GetAlpha(m), Math.Pow(GetAlpha(m) * i * h, GetBeta(m)));
                var xi = s0 + i * h;
                // var hmh = xi - (S0 + Math.Pow(GetAlpha(i,m) * (i - 1) * h, GetBeta(m))); // h_{i-1/2}
                // var hph = S0 + Math.Pow(GetAlpha(i,m) * (i + 1) * h, GetBeta(m)) - xi; // h_{i+1/2}
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

                var f = GetF(sigma, r, GetK(), h, xi, tau, m, GetT(), s0, i);
                rp[i] = (hmh + hph) / 2d * f + (hmh + hph) / 2d * (betaIm1 * vk1[i - 1] + betaI * vk1[i] + betaIp1 * vk1[i + 1]);
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

        private double[] CalculateV(double s0Old, IReadOnlyList<double> vk1, double tau, int m, double h, out double[] rp)
        {
            rp = CalculateRightPart(s0Old, vk1, h, tau, m);
            double[] bT = GetB(GetN1(), s0Old, h, GetSquaredSigma(), tau);
            double[] cT = GetC(GetN1(), s0Old, h, GetSquaredSigma(), tau, GetR());
            double[] dT = GetD(GetN1(), s0Old, h, GetSquaredSigma(), tau);
            double[] vk = ThomasAlgorithmCalculator.Calculate(bT, cT, dT, rp);

            // PrintThomasArraysToConsole(b_t, c_t, d_t);
            return vk;
        }

        private void CheckS0(double s0Prev, double s0Current)
        {
            if (s0Prev <= 0d)
                throw new Exception("S0New <= 0d");

            if (s0Prev >= GetK())
                throw new Exception($"S0New >= K: S0New = {s0Prev} K = {GetK()}");

            if (s0Prev >= s0Current)
                throw new Exception($"S0CurrentPrev: {s0Prev.ToString("e8", CultureInfo.InvariantCulture)} >= S0Current: {s0Current.ToString("e8", CultureInfo.InvariantCulture)}");
        }

        private void CreateFileWithStatistics()
        {
            if (!_allowOutputFile)
                return;

            if (File.Exists(Path.Combine(_outputPathStat, "stat.txt")))
                File.Delete(Path.Combine(_outputPathStat, "stat.txt"));
        }

        private double Get_a0(double T,  double k, double smooth)
        {
            var calculateA0 = k / (smooth * T);
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
        private double[] GetB(int n, double s0, double h, double sigmaSq, double tau)
        {
            var b = new double[n];
            for (var i = 1; i <= n - 1; ++i)
            {
                var si = s0 + i * h;
                var hmh = s0 + i * h - (s0 + (i - 1) * h); // h_{i - 1/2}
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
        private double[] GetC(int n, double s0, double h, double sigmaSq, double tau, double r)
        {
            var c = new double[n];
            for (var i = 1; i <= n - 1; ++i)
            {
                var si = s0 + i * h;
                var hmh = si - (s0 + (i - 1) * h); // h_{i-1/2}
                var hph = s0 + (i + 1) * h - si; // h_{i+1/2}
                if (hmh * hmh > 4d * tau * sigmaSq * si * si)
                    throw new ArgumentException("hmh is invalid");

                c[i] = (sigmaSq * si * si) / (2d * hmh) + (sigmaSq * si * si) / (2d * hph) + (hmh + hph) * (1d / (4d * tau) + r / 2d);
            }

            // left boundary condition
            var si0 = s0 + 0 * h;
            var hph0 = s0 + (0 + 1) * h - si0; // h_{i+1/2}
            var hmh0 = si0 - (s0 + (0d - 1d) * h); // h_{i-1/2}
            var sph0 = si0 + 0.5d * hph0; // s_{i+1/2}
            var smh0 = si0 - 0.5d * hmh0; // s_{i-1/2}
            CheckHCorrectness(hph0, tau, sph0, r);
            CheckHCorrectness(hmh0, tau, smh0, r);

            c[0] = sigmaSq * si0 * si0 / (2d * hph0) + hph0 * r / 2d + hph0 / (4d * tau);

            // right boundary condition
            c[n - 1] = 0d;

            return c;
        }

        // диагональ, лежащая над главной (нумеруется: [0;n-2])
        // pass n = n + 1!
        private double[] GetD(int n, double s0, double h, double sigmaSq, double tau)
        {
            var d = new double[n];
            for (var i = 0; i <= n - 2; ++i)
            {
                var si = s0 + i * h;
                var hph = s0 + (i + 1) * h - si; // h_{i+1/2}

                d[i] = hph / (4d * tau) - sigmaSq * si * si / (2d * hph);
            }

            // left boundary condition
            var si0 = s0 + 0 * h;
            var hph0 = s0 + (0 + 1) * h - si0; // h_{i+1/2}

            d[0] = hph0 / (4d * tau) - sigmaSq * si0 * si0 / (2d * hph0);

            // right boundary condition
            d[n - 2] = 0d;
            return d;
        }

        private Point[] GetExactSolution(double tau, int m, double s0)
        {
            double t = tau * m;
            // UpdateH(S0);
            
            var v = new Point[GetN1()];
            var r = GetR();
            var k = GetK();
            var T = GetT();
            var sigma = GetSquaredSigma();

            var list = new List<Point>();
            var superH = GetSuperH();
            int j = 0;
            while (0d + j * superH < s0)
            {
                var s = 0d + j * superH;
                if (s > s0) break;
                list.Add(new Point(s, GetK()-s));
                j++;
            }
            
            for (int i = 0; i < GetN1() ; i++)
            {
                var s = s0 + i * GetH();
                // Vi
                var vs = (T - t)
                         * (sigma / (2d * r))
                         * Math.Pow(k / (1d + sigma / (2d * r)), (2d * r + sigma) / sigma)
                         * Math.Pow(s - Get_a0(T, k, GetSmoothness()) * t, -2d * r / sigma);

                v[i] = new Point(s, vs);
            }

            foreach (var point in v)
                list.Add(point);

            // for (int i = 0; i < GetN1() ; i++)
            // {
            //     // the first version (based on t*(exact_soltion_perp_am_option)) that is not working
            //     // because it is not time-dependent
            //     // var si = S0 + i * GetH();
            //     // var p1 = GetSquaredSigma() / (2d * GetR());
            //     //
            //     // var arg = GetK() / (1 + GetSquaredSigma() / (2d * GetR()));
            //     // var pow = (2d * GetR() + GetSquaredSigma()) / GetSquaredSigma();
            //     // var p2 = Math.Pow(arg, pow);
            //     // var arg2 = si;
            //     // var pow2 = -2d * GetR() / GetSquaredSigma();
            //     // var p3 = Math.Pow(arg2, pow2);
            //     //
            //     // var v = p1 * p2 * p3;
            //     // V[i] = tk * v;
            //     
            //     // the second version
            //     // the second version (check validation-equation-am_option.nb and validation-equation-am_option.docx)
            //     // missed
            //     var S = S0 + i * GetH();
            //
            //     // Vi
            //     var VS = (T - t)
            //              * (sigma / (2d * r))
            //              * Math.Pow(K / (1d + sigma / (2d * r)), (2d * r + sigma) / sigma)
            //              * Math.Pow(S - Get_a0(T, K, GetSmoothness()) * t, -2d * r / sigma);
            //              // * Math.Pow(S - S0 + (K - S0) * ((2d * r) / (sigma_sq)), (-2d * r) / sigma_sq);
            //     
            //     // the third version 
            //     // var S = S0 + i * GetH();
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
            // UpdateH(s0);

            var h = GetH();
            var v = new List<Point>();
            var r = GetR();
            var k = GetK();
            var T = GetT();
            var sigma = GetSquaredSigma();
            var i = 0;
            while (i * GetH() < s0)
            {
                v.Add(new Point(
                    Geta() + i * GetH(), 
                    GetK() - i * GetH()));
                i++;
            }

            startPosOfS0 = i;
            i = 0;
            while (s0 + i * h < GetK())
            {
                var s = s0 + i * h;
                var vs = (T - t)
                        * (sigma / (2d * r))
                        * Math.Pow(k / (1d + sigma / (2d * r)), (2d * r + sigma) / sigma)
                        * Math.Pow(s - Get_a0(T, k, _smoothness) * t, -2d * r / sigma);
                v.Add(new Point(s, vs));
                
                i++;
            }
            
            while (s0 + i * h < Getb())
            {
                v.Add(new Point());
                i++;
            }

            return v.ToArray();
        }

        [SuppressMessage("ReSharper", "UnusedParameter.Local")]
        private double GetF(double sigma, double r, double k, double h, double s,
             double tau, int m, double T, double s0, int i)
        {
            return 0d;
            
            //// var t = tau * k; // => k = (M-m) -> t=(M-m)*tau
            
            //// dynamic t mesh
            //var t = Math.Pow(GetAlpha(i,m) * m * tau, GetBeta(m)); // => k = (M-m) -> t=(alpha*(M-m)*tau)^beta

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
            
            //var a0 = Get_a0(T, K, GetSmoothness());
            
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
            return _smoothness;
        }

        private double[] GetVst()
        {
            // UpdateH(GetK());
            var arr = new double[GetN1()];
            for (var i = 0; i < arr.Length; ++i)
            {
                var s = GetK() + i * GetH();
                if (s <= GetK())
                    arr[i] = GetK() - s;
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

        private void PrintData(double s0Current, TecplotPrinter tecplotPrinter, double[] rp, int iter, int m, double s0Next, IReadOnlyList<double> vk)
        {
            PrintRpToTecplot(tecplotPrinter, rp, s0Current);
            PrintStatistics(iter, m, GetH(), s0Next, s0Current);
            PrintValuesToConsole(iter, GetH(), s0Next, vk, s0Current);
        }

        private void PrintHeader(IReadOnlyList<double> st)
        {
            if (!_allowOutputConsole)
                return;

            Console.WriteLine($"Number of time steps = {GetM()} h = {GetH()} S(T) = {st[st.Count - 1]}");
            Console.WriteLine("--------------------------------------------------");
        }

        private void PrintRpToTecplot(TecplotPrinter tecplotPrinter, double[] rp, double s0New)
        {
            if (_allowOutputFile)
                tecplotPrinter.PrintXY(Path.Combine(_outputPathRp, "temporal-rp"), 0d, 0, GetH(), rp, s0New);
        }

        private void PrintStatistics(int iter, int m, double hOld, double s0Old, double s0New)
        {
            if (!_allowOutputFile)
                return;

            using var streamWriter = File.AppendText(Path.Combine(_outputPathStat, "stat.txt"));
            streamWriter.WriteLine(new string(' ', 2) + " Time step = {0} Iteration = " + iter + " h = {1} S0 = {2} Abs(S0New-S0Old)={3} S0Eps={4} Cnd={5}", m, hOld, s0Old,
                Math.Abs(s0New - s0Old),
                GetS0Eps(),
                Math.Abs(s0New - s0Old) > GetS0Eps());
        }

        private void PrintValuesToConsole(int iter, double hOld, double s0Old, IReadOnlyList<double> vk, double s0New)
        {
            if (_allowOutputConsole)
                Console.WriteLine(new string(' ', 2) + "Iteration = " + iter + " h = {0} S0_old = {1} Vk[0] = {2} S0_new = {3}", hOld, s0Old, vk[0], s0New);
        }

        private void PrintVst(TecplotPrinter tecplotPrinter, double[] vk1, IReadOnlyList<double> st)
        {
            if (_allowOutputFile)
                tecplotPrinter.PrintXY(_outputPath + "V" + GetM(), GetTau(), GetM(), GetH(), vk1, st[st.Count - 1]);
        }

        private void SaveNumericSolution(IReadOnlyList<double> v, int m, double s0)
        {
            var temp = new Point[v.Count];
            for (var i = 0; i < v.Count; i++)
                temp[i] = new Point(s0+i*GetH(),v[i]);

            _solutions.Add(new SolutionData(s0, m)
            {
                Solution = temp
            });
        }
    }
}