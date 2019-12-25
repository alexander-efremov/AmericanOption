namespace PerpetualAmericanOptions
{
    using System;
    using System.Collections.Generic;
    using System.Diagnostics.CodeAnalysis;
    using System.IO;

    using CoreLib;

    using TemporalAmericanOption;

    // from new presentation with FEM
    [SuppressMessage("ReSharper", "CommentTypo", Justification = "OK.")]
    public class TemporalAmericanOptionCalculator : AmericanOptionCalculator
    {
        private readonly bool allowOutputFile;

        private readonly bool allowOutputConsole;

        private readonly string outputPath;

        private readonly string outputPathStat;

        private readonly string outputPathRp;

        private readonly List<double[]> solutions = new List<double[]>();

        private readonly bool saveSolutions;

        private readonly int M;

        private readonly double T;

        public TemporalAmericanOptionCalculator(TemporalParameters parameters, bool allowOutputFile, bool allowOutputConsole)
            : base(parameters)
        {
            this.allowOutputFile = allowOutputFile;
            this.allowOutputConsole = allowOutputConsole;
            this.outputPath = parameters.WorkDir;
            if (string.IsNullOrEmpty(this.outputPath))
            {
                this.outputPath = this.GetWorkDir() + "AO/";
            }

            this.M = parameters.M;
            this.T = parameters.T;
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
            return this.M;
        }

        public double GetT()
        {
            return this.T;
        }

        public List<double[]> GetNumericSolutions()
        {
            return this.solutions;
        }

        public double[] Solve()
        {
            var tecplotPrinter = new TecplotPrinterSpecial(this.GetN1(), 0d, this.GetRightBoundary(), this.GetTau());

            // var printer = new ThomasArrayPrinter();
            var Vk1 = this.GetVST(); // V(S, T) = (K - S)+

            // tecplotPrinter.PrintXY(Path.Combine(_outputPath, "VST"), 0d, GetH(), Vk1);
            this.PushToSolutions(Vk1);
            var St = new double[this.GetM()];
            St[St.Length - 1] = this.GetK();

            this.PrintHeader(St);
            this.CreateStatFile();
            this.PrintVST(tecplotPrinter, Vk1, St);

            for (int k = this.GetM() - 2; k >= 0; --k)
            {
                var ts = k + 1;
                if (this.allowOutputConsole)
                {
                    Console.WriteLine("Time step = " + ts);
                }

                double[] Vk;
                int iter = 0;
                double S0Old;
                double S0New = St[St.Length - (this.GetM() - ts)];
                do
                {
                    S0Old = S0New;
                    this.UpdateH(S0Old);
                    double[] rp = this.CalculateVk(S0Old, Vk1, this.GetTau(), out Vk);
                    S0New = this.GetK() - Vk[0];

                    this.PrintRpToTecplot(tecplotPrinter, rp, S0New);
                    iter++;
                    this.PrintStatistics(iter, k, this.GetH(), S0Old, S0New);
                    this.PrintValuesToConsole(iter, this.GetH(), S0Old, Vk, S0New);
                    this.CheckS0NewValidity(S0New);
                    if (this.allowOutputConsole)
                    {
                        Console.WriteLine(
                            new string(' ', 2) + "Iteration = " + iter + " h = {0} S0_old = {1} Vk[0] = {2} S0_new = {3}",
                            this.GetH(),
                            S0Old,
                            Vk[0],
                            S0New);
                    }

                    if (S0New <= 0d) 
                    {
                        throw new Exception("S0New <= 0d");
                    }

                    if (S0New >= this.GetK())
                    {
                        throw new Exception("S0New >= GetK()");
                    }
                }
                while (Math.Abs(S0New - S0Old) > this.GetS0Eps());

                this.PushToSolutions(Vk);
                St[St.Length - (this.GetM() - ts) - 1] = S0New;
                for (int i = 0; i < Vk.Length; i++)
                {
                    Vk1[i] = Vk[i];
                }

                if (this.allowOutputFile)
                {
                    tecplotPrinter.PrintXY(Path.Combine(this.outputPath, "V" + ts), this.GetTau() * ts, this.GetH(), Vk, S0New);
                }

                if (this.allowOutputConsole)
                {
                    Console.WriteLine("--------------------------------------------------");
                }
            }

            return St;
        }

        public TecplotPrinterSpecial GetTecplotPrinter()
        {
            var tecplotPrinter = new TecplotPrinterSpecial(this.GetN1(), 0d, this.GetRightBoundary(), this.GetTau());
            return tecplotPrinter;
        }

        

        // private void PrintXY(string filename, double h, double[] data, double start, int tl)
        // {
        // var name = string.Format("{0}_t={1}.dat", filename, tl);
        // using (var writer = new StreamWriter(name, false))
        // {
        // writer.WriteLine("TITLE = \"V(S): from M to 1\"");
        // writer.WriteLine("VARIABLES = S V");
        // writer.WriteLine("DATASETAUXDATA N=\"{0}\"", GetN());
        // writer.WriteLine("DATASETAUXDATA M=\"{0}\"", GetM());
        // writer.WriteLine("DATASETAUXDATA tau=\"{0}\"", GetTau());
        // writer.WriteLine("DATASETAUXDATA K=\"{0}\"", GetK());
        // writer.WriteLine("DATASETAUXDATA r=\"{0}\"", GetR());
        // writer.WriteLine("DATASETAUXDATA sigma_sq=\"{0}\"", GetSquaredSigma());
        // writer.WriteLine("DATASETAUXDATA s0_eps=\"{0}\"", GetS0Eps());
        // writer.WriteLine("DATASETAUXDATA tl=\"{0}\"", tl);
        // writer.WriteLine("ZONE T='SubZone'");
        // writer.WriteLine("I={0} K={1} ZONETYPE=Ordered", data.Length, 1);
        // writer.WriteLine("DATAPACKING=POINT");
        // writer.WriteLine("DT=(DOUBLE DOUBLE)");
        // for (var i = 0; i < data.Length; i++)
        // {
        // writer.WriteLine("{0:e8} {1:e8}", start + i * h, data[i]);
        // }
        // }
        // }
        public List<double[]> GetExactSolutions(double[] s0)
        {
            var list = new List<double[]>();
            if (Math.Abs(this.GetM() * this.GetTau() - this.GetT()) > double.Epsilon)
            {
                throw new Exception("GetExactSolutions");
            }

            list.Add(this.GetVST()); // V(S, T) = (K - S)+
            for (int k = this.GetM() - 1; k >= 1; k--)
            {
                double[] solution = this.GetExactSolution(k, s0);
                var tecplotPrinter = this.GetTecplotPrinter();
                tecplotPrinter.PrintXY(Path.Combine(this.outputPath, "solution"), 0d, this.GetH(), solution);
                list.Add(solution);
            }

            return list;
        }

        private void PrintVST(TecplotPrinterSpecial tecplotPrinter, double[] Vk1, double[] St)
        {
            if (this.allowOutputFile)
            {
                tecplotPrinter.PrintXY(this.outputPath + "V" + this.GetM(), this.GetTau() * this.GetM(), this.GetH(), Vk1, St[St.Length - 1]);
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

        private void PrintHeader(double[] St)
        {
            if (this.allowOutputConsole)
            {
                Console.WriteLine("Time step = " + this.GetM() + string.Format(" h = {0} S(T) = {1}", this.GetH(), St[St.Length - 1]));
                Console.WriteLine("--------------------------------------------------");
            }
        }

        private void CreateStatFile()
        {
            if (this.allowOutputFile)
            {
                if (File.Exists(Path.Combine(this.outputPathStat, "stat.txt")))
                {
                    File.Delete(Path.Combine(this.outputPathStat, "stat.txt"));
                }
            }
        }

        private void PushToSolutions(double[] Vk1)
        {
            if (!this.saveSolutions)
            {
                return;
            }

            var r = new double[Vk1.Length];
            for (int i = 0; i < Vk1.Length; i++)
            {
                r[i] = Vk1[i];
            }

            this.solutions.Add(r);
        }

        private void CheckS0NewValidity(double s0New)
        {
            if (s0New <= 0d)
            {
                throw new Exception("S0New <= 0d");
            }

            if (s0New >= this.GetK())
            {
                throw new Exception("S0New >= GetK()");
            }
        }

        private double[] CalculateVk(double s0Old, double[] Vk1, double tau1, out double[] Vk)
        {
            var h_value = this.GetH();
            var rp = this.CalculateRightPart(s0Old, Vk1, h_value, tau1);
            var b_t = this.GetB(this.GetN1(), s0Old, h_value, this.GetSquaredSigma(), tau1);
            var c_t = this.GetC(this.GetN1(), s0Old, h_value, this.GetSquaredSigma(), tau1, this.GetR());
            var d_t = this.GetD(this.GetN1(), s0Old, h_value, this.GetSquaredSigma(), tau1);
            Vk = this.ThomasAlgorithmCalculator.Calculate(b_t, c_t, d_t, rp);

            // PrintThomasArraysToConsole(b_t, c_t, d_t);
            return rp;
        }

        private void PrintValuesToConsole(int iter, double h_old, double S0Old, double[] Vk, double S0New)
        {
            if (this.allowOutputConsole)
            {
                Console.WriteLine(new string(' ', 2) + "Iteration = " + iter + " h = {0} S0_old = {1} Vk[0] = {2} S0_new = {3}", h_old, S0Old, Vk[0], S0New);
            }
        }

        private void PrintStatistics(int iter, int k, double h_old, double S0Old, double S0New)
        {
            if (this.allowOutputFile)
            {
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
        }

        private void PrintRpToTecplot(TecplotPrinterSpecial tecplotPrinter, double[] rp, double S0New)
        {
            if (this.allowOutputFile)
            {
                tecplotPrinter.PrintXY(Path.Combine(this.outputPathRp, "temporal-rp"), 0d, this.GetH(), rp, S0New);
            }
        }

        private double[] CalculateRightPart(double S0, double[] VK1, double h, double tau)
        {
            var rp = new double[this.GetN1()];

            for (var i = 1; i < rp.Length - 1; ++i)
            {
                var si = S0 + i * h;
                var hmh = si - (S0 + (i - 1) * h); // h_{i-1/2}
                var hph = S0 + (i + 1) * h - si; // h_{i+1/2}
                var smh = si - 0.5d * hmh; // s_{i-1/2}
                var sph = si + 0.5d * hph; // s_{i+1/2}
                CheckHCorrectness(hph, tau, sph, this.GetR());
                CheckHCorrectness(hmh, tau, sph, this.GetR());
                var beta1 = 1d / (8d * tau) * (1d + 2d * tau * this.GetR() * smh / hmh) * (1d + 2d * tau * this.GetR() * smh / hmh);

                var beta2 = 1d / (8d * tau) * (3d - 2d * tau * this.GetR() * smh / hmh) * (1d + 2d * tau * this.GetR() * smh / hmh)
                            + 1d / (8d * tau) * (3d + 2d * tau * this.GetR() * sph / hph) * (1d - 2d * tau * this.GetR() * sph / hph);

                var beta3 = 1d / (8d * tau) * (1d - 2d * tau * this.GetR() * sph / hph) * (1d - 2d * tau * this.GetR() * sph / hph);

                var f = this.GetF(this.GetSquaredSigma(), this.GetR(), this.GetK(), i, S0, h);
                rp[i] = ((hmh + hph) / 2d) * f + ((hmh + hph) / 2d) * (beta1 * VK1[i - 1] + beta2 * VK1[i] + beta3 * VK1[i + 1]);
            }

            var si0 = S0 + 0 * h;
            var hph0 = S0 + (0 + 1) * h - si0; // h_{i+1/2}
            var sph0 = si0 + 0.5d * hph0; // s_{i+1/2}
            CheckHCorrectness(hph0, tau, sph0, this.GetR());
            var beta20 = 1d / (8d * tau) * (3d + 2d * tau * this.GetR() * sph0 / hph0) * (1d - 2d * tau * this.GetR() * sph0 / hph0);
            var beta30 = 1d / (8d * tau) * (1d - (2d * tau * this.GetR() * sph0) / hph0) * (1d - (2d * tau * this.GetR() * sph0) / hph0);
            var f0 = this.GetF(this.GetSquaredSigma(), this.GetR(), this.GetK(), 0, S0, h);

            rp[0] = ((-hph0 / 2d) * f0) + ((this.GetSquaredSigma() * si0 * si0) / 2d) + hph0 * (beta20 * VK1[0] + beta30 * VK1[1]);

            rp[rp.Length - 1] = 0d;

            return rp;
        }

        private double GetF(double sigma_sq, double r, double K, int i, double S0, double h)
        {
            // return 0;
            var si = S0 + i * h;
            var p1 = sigma_sq / (2d * r);

            var arg = K / (1 + sigma_sq / (2d * r));
            var pow = (2d * r + sigma_sq) / sigma_sq;
            var p2 = Math.Pow(arg, pow);

            var arg2 = si;
            var pow2 = -2d * r / sigma_sq;
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
        private double[] GetB(int n, double S0, double h, double sigmaSq, double tau)
        {
            var b = new double[n];
            for (var i = 1; i <= n - 1; ++i)
            {
                var si = S0 + i * h;
                var hmh = (S0 + i * h) - (S0 + (i - 1) * h); // h_{i - 1/2}
                if (hmh * hmh > 2d * tau * sigmaSq * si * si)
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
                if (hmh * hmh > 2d * tau * sigma_sq * si * si)
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
                // throw new ArgumentOutOfRangeException("tau/h");
            }
        }

        private double[] GetExactSolution(int tl, double[] s0)
        {
            var V = new double[this.GetN1()];
            for (int i = 0; i < V.Length; i++)
            {
                var si = s0[this.GetM() - 1 - tl] + i * this.GetH();
                var p1 = this.GetSquaredSigma() / (2d * this.GetR());

                var arg = this.GetK() / (1 + this.GetSquaredSigma() / (2d * this.GetR()));
                var pow = (2d * this.GetR() + this.GetSquaredSigma()) / this.GetSquaredSigma();
                var p2 = Math.Pow(arg, pow);

                var arg2 = si;
                var pow2 = -2d * this.GetR() / this.GetSquaredSigma();
                var p3 = Math.Pow(arg2, pow2);

                var v = p1 * p2 * p3;
                var t = tl * this.GetTau();
                V[i] = t * v;
            }

            return V;
        }
    }
}