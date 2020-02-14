namespace PerpetualAmericanOptions
{
    using System;
    using System.Collections.Generic;
    using System.Diagnostics;
    using System.Globalization;
    using System.IO;
    using NUnit.Framework;
    using TemporalAmericanOption;

    public class TemporalAmericanOptionTests : UnitTestBase
    {
        [Test]
        public void TemporalAmericanOption()
        {
            var parameters = this.GetParameters(true, this.GetWorkingDir() + "AO/");
            var calculator = new TemporalAmericanOptionCalculator(parameters, true, true);

            PrintParameters(calculator);

            var s0 = calculator.Solve();

            Console.WriteLine("Numeric S0:");
            for (var i = s0.Length - 1; i >= 0; i--)
            {
                Console.WriteLine("k = " + (i + 1) + " -> " + s0[i]);
            }

            var exactSolutions = calculator.GetExactSolutions(s0);
            var numericSolutions = calculator.GetNumericSolutions();
            
            for (var i = 0; i < numericSolutions.Count; i++)
            {
                var exactSolution = exactSolutions[i];
                var numericSolution = numericSolutions[i];
                var error = Utils.GetAbsError(exactSolution, numericSolution);
                
                var printer = calculator.GetTecplotPrinter();
                printer.PrintXY(Path.Combine(parameters.WorkDir, "exactSolution"), 0d, calculator.GetH(), exactSolution);
                printer.PrintXY(Path.Combine(parameters.WorkDir, "numericSolution"), 0d, calculator.GetH(), numericSolution);
                printer.PrintXY(Path.Combine(parameters.WorkDir, "error"), 0d, calculator.GetH(), error);
            }
            
            Assert.AreEqual(exactSolutions.Count, numericSolutions.Count);
            for (var i = 0; i < numericSolutions.Count; i++)
            {
                var exactSolution = exactSolutions[i];
                var numericSolution = numericSolutions[i];
                var error = Utils.GetAbsError(exactSolution, numericSolution);
                
                Assert.AreEqual(exactSolution.Length, numericSolution.Length);
                for (var j = 0; j < exactSolution.Length; j++)
                {
                    Assert.AreEqual(exactSolution[j], numericSolution[j], error[j].ToString(CultureInfo.InvariantCulture));
                }
            }

            this.Print(Path.Combine(parameters.WorkDir, "s0"), s0, calculator.GetTau());
        }

        [Test]
        public void TestSeries()
        {
            var startN = 150;
            var Ksteps = 3;
            var Nsteps = 6;
            var tls = 3;
            var startK = 5d;
            var a = 0d;
            var b = 50d;
            var r = 0.1d;
            var sigmaSq = 0.2d;
            var K = startK;
            var M = 20;
            var S0eps = 1e-5;
            var T = 1d;
            var tau = (T - 0d) / M;
            this.PrintParamsForSeries(r, sigmaSq, K, M, tau);
            PrintTableHeader(Nsteps, startN);

            var table = new Dictionary<string, List<double>>();

            for (double ki = 1; ki <= Ksteps; ki++)
            {
                var Ki = K * ki;
                for (var i = 0; i < Nsteps; i++)
                {
                    var n = (int)Math.Pow(2, i) * startN;
                    var folderPath = this.CreateOutputFolder(Ki, n, string.Empty);
                    var parameters = this.GetSeriesParameters(n, T, Ki, M, tau, a, b, r, sigmaSq, S0eps, folderPath);
                    var calculator = new TemporalAmericanOptionCalculator(parameters, false, false);
                    double[] S0Arr = calculator.Solve();

                    for (var tl = tls; tl >= 1; tl--)
                    {
                        var time = tau * (5d * tl);
                        var timeIndex = (int)(time / tau);
                        var value = S0Arr[timeIndex];
                        var s1 = value;
                        var format = Ki + new string(' ', Ki > 5 ? 9 : 10) + time.ToString("0.00") +
                                     new string(' ', 10);
                        if (table.ContainsKey(format))
                        {
                            table[format].Add(s1);
                        }
                        else
                        {
                            table[format] = new List<double> {s1};
                        }
                    }
                }
            }

            PrintTable(table);
            Console.WriteLine();
            this.PrintCsv(table);
        }

        [Test]
        public void TestSeriesConvergence()
        {
            var startN = 150;

            // int Ksteps = 3;
            var Nsteps = 7;

            // int tls = 3;
            var startK = 5d;
            var a = 0d;
            var b = 50d;
            var r = 0.1d;
            var sigmaSq = 0.2d;
            var K = startK;
            var S0eps = 1e-5;
            var T = 1d;
            this.PrintParamsForSeries2(r, sigmaSq, K);
            PrintTableHeader(Nsteps, startN);

            var w = new Stopwatch();
            w.Reset();
            int M;
            double tau;
            double[] S0ArrGold;
            Console.WriteLine("Started Nsteps - 1");
            w.Reset();
            w.Start();
            {
                M = (int)(10 * Math.Pow(4, Nsteps - 1));
                tau = (T - 0d) / M;
                var n = (int)Math.Pow(2, Nsteps - 1) * startN;
                var folderPath = this.GetTrashFolder();

                /*
                                                if (allowOutputFile)
                                                {
                                                    folderPath = CreateOutputFolder(K, n, "convergence");
                                                }
                                */
                var parameters = this.GetSeriesParameters(n, T, K, M, tau, a, b, r, sigmaSq, S0eps, folderPath);
                var calculator = new TemporalAmericanOptionCalculator(parameters, false, false);
                S0ArrGold = calculator.Solve();
            }

            Console.WriteLine("Finished Nsteps - 1");
            var dic = new Dictionary<int, double>();
            Console.WriteLine("Started from 0 to Nsteps - 1");
            for (var i = 0; i < Nsteps - 1; i++)
            {
                Console.WriteLine("Started step " + i);
                M = (int)(10 * Math.Pow(4, i));
                tau = (T - 0d) / M;
                var n = (int)Math.Pow(2, i) * startN;
                var folderPath = this.GetTrashFolder();

                /*
                                                if (allowOutputFile)
                                                {
                                                    folderPath = CreateOutputFolder(K, n, "convergence");
                                                }
                                */
                var parameters = this.GetSeriesParameters(n, T, K, M, tau, a, b, r, sigmaSq, S0eps, folderPath);
                var calculator = new TemporalAmericanOptionCalculator(parameters, false, false);
                double[] S0Arr = calculator.Solve();
                dic[n] = this.GetErrorLInf(S0ArrGold, S0Arr);
                Console.WriteLine("Finished step = " + i);
            }

            w.Stop();
            Console.WriteLine("Elapsed = " + w.Elapsed.TotalSeconds + " s.");
            var l = new List<double>();
            foreach (KeyValuePair<int, double> pair in dic)
            {
                l.Add(pair.Value);
                Console.WriteLine(pair.Key + " " + pair.Value);
            }

            Console.WriteLine();
            Console.WriteLine("Log2");
            for (var i = 1; i < l.Count; i++)
            {
                var value = l[i - 1] / l[i];
                Console.WriteLine(value + " " + Math.Log(value, 2));
            }

            // PrintTable(table);
            // Console.WriteLine();
            // PrintCsv(table);
        }

        private string GetWorkingDir()
        {
            return Environment.GetFolderPath(Environment.SpecialFolder.Desktop) + Path.DirectorySeparatorChar;
        }

        private double GetErrorLInf(double[] gold, double[] sol)
        {
            var stride = (int)((double)gold.Length / sol.Length);
            var nSol = new double[gold.Length];
            var error = new double[gold.Length];
            for (int i = 0, k = 0; i < gold.Length; i += stride, k++)
            {
                nSol[i] = sol[k];
                error[i] = Math.Abs(gold[i] - nSol[i]);
            }

            // var error = Utils.GetAbsError(gold, nSol);
            return Utils.GetLInf(error);
        }

        private string GetTrashFolder()
        {
            return this.GetWorkingDir() + "AO/trash";
        }

        private void PrintCsv(Dictionary<string, List<double>> table)
        {
            var dictionary = new Dictionary<string, string>();
            foreach (KeyValuePair<string, List<double>> pair in table)
            {
                string[] strings = pair.Key.Split(new[] {' '}, StringSplitOptions.RemoveEmptyEntries);
                Tuple<string, string> tuple = Tuple.Create(strings[0], strings[1]);
                var tupleItem1 = tuple.Item1 + ";" + tuple.Item2;
                foreach (var tt1 in pair.Value)
                {
                    tupleItem1 += ";" + tt1.ToString("0.0000000");
                }

                dictionary[pair.Key] = tupleItem1;
            }

            foreach (KeyValuePair<string, string> pair in dictionary)
            {
                Console.WriteLine(pair.Value);
            }
        }

        private void PrintParamsForSeries(double r, double sigmaSq, double K, int M, double tau)
        {
            Console.Write("r = " + r);
            Console.Write(" M = " + M);
            Console.Write(" tau = " + tau);
            Console.Write(" sigma_sq = " + sigmaSq);
            Console.WriteLine(" K = " + K);
        }

        private void PrintParamsForSeries2(double r, double sigmaSq, double K)
        {
            Console.Write("r = " + r);
            Console.Write(" sigma_sq = " + sigmaSq);
            Console.WriteLine(" K = " + K);
        }

        private string CreateOutputFolder(double Ki, int n, string subfolder)
        {
            var folderPath = this.GetWorkingDir() + "/AO/" + subfolder + "/" + Ki + "_" + "_" + n + "/";
            if (!Directory.Exists(folderPath))
            {
                Directory.CreateDirectory(folderPath);
            }

            return folderPath;
        }

        private static void PrintTableHeader(int Nsteps, int startN)
        {
            var whitespaceNumber = 25 * Nsteps;
            Console.WriteLine(new string('=', whitespaceNumber));
            Console.Write("K" + new string(' ', 10));
            Console.Write("t" + new string(' ', 13));
            for (var i = 0; i < Nsteps; i++)
            {
                var n = (int)Math.Pow(2, i) * startN;
                Console.Write("N(" + n + ")" + new string(' ', 15));
            }

            Console.WriteLine();
            Console.WriteLine(new string('=', whitespaceNumber));
        }

        private static void PrintTable(Dictionary<string, List<double>> table)
        {
            foreach (var pair in table)
            {
                var s = pair.Key;
                foreach (var t in pair.Value)
                {
                    s += t.ToString("0.0000000") + new string(' ', t > 10d ? 10 : 11);
                }

                Console.WriteLine(s);
            }
        }

        private void Print(string filename, double[] St, double tau)
        {
            var name = $"{filename}_nx={St.Length}_tau={tau}.dat";
            using (var writer = new StreamWriter(name, false))
            {
                writer.WriteLine("TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA'\nVARIABLES = S {0}\nZONE T='{1}'", "t", "SubZone");
                writer.WriteLine($"I={St.Length} K={1} ZONETYPE=Ordered");
                writer.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE)");
                for (var i = 0; i < St.Length; i++)
                {
                    writer.WriteLine("{0:e8} {1:e8}", St[i], i);
                }
            }
        }

        private TemporalParameters GetParameters(bool saveSolutions, string workPath)
        {
            var a = 0d;
            var b = 10d;
            var r = 0.1d;
            var n = 300;
            var tau = 2e-5;
            var sigmaSq = 0.2d;
            var eps = 1e-5;
            var K = 5d;
            var M = 2;
            var T = M * tau;
            return new TemporalParameters(a, b, n, r, tau, sigmaSq, K, eps, M, T, workPath)
                       {
                           SaveVSolutions = saveSolutions
                       };
        }

        private TemporalParameters GetSeriesParameters(
            int n,
            double T,
            double K,
            int M,
            double tau,
            double a,
            double b,
            double r,
            double sigma,
            double S0eps,
            string workDir)
        {
            return new TemporalParameters(a, b, n, r, tau, sigma, K, S0eps, M, T, workDir);
        }

        private static void PrintParameters(TemporalAmericanOptionCalculator calculator)
        {
            Console.WriteLine("a = " + calculator.GetLeftBoundary());
            Console.WriteLine("b = " + calculator.GetRightBoundary());
            Console.WriteLine("r = " + calculator.GetR());
            Console.WriteLine("N = " + calculator.GetN());
            Console.WriteLine("N_1 = " + calculator.GetN1());
            Console.WriteLine("tau = " + calculator.GetTau());
            Console.WriteLine("sigma_sq = " + calculator.GetSquaredSigma());
            Console.WriteLine("K = " + calculator.GetK());
            Console.WriteLine("M = " + calculator.GetM());
            Console.WriteLine("T = " + calculator.GetT());
            Console.WriteLine("S0 Eps = " + calculator.GetS0Eps());
            Console.WriteLine();
        }
    }
}