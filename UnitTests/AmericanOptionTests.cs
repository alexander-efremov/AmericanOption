namespace UnitTests
{
    using System;
    using System.Collections.Generic;
    using System.Diagnostics;
    using System.Globalization;
    using System.IO;
    using System.Linq;
    using System.Text;
    using NUnit.Framework;
    using AmericanOption;
    using AmericanOptions;
    using CoreLib;
    using CoreLib.Utils;

    [TestFixture]
    public class AmericanOptionTests : UnitTestBase
    {
        [Test]
        public void AmericanOption()
        {
            var parameters = this.GetParameters(true, this.GetWorkingDir() + "AO/");
            var calculator = new AmericanOptionCalculator(parameters, true, true);

            PrintParameters(calculator);

            // get S0t in direct order 
            double[] S0t = calculator.Solve();
            // get V(S(t), S0(t)) in direct order 
            List<SolutionData> numericSolutions = calculator.GetNumericSolutions();
            
            Console.WriteLine("Numeric S0:");
            for (var i = S0t.Length - 1; i >= 0; i--)
            {
                Console.WriteLine("k = " + (i + 1) + " -> " + S0t[i]);
            }

            List<SolutionData> exactSolutions = calculator.GetExactSolutions(S0t);

            Assert.AreEqual(exactSolutions.Count, numericSolutions.Count);

            Console.WriteLine();
            Console.WriteLine("Print Exact Sols");
            for (var index = 0; index < exactSolutions.Count; index++)
            {
                var exactSolution = exactSolutions[index];
                var sb = new StringBuilder();
                foreach (var d in exactSolution.Solution.Take(5))
                {
                    sb.Append(d.ToString("e8", CultureInfo.InvariantCulture) + " ");
                }

                Console.WriteLine(exactSolution.k + " " + exactSolution.S0.ToString("e8", CultureInfo.InvariantCulture) + " " + sb);
                if (index > 10)
                {
                    break;
                }
            }

            Console.WriteLine();
            Console.WriteLine("Print Numeric Sols");
            for (var index = 0; index < numericSolutions.Count; index++)
            {
                var numericSolution = numericSolutions[index];
                var sb = new StringBuilder();
                foreach (var d in numericSolution.Solution.Take(5))
                {
                    sb.Append(d.ToString("e8", CultureInfo.InvariantCulture) + " ");
                }

                Console.WriteLine(numericSolution.k + " " + numericSolution.S0.ToString("e8", CultureInfo.InvariantCulture) + " " + sb);
                if (index > 10)
                {
                    break;
                }
            }
            
            Console.WriteLine();
            Console.WriteLine("Print Errors");
            for (var index = 0; index < numericSolutions.Count; index++)
            {
                var exactSolution = exactSolutions[index];
                var numericSolution = numericSolutions[index];
                double[] error = Utils.GetAbsError(exactSolution.Solution, numericSolution.Solution);
                var sb = new StringBuilder();
                foreach (var d in error.Take(5))
                {
                    sb.Append(d.ToString("e8", CultureInfo.InvariantCulture) + " ");
                }

                Console.WriteLine(numericSolution.k + " " + numericSolution.S0.ToString("e8", CultureInfo.InvariantCulture) + " " + sb);
                if (index > 10)
                {
                    break;
                }
            }

            ClearData(parameters);

            var printer = calculator.GetTecplotPrinter();
            var solutionsCount = numericSolutions.Count;
            for (var k = 0; k < solutionsCount; k++)
            {
                var exactSolution = exactSolutions[k];
                var numericSolution = numericSolutions[k];
                double[] error = Utils.GetAbsError(exactSolution.Solution, numericSolution.Solution);

                printer.PrintXY(Path.Combine(parameters.WorkDir , "exactSolution"), calculator.GetTau() * k, calculator.GetH(),
                    exactSolution.Solution, k, "exact_" + k);
                printer.PrintXY(Path.Combine(parameters.WorkDir, "numericSolution"), calculator.GetTau() * k, calculator.GetH(),
                    numericSolution.Solution, k, "numeric_" + k);
                printer.PrintXY(Path.Combine(parameters.WorkDir, "error"), calculator.GetTau() * k, calculator.GetH(), error);
            }
            
            for (var i = 0; i < solutionsCount; i++)
            {
                var exactSolution = exactSolutions[i];
                var numericSolution = numericSolutions[i];
                Assert.AreEqual(exactSolution.Solution.Length, numericSolution.Solution.Length);
                
                double[] error = Utils.GetAbsError(exactSolution.Solution, numericSolution.Solution);
                
                for (var j = 0; j < exactSolution.Solution.Length; j++)
                {
                    try
                    {
                        Assert.AreEqual(
                            exactSolution.k, 
                            numericSolution.k,
                            "k is wrong");
                        Assert.AreEqual(
                            exactSolution.S0, 
                            numericSolution.S0,
                            "S0 is wrong");
                        Assert.AreEqual(
                            exactSolution.Solution[j], 
                            numericSolution.Solution[j],
                            // 10e-6,
                            error[j].ToString(CultureInfo.InvariantCulture));
                    }
                    catch (Exception)
                    {
                        Console.WriteLine("i = " + i + " j = " + j);
                        throw;
                    }
                }
            }

            this.Print(Path.Combine(parameters.WorkDir, "s0"), S0t, calculator.GetTau());
        }

        private static void ClearData(Parameters parameters)
        {
            var di = new DirectoryInfo(parameters.WorkDir);
            var files = di.GetFiles("*.dat")
                .Where(p => p.Extension == ".dat").ToArray();
            foreach (var file in files)
            {
                try
                {
                    file.Attributes = FileAttributes.Normal;
                    File.Delete(file.FullName);
                }
                catch
                {
                    // ignored
                }
            }
        }

        [Test]
        public void TestSeries()
        {
            const int startN = 150;
            const int Ksteps = 3;
            const int Nsteps = 6;
            const int tls = 3;
            const double startK = 5d;
            const double a = 0d;
            const double b = 50d;
            const double r = 0.1d;
            const double sigmaSq = 0.2d;
            const double K = startK;
            const int M = 20;
            const double S0eps = 1e-5;
            const double T = 1d;
            const double tau = (T - 0d) / M;
            this.PrintParamsForSeries(r, sigmaSq, K, M, tau);
            PrintTableHeader(Nsteps, startN);

            var table = new Dictionary<string, List<double>>();

            for (double ki = 1; ki <= Ksteps; ki++)
            {
                var Ki = K * ki;
                for (var i = 0; i < Nsteps; i++)
                {
                    var n = (int) Math.Pow(2, i) * startN;
                    var folderPath = this.CreateOutputFolder(Ki, n, string.Empty);
                    var parameters = this.GetSeriesParameters(n, Ki, M, tau, a, b, r, sigmaSq, S0eps, folderPath);
                    var calculator = new AmericanOptionCalculator(parameters, false, false);
                    double[] S0Arr = calculator.Solve();

                    for (var tl = tls; tl >= 1; tl--)
                    {
                        var time = tau * (5d * tl);
                        var timeIndex = (int) (time / tau);
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
            const int startN = 150;
            const int Nsteps = 7; // 3;

            const double startK = 5d;
            const double  a = 0d;
            const double  b = 50d;
            const double r = 0.1d;
            const double sigmaSq = 0.2d;
            const double K = startK;
            const double S0eps = 1e-5;
            const double T = 1d;
            this.PrintParamsForSeries2(r, sigmaSq, K);
            PrintTableHeader(Nsteps, startN);

            var watch = new Stopwatch();
            watch.Reset();
            int M;
            double tau;
            double[] S0ArrGold;
            Console.WriteLine("Started Nsteps - 1");
            watch.Reset();
            watch.Start();
            {
                M = (int) (10 * Math.Pow(4, Nsteps - 1));
                tau = (T - 0d) / M;
                var n = (int) Math.Pow(2, Nsteps - 1) * startN;
                var folderPath = this.GetTrashFolder();

                /*
                    if (allowOutputFile)
                    {
                        folderPath = CreateOutputFolder(K, n, "convergence");
                    }
                */
                var parameters = this.GetSeriesParameters(n, K, M, tau, a, b, r, sigmaSq, S0eps, folderPath);
                var calculator = new AmericanOptionCalculator(parameters, false, false);
                S0ArrGold = calculator.Solve();
            }

            Console.WriteLine("Finished Nsteps - 1");
            var dictionary = new Dictionary<int, double>();
            Console.WriteLine("Started from 0 to Nsteps - 1");
            for (var i = 0; i < Nsteps - 1; i++)
            {
                Console.WriteLine("Started step " + i);
                M = (int) (10 * Math.Pow(4, i));
                tau = (T - 0d) / M;
                var n = (int) Math.Pow(2, i) * startN;
                var folderPath = this.GetTrashFolder();

                /*
                    if (allowOutputFile)
                    {
                        folderPath = CreateOutputFolder(K, n, "convergence");
                    }
                */
                var parameters = this.GetSeriesParameters(n, K, M, tau, a, b, r, sigmaSq, S0eps, folderPath);
                var calculator = new AmericanOptionCalculator(parameters, false, false);
                var S0Arr = calculator.Solve();
                dictionary[n] = this.GetErrorLInf(S0ArrGold, S0Arr);
                Console.WriteLine("Finished step = " + i);
            }

            watch.Stop();
            Console.WriteLine("Elapsed = " + watch.Elapsed.TotalSeconds + " s.");
            var list = new List<double>();
            foreach (KeyValuePair<int, double> pair in dictionary)
            {
                list.Add(pair.Value);
                Console.WriteLine(pair.Key + " " + pair.Value);
            }

            Console.WriteLine();
            Console.WriteLine("Log2");
            for (var i = 1; i < list.Count; i++)
            {
                var value = list[i - 1] / list[i];
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

        private double GetErrorLInf(IReadOnlyList<double> gold, IReadOnlyList<double> sol)
        {
            var stride = (int) ((double) gold.Count / sol.Count);
            var nSol = new double[gold.Count];
            var error = new double[gold.Count];
            for (int i = 0, k = 0; i < gold.Count; i += stride, k++)
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
            foreach (var pair in table)
            {
                string[] strings = pair.Key.Split(new[] {' '}, StringSplitOptions.RemoveEmptyEntries);
                var (item1, item2) = Tuple.Create(strings[0], strings[1]);
                var tupleItem1 = item1 + ";" + item2;
                tupleItem1 = pair.Value.Aggregate(tupleItem1, (current, tt1) => current + (";" + tt1.ToString("0.0000000")));

                dictionary[pair.Key] = tupleItem1;
            }

            foreach (var pair in dictionary)
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
                var n = (int) Math.Pow(2, i) * startN;
                Console.Write("N(" + n + ")" + new string(' ', 15));
            }

            Console.WriteLine();
            Console.WriteLine(new string('=', whitespaceNumber));
        }

        private static void PrintTable(Dictionary<string, List<double>> table)
        {
            foreach (KeyValuePair<string, List<double>> pair in table)
            {
                var value = pair.Value.Aggregate(pair.Key, (current, t) => current + (t.ToString("0.0000000") + new string(' ', t > 10d ? 10 : 11)));
                Console.WriteLine(value);
            }
        }

        private void Print(string filename, IReadOnlyList<double> St, double tau)
        {
            var name = $"{filename}_nx={St.Count}_tau={tau}.dat";
            using (var writer = new StreamWriter(name, false))
            {
                writer.WriteLine("TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA'\nVARIABLES = S {0}\nZONE T='{1}'",
                    "t", "SubZone");
                writer.WriteLine($"I={St.Count} K={1} ZONETYPE=Ordered");
                writer.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE)");
                for (var i = 0; i < St.Count; i++)
                {
                    writer.WriteLine("{0:e8} {1:e8}", St[i], i);
                }
            }
        }

        private AmericanOptionParameters GetParameters(bool saveSolutions, string workPath)
        {
            const double a = 0d;
            const double b = 8d;
            const double r = 0.1d;
            const int n = 300;
            const double tau = 1e-5d;
            const double sigmaSq = 0.2d;
            const double eps = 1e-5d;
            const double K = 5d;
            const int M = 365;
            return new AmericanOptionParameters(a, b, n, r, tau, sigmaSq, K, eps, M, workPath, saveSolutions);
        }

        private AmericanOptionParameters GetSeriesParameters(
            int n,
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
            return new AmericanOptionParameters(a, b, n, r, tau, sigma, K, S0eps, M, workDir, false);
        }

        private static void PrintParameters(AmericanOptionCalculator calculator)
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