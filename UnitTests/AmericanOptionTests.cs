namespace UnitTests
{
    using System;
    using System.Collections.Generic;
    using System.Diagnostics;
    using System.Globalization;
    using System.IO;
    using System.Linq;
    using System.Text;
    using AmericanOption;
    using AmericanOptions;
    using CoreLib;
    using CoreLib.Utils;
    using NUnit.Framework;

    [TestFixture]
    public class AmericanOptionTests : UnitTestBase
    {
        [Test]
        public void AmericanOption()
        {
            var parameters = this.GetParameters(true, this.GetWorkingDir() + "AO/");
            var calculator = new AmericanOptionCalculator(parameters, true, false);

            PrintParameters(calculator);

            double[] S0t = calculator.Solve();
            this.ConvexityCheck(S0t);
            List<SolutionData> numericSolutions = calculator.GetNumericSolutions();

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
                    sb.Append(d.VS.ToString("e8", CultureInfo.InvariantCulture) + " ");
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
                    sb.Append(d.VS.ToString("e8", CultureInfo.InvariantCulture) + " ");
                }

                Console.WriteLine(numericSolution.k + " " + numericSolution.S0.ToString("e8", CultureInfo.InvariantCulture) + " " + sb);
                if (index > 10)
                {
                    break;
                }
            }

            // Console.WriteLine();
            // Console.WriteLine("Print Errors");
            // for (var index = 0; index < numericSolutions.Count; index++)
            // {
            //     var exactSolution = exactSolutions[index];
            //     var numericSolution = numericSolutions[index];
            //     Point[] error = Utils.GetAbsError(exactSolution.Solution, numericSolution.Solution);
            //     var sb = new StringBuilder();
            //     foreach (var d in error.Take(5))
            //     {
            //         sb.Append(d.VS.ToString("e8", CultureInfo.InvariantCulture) + " ");
            //     }
            //
            //     Console.WriteLine(numericSolution.k + " " + numericSolution.S0.ToString("e8", CultureInfo.InvariantCulture) + " " + sb);
            //     if (index > 10)
            //     {
            //         break;
            //     }
            // }

            ClearData(parameters);

            var printer = calculator.GetTecplotPrinter();

            // checks solutions from the T to 0 time
            for (var i = 0; i < numericSolutions.Count; i++)
            {
                var exactSolution = exactSolutions[i];
                var numericSolution = numericSolutions[i];
                // Point[] error = Utils.GetAbsError(exactSolution.Solution, numericSolution.Solution);

                printer.PrintXY(
                    Path.Combine(parameters.WorkDir + "/exact/", "exactSolution"),
                    calculator.GetTau() * (calculator.GetM() - i),
                    calculator.GetH(),
                    exactSolution.Solution,
                    calculator.GetM() - i,
                    "exact_" + (calculator.GetM() - i));
                printer.PrintXY(
                    Path.Combine(parameters.WorkDir + "/numeric/", "numericSolution"),
                    calculator.GetTau() * (calculator.GetM() - i),
                    calculator.GetH(),
                    numericSolution.Solution,
                    calculator.GetM() - i,
                    "numeric_" + (calculator.GetM() - i));
                // printer.PrintXY(
                //     Path.Combine(parameters.WorkDir + "/error/", "error"),
                //     calculator.GetTau() * (calculator.GetM() - i),
                //     calculator.GetH(),
                //     error,
                //     calculator.GetM() - i,
                //     "error_" + (calculator.GetM() - i));
                // if (i > 2)
                // {
                //     break;
                // }
            }

            printer.PrintXYZ(
                Path.Combine(parameters.WorkDir + "/exact/", "exactSolution3d"),
                calculator.GetH(),
                calculator.GetTau(),
                exactSolutions,
                0,
                "exact");
            printer.PrintXYZ(
                Path.Combine(parameters.WorkDir + "/numeric/", "numericSolution3d"),
                calculator.GetH(),
                calculator.GetTau(),
                numericSolutions,
                0,
                "numeric");

            // // checks solutions from the T to 0 time
            // for (var i = 0; i < numericSolutions.Count; i++)
            // {
            //     var exactSolution = exactSolutions[i];
            //     var numericSolution = numericSolutions[i];
            //     Point[] error = Utils.GetAbsError(exactSolution.Solution, numericSolution.Solution);
            //
            //     Assert.AreEqual(exactSolution.Solution.Length, numericSolution.Solution.Length);
            //
            //     for (var j = 0; j < exactSolution.Solution.Length; j++)
            //     {
            //         try
            //         {
            //             Assert.AreEqual(
            //                 exactSolution.k,
            //                 numericSolution.k,
            //                 "k is wrong");
            //             Assert.AreEqual(
            //                 exactSolution.S0,
            //                 numericSolution.S0,
            //                 "S0 is wrong");
            //             Assert.AreEqual(
            //                 exactSolution.Solution[j],
            //                 numericSolution.Solution[j],
            //                 // 10e-6,
            //                 "ex: " + exactSolution.Solution[j].VS.ToString("e8", CultureInfo.InvariantCulture)
            //                        + " num: " + numericSolution.Solution[j].VS.ToString("e8", CultureInfo.CurrentCulture)
            //                        + " err: " + error[j].VS.ToString("e8", CultureInfo.CurrentCulture)
            //                        + " rel: " + numericSolution.Solution[j].VS / exactSolution.Solution[j].VS);
            //         }
            //         catch (Exception)
            //         {
            //             Console.WriteLine("i = " + i + " j = " + j);
            //             throw;
            //         }
            //     }
            // }

            this.Print(Path.Combine(parameters.WorkDir, "s0"), S0t, calculator.GetTau());
        }

        [Test]
        public void AmericanOptionNumericSolutionDraw()
        {
            var parameters = this.GetParameters(true, this.GetWorkingDir() + "AO/");

            var calculator = new AmericanOptionCalculator(parameters, true, false);
            calculator.Solve();
            List<SolutionData> numericSolutions = calculator.GetNumericSolutions();
            ClearData(parameters);

            var printer = calculator.GetTecplotPrinter();
            for (var index = 0; index < numericSolutions.Count; index++)
            {
                var solution = numericSolutions[index];
                // calculator.UpdateH(solution.S0);
                var n0toS0 = (int)(Math.Floor(solution.S0 - calculator.Geta()) / calculator.GetH());
                var sol = new List<Point>();

                for (var i = 0; i <= n0toS0; i++)
                {
                    sol[i] = new Point(calculator.GetK() - (calculator.Geta() + i * calculator.GetH()), calculator.Geta() + i * calculator.GetH());
                }

                sol.AddRange(solution.Solution);

                printer.PrintXY(
                    Path.Combine(parameters.WorkDir + "/numeric/", "numericSolution"),
                    calculator.GetTau() * (calculator.GetM() - index),
                    calculator.GetH(),
                    sol.ToArray(),
                    calculator.GetM() - index,
                    "numeric_" + (calculator.GetM() - index),
                    "here the first point of numeric V(S,t)",
                    n0toS0);
            }
        }

        [Test]
        public void AmericanOptionSolutionsDrawing()
        {
            var parameters = this.GetParameters(true, this.GetWorkingDir() + "AO/");
            var calculator = new AmericanOptionCalculator(parameters, true, false);
            double[] S0t = calculator.Solve();
            List<SolutionData> numericSolutions = calculator.GetNumericSolutions();
            ClearData(parameters);

            PrintParameters(calculator);
            ClearData(parameters);

            List<SolutionData> exactSolutions = calculator.GetExactSolutions2(S0t);
            var printer = calculator.GetTecplotPrinter();
            for (var index = 0; index < exactSolutions.Count; index++)
                //for (var index = 0; index < 2; index++)
            {
                var exactSolution = exactSolutions[index];
                // calculator.UpdateH(exactSolution.S0);

                printer.PrintXY(
                    Path.Combine(parameters.WorkDir + "/exact/", "exactSolution"),
                    calculator.GetTau() * (calculator.GetM() - index),
                    calculator.GetH(),
                    exactSolution.Solution,
                    calculator.GetM() - index,
                    "exact_" + (calculator.GetM() - index),
                    "here the first point of exact V(S,t)",
                    exactSolution.StartPosOfS0);
            }

            for (var index = 0; index < numericSolutions.Count; index++)
                //for (var index = 0; index < 2; index++)
            {
                var solution = numericSolutions[index];
                // calculator.UpdateH(solution.S0);

                var sol = new List<Point>();
                var k = 0;
                while (k * calculator.GetH() < solution.S0)
                {
                    var point = new Point(k * calculator.GetH(), calculator.GetK() - k * calculator.GetH());
                    sol.Add(point);
                    k++;
                }

                sol.AddRange(solution.Solution);

                printer.PrintXY(
                    Path.Combine(parameters.WorkDir + "/numeric/", "numericSolution"),
                    calculator.GetTau() * (calculator.GetM() - index),
                    calculator.GetH(),
                    sol.ToArray(),
                    calculator.GetM() - index,
                    "numeric_" + (calculator.GetM() - index),
                    "here the first point of numeric V(S,t)",
                    k);
            }
        }

        [Test]
        public void TestSeries()
        {
            const int startN = 4800;
            const int Nsteps = 5;
            const double startK = 5d;

            const double a = 0d;
            const double b = 10d;
            const double alpha = 1d;
            const double beta = 1d;
            const double r = 0.1d;

            const double sigmaSq = 0.2d;
            const double S0Eps = 1e-5d;
            const double K = startK;
            const double h = (b - K) / startN;
            const int startM = 365;
            const double tauStart = 1e-6d;
            const double smoothness = 1000d;
            var T = startM * tauStart;
            
            this.PrintParamsForSeries(r, sigmaSq, K, startM, tauStart);
            PrintTableHeader(Nsteps, startN);

            var table = new Dictionary<string, List<double>>();

            List<double[]> list = new List<double[]>();
            for (var i = 0; i < Nsteps; i++)
            {
                var n = (int)Math.Pow(2, i) * startN;
                double tau = tauStart  ;
                var d = (2 * (i + 1));
                if (d > 0)
                {
                    tau = tau / d;
                }
               
                var M = (int)(T/tau);

                var folderPath = this.CreateOutputFolder(M, n, string.Empty);
                var parameters = this.GetSeriesParameters(alpha, beta, n, K, M, tau, a, b, r, sigmaSq, S0Eps, h, smoothness, folderPath);
                PrintParameters(parameters);
                var calculator = new AmericanOptionCalculator(parameters, false, false);
                double[] S0Arr = calculator.Solve();
                
                 list.Add(S0Arr);
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
            const double a = 0d;
            const double b = 50d;
            const double alpha = 1d;
            const double beta = 1d;
            const double r = 0.1d;
            const double sigmaSq = 0.2d;
            const double K = startK;
            const double S0eps = 1e-5;
            const double h = b / (2 * startN);
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
                var parameters = this.GetSeriesParameters(alpha, beta, n, K, M, tau, a, b, r, sigmaSq, S0eps, h, 1000d, folderPath);
                var calculator = new AmericanOptionCalculator(parameters, false, false);
                S0ArrGold = calculator.Solve();
            }

            Console.WriteLine("Finished Nsteps - 1");
            var dictionary = new Dictionary<int, double>();
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
                var parameters = this.GetSeriesParameters(alpha, beta, n, K, M, tau, a, b, r, sigmaSq, S0eps, h, 1000d, folderPath);
                var calculator = new AmericanOptionCalculator(parameters, false, false);
                double[] S0Arr = calculator.Solve();
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
        
        [Test]
        public void TestSeriesConvergence2()
        {
            const int startN = 4800;
            const int Nsteps = 5;
            const double startK = 5d;

            const double a = 0d;
            const double b = 10d;
            const double alpha = 1d;
            const double beta = 1d;
            const double r = 0.1d;

            const double sigmaSq = 0.2d;
            const double S0Eps = 1e-5d;
            const double K = startK;
            const double h = (b - K) / startN;
            const int startM = 365;
            const double tauStart = 1e-6d;
            const double smoothness = 1000d;
            this.PrintParamsForSeries2(r, sigmaSq, K);
            PrintTableHeader(Nsteps, startN);
            var T = startM * tauStart;
            double[] S0ArrGold;
            Console.WriteLine("Started Nsteps - 1");
            {
                 var n = (int)Math.Pow(2, (Nsteps-1)) * startN;
                 double tau = tauStart  ;
                 var d = (2 * ((Nsteps-1) + 1));
                 if (d > 0)
                 {
                     tau = tau / d;
                 }
               
                 var M = (int)(T/tau);
                 var folderPath = this.GetTrashFolder();

                /*
                    if (allowOutputFile)
                    {
                        folderPath = CreateOutputFolder(K, n, "convergence");
                    }
                */
                var parameters = this.GetSeriesParameters(alpha, beta, n, K, M, tau, a, b, r, sigmaSq, S0Eps, h, smoothness, folderPath);
                PrintParameters(parameters);
                var calculator = new AmericanOptionCalculator(parameters, false, false);
                S0ArrGold = calculator.Solve();
            }

            Console.WriteLine("Finished Nsteps - 1");
            
            var dictionary = new Dictionary<int, double>();
            Console.WriteLine("Started from 0 to Nsteps - 1");
            for (var i = 0; i < Nsteps - 1; i++)
            {
                var n = (int)Math.Pow(2, i) * startN;
                double tau = tauStart  ;
                var d = (2 * (i + 1));
                if (d > 0)
                {
                    tau = tau / d;
                }
               
                var M = (int)(T/tau);
                 
                var folderPath = this.GetTrashFolder();

                /*
                    if (allowOutputFile)
                    {
                        folderPath = CreateOutputFolder(K, n, "convergence");
                    }
                */
                var parameters = this.GetSeriesParameters(alpha, beta, n, K, M, tau, a, b, r, sigmaSq, S0Eps, h, smoothness, folderPath);
                PrintParameters(parameters);
                var calculator = new AmericanOptionCalculator(parameters, false, false);
                double[] S0Arr = calculator.Solve();
                dictionary[n] = this.GetErrorLInf(S0ArrGold, S0Arr);
                Console.WriteLine("Finished step = " + i);
            }

             
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

        private static void ClearData(Parameters parameters)
        {
            var di = new DirectoryInfo(parameters.WorkDir);
            FileInfo[] files = di.GetFiles("*.dat")
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

            di = new DirectoryInfo(parameters.WorkDir + "numeric/");
            files = di.GetFiles("*.dat")
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

            di = new DirectoryInfo(parameters.WorkDir + "exact/");
            files = di.GetFiles("*.dat")
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

            di = new DirectoryInfo(parameters.WorkDir + "error/");
            files = di.GetFiles("*.dat")
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
            Console.WriteLine("h = " + calculator.GetH());
            Console.WriteLine("S0 Eps = " + calculator.GetS0Eps());
            Console.WriteLine();
        }

        private static void PrintParameters(AmericanOptionParameters calculator)
        {
            Console.WriteLine("a = " + calculator.A);
            Console.WriteLine("b = " + calculator.B);
            Console.WriteLine("r = " + calculator.R);
            Console.WriteLine("N = " + calculator.N);
            Console.WriteLine("N_1 = " + calculator.N1);
            Console.WriteLine("tau = " + calculator.Tau);
            Console.WriteLine("sigma_sq = " + calculator.SigmaSq);
            Console.WriteLine("K = " + calculator.K);
            Console.WriteLine("M = " + calculator.M);
            Console.WriteLine("T = " + calculator.T);
            //Console.WriteLine("h = " + calculator.h);
            Console.WriteLine("S0 Eps = " + calculator.S0Eps);
            Console.WriteLine();
        }

        private static void PrintTable(Dictionary<string, List<double>> table)
        {
            foreach (var value in
                table.Select(pair => pair.Value.Aggregate(pair.Key, (current, t) => current + t.ToString("0.0000000") + new string(' ', t > 10d ? 10 : 11))))
            {
                Console.WriteLine(value);
            }
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

        private void ConvexityCheck(IReadOnlyList<double> S0t)
        {
            for (var i = S0t.Count - 1; i >= 1; i--)
            {
                Assert.IsTrue(S0t[i] > S0t[i - 1]);
            }

            Console.WriteLine("Numeric S0 (first 10):");
            for (var i = S0t.Count - 1; i >= S0t.Count - 10; i--)
            {
                Console.WriteLine("k = " + (i + 1) + " -> " + S0t[i]);
            }
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

        private double GetErrorLInf(IReadOnlyList<double> gold, IReadOnlyList<double> sol)
        {
            var stride = (int)((double)gold.Count / sol.Count);
            var nSol = new double[gold.Count];
            var error = new double[gold.Count];
            for (int i = 0, k = 0; i < gold.Count; i += stride, k++)
            {
                try
                {
                    nSol[i] = sol[k];
                    error[i] = Math.Abs(gold[i] - nSol[i]);
                }
                catch (Exception e)
                {
                    Console.WriteLine(e);
                    throw;
                }
            }

            // var error = Utils.GetAbsError(gold, nSol);
            return Utils.GetLInf(error);
        }

        private AmericanOptionParameters GetParameters(bool saveSolutions, string workPath)
        {
            const double a = 0d;
            const double b = 10d;
            const double alpha = 1d;
            const double beta = 1d;
            const double r = 0.1d;
            const int n = 1200;
            const double tau = 1e-5d;
            const double sigmaSq = 0.2d;
            const double S0Eps = 1e-5d;
            const double K = 5d;
            const double h = (b - K) / n;
            const int M = 365;
            const double smoothness = 1000d;
            return new AmericanOptionParameters(alpha, beta, a, b, n, r, tau, sigmaSq, K, S0Eps, h, M, workPath, saveSolutions, smoothness);
        }

        private AmericanOptionParameters GetSeriesParameters(
            double alpha,
            double beta,
            int n,
            double K,
            int M,
            double tau,
            double a,
            double b,
            double r,
            double sigma,
            double S0eps,
            double h,
            double smoothness,
            string workDir)
        {
            return new AmericanOptionParameters(alpha, beta, a, b, n, r, tau, sigma, K, S0eps, h, M, workDir, false, smoothness);
        }

        private string GetTrashFolder()
        {
            return this.GetWorkingDir() + "AO/trash";
        }

        private string GetWorkingDir()
        {
            return Environment.GetFolderPath(Environment.SpecialFolder.Desktop) + Path.DirectorySeparatorChar;
        }

        private void Print(string filename, IReadOnlyList<double> St, double tau)
        {
            var name = $"{filename}_nx={St.Count}_tau={tau}.dat";
            using (var writer = new StreamWriter(name, false))
            {
                writer.WriteLine(
                    "TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA'\nVARIABLES = S {0}\nZONE T='{1}'",
                    "t",
                    "SubZone");
                writer.WriteLine($"I={St.Count} K={1} ZONETYPE=Ordered");
                writer.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE)");
                for (var i = 0; i < St.Count; i++)
                {
                    writer.WriteLine("{0:e8} {1:e8}", St[i], i);
                }
            }
        }

        private void PrintCsv(Dictionary<string, List<double>> table)
        {
            var dictionary = new Dictionary<string, string>();
            foreach (KeyValuePair<string, List<double>> pair in table)
            {
                string[] strings = pair.Key.Split(new[] {' '}, StringSplitOptions.RemoveEmptyEntries);
                var (item1, item2) = Tuple.Create(strings[0], strings[1]);
                var tupleItem1 = item1 + ";" + item2;
                tupleItem1 = pair.Value.Aggregate(tupleItem1, (current, tt1) => current + ";" + tt1.ToString("0.0000000"));

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
    }
}