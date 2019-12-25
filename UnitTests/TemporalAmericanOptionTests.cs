using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using NUnit.Framework;
using TemporalAmericanOption;

namespace PerpetualAmericanOptions
{
    public class TemporalAmericanOptionTests : UnitTestBase
    {
        [Test]
        public void TemporalAmericanOption()
        {
            var parameters = GetParameters();
            var calculator = new TemporalAmericanOptionCalculator(parameters, true, true, Environment.GetFolderPath(Environment.SpecialFolder.Desktop) + "\\AO\\");

            PrintParameters(calculator);
            Console.WriteLine();

            double[] S0Arr = calculator.Solve();

            Console.WriteLine("Calculated S0");
            for (var i = S0Arr.Length - 1; i >= 0; i--) Console.WriteLine("k = " + (i + 1) + " -> " + S0Arr[i]);

            Print(Environment.GetFolderPath(Environment.SpecialFolder.Desktop) + "\\AO\\" + "s0", S0Arr,
                calculator.GetTau());

            Console.WriteLine();
        }

        [Test]
        public void TestSeries()
        {
            var allowOutputFile = false;
            var allowOutputConsole = false;
            var startN = 150;
            int Ksteps = 3;
            int Nsteps = 6;
            int tls = 3;
            double startK = 5d;
            var a = 0d;
            var b = 50d;
            var r = 0.1d;
            var sigmaSq = 0.2d;
            var K = startK;
            var M = 20;
            var S0eps = 1e-5;
            double T = 1d;
            var tau = (T - 0d) / M;
            PrintParamsForSeries(r, sigmaSq, K, M, tau);
            PrintTableHeader(Nsteps, startN);

            var table = new Dictionary<string, List<double>>();
            
            for (double ki = 1; ki <= Ksteps; ki++)
            {
                var Ki = K * ki;
                for (var i = 0; i < Nsteps; i++)
                {
                    var n = (int) Math.Pow(2, i) * startN;
                    var folderPath = CreateOutputFolder(Ki, n, string.Empty);
                    var parameters = GetSeriesParameters(n, T, Ki, M, tau, a, b, r, sigmaSq, S0eps);
                    var calculator = new TemporalAmericanOptionCalculator(parameters, allowOutputFile, allowOutputConsole, folderPath);
                    double[] S0Arr = calculator.Solve();
                    
                    for (int tl = tls; tl >= 1; tl--)
                    {
                        double time = tau * (5d * tl);
                        var timeIndex = (int) (time / tau);
                        var value = S0Arr[timeIndex];
                        var s1 = value;
                        var format = Ki + new string(' ', Ki > 5 ? 9 : 10) +
                                     (time.ToString("0.00") + new string(' ', 10));
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
            PrintCsv(table);
        }
        
        [Test]
        public void TestSeriesConvergence()
        {
            var allowOutputFile = false;
            var allowOutputConsole = false;
            var startN = 150;
//            int Ksteps = 3;
            int Nsteps = 7;
//            int tls = 3;
            double startK = 5d;
            var a = 0d;
            var b = 50d;
            var r = 0.1d;
            var sigmaSq = 0.2d;
            var K = startK;
            var S0eps = 1e-5;
            double T = 1d;
            PrintParamsForSeries2(r, sigmaSq, K);
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
                M = (int) (10 * Math.Pow(4, (Nsteps-1)));
                tau = (T - 0d) / M;
                var n = (int) Math.Pow(2, Nsteps - 1) * startN;
                string folderPath = GetTrashFolder();
/*
                if (allowOutputFile)
                {
                    folderPath = CreateOutputFolder(K, n, "convergence");
                }
*/

                var parameters = GetSeriesParameters(n, T, K, M, tau, a, b, r, sigmaSq, S0eps);
                var calculator =
                    new TemporalAmericanOptionCalculator(parameters, allowOutputFile, allowOutputConsole, folderPath);
                S0ArrGold = calculator.Solve();
            }
            Console.WriteLine("Finished Nsteps - 1");
            var dic = new Dictionary<int, double>();
            Console.WriteLine("Started from 0 to Nsteps - 1");
            for (var i = 0; i < Nsteps - 1; i++)
            {
                Console.WriteLine("Started step " + i);
                M = (int) (10 * Math.Pow(4, i));
                tau = (T - 0d) / M;
                var n = (int) Math.Pow(2, i) * startN;
                string folderPath = GetTrashFolder();
                if (allowOutputFile)
                {
                    folderPath = CreateOutputFolder(K, n, "convergence");
                }
                var parameters = GetSeriesParameters(n, T, K, M, tau, a, b, r, sigmaSq, S0eps);
                var calculator =
                    new TemporalAmericanOptionCalculator(parameters, allowOutputFile, allowOutputConsole, folderPath);
                double[] S0Arr = calculator.Solve();
                double lInf = GetErrorLInf(S0ArrGold, S0Arr);
                dic[n] = lInf;
                Console.WriteLine("Finished step = "  + i);
            }
            w.Stop();
            Console.WriteLine("Elapsed = " + w.Elapsed.TotalSeconds + " s.");
            var l = new List<double>();
            foreach (KeyValuePair<int,double> pair in dic)
            {
                l.Add(pair.Value);
                Console.WriteLine(pair.Key + " " + pair.Value);
            }

            Console.WriteLine();
            Console.WriteLine("Log2");
            for (int i = 1; i < l.Count; i++)
            {
                var value = l[i-1]/l[i];
                Console.WriteLine(value + " " + Math.Log(value, 2));
            }

//            PrintTable(table);
//            Console.WriteLine();
//            PrintCsv(table);
        }

        private double GetErrorLInf(double[] gold, double[] sol)
        {
            var stride = (int)((double)gold.Length / sol.Length);
            var nSol = new double[gold.Length];
            var error = new double[gold.Length];
            for (int i = 0, k = 0; i < gold.Length; i+=stride, k++)
            {
                nSol[i] = sol[k];
                error[i] = Math.Abs(gold[i] - nSol[i]);
            }

//            double[] error = Utils.GetAbsError(gold, nSol);
            return Utils.GetLInf(error);
        }

        private string GetTrashFolder()
        {
            return Environment.GetFolderPath(Environment.SpecialFolder.Desktop) + "\\" + "AO/trash";
        }

        private void PrintCsv(Dictionary<string, List<double>> table)
        {
            var dictionary = new Dictionary<string, string>();
            foreach (KeyValuePair<string,List<double>> pair in table)
            {
                string[] strings = pair.Key.Split(new []{' '},StringSplitOptions.RemoveEmptyEntries);
                Tuple<string, string> tuple = Tuple.Create(strings[0], strings[1]);
                var tupleItem1 = tuple.Item1+";"+tuple.Item2;
                foreach (var tt1 in pair.Value)
                {
                    tupleItem1 +=";" + tt1.ToString("0.0000000");
                }
                dictionary[pair.Key] = tupleItem1;
            }

            foreach (KeyValuePair<string,string> pair in dictionary)
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
        
        private void PrintParamsForSeries2(double r, double sigmaSq, double K )
        {
            Console.Write("r = " + r); 
            Console.Write(" sigma_sq = " + sigmaSq);
            Console.WriteLine(" K = " + K);
        }

        private static string CreateOutputFolder(double Ki, int n, string subfolder)
        {
            var folderPath = Environment.GetFolderPath(Environment.SpecialFolder.Desktop) + "/AO/" + subfolder + "/" + Ki + "_" + "_" + n + "/";
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
            var name = string.Format("{0}_nx={1}_tau={2}.dat", filename, St.Length, tau);
            using (var writer = new StreamWriter(name, false))
            {
                writer.WriteLine(
                    "TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA'\nVARIABLES = S {0}\nZONE T='{1}'", "t",
                    "SubZone");
                writer.WriteLine("I={0} K={1} ZONETYPE=Ordered", St.Length, 1);
                writer.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE)");
                for (var i = 0; i < St.Length; i++) writer.WriteLine("{0:e8} {1:e8}", St[i], i);
            }
        }

        private TemporalParameters GetParameters()
        {
            var a = 0d;
            var b = 50d;
            var r = 0.1d;
            var sigmaSq = 0.2d;
            var K = 5d;
            var M = 20;
            var S0eps = 1e-5;
            var T = 1d;
            var tau = (T - 0d) / M;
            var n = 100;
            return new TemporalParameters(a, b, n, r, tau, sigmaSq, K, S0eps, M, T);
        }

        private TemporalParameters GetSeriesParameters(int n, double T, double K, int M, double tau, double a, double b,
            double r, double sigma, double S0eps)
        {
            return new TemporalParameters(a, b, n, r, tau, sigma, K, S0eps, M, T);
        }

        private static void PrintParameters(TemporalAmericanOptionCalculator calculator)
        {
            Console.WriteLine("b = " + calculator.GetRightBoundary());
            Console.WriteLine("r = " + calculator.GetR());
            Console.WriteLine("N = " + calculator.GetN());
            Console.WriteLine("N_1 = " + calculator.GetN1());
            Console.WriteLine("tau = " + calculator.GetTau());
            Console.WriteLine("sigma_sq = " + calculator.GetSquaredSigma());
            Console.WriteLine("K = " + calculator.GetK());
            Console.WriteLine("M = " + calculator.GetM());
            Console.WriteLine("S0 Eps = " + calculator.GetS0Eps());
        }

        protected override string SetWorkingDir()
        {
            return Environment.GetFolderPath(Environment.SpecialFolder.Desktop) + "\\";
        }
    }
}