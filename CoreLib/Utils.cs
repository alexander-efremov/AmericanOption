namespace CoreLib.Utils
{
    using System;
    using System.Collections.Generic;
    using System.Diagnostics.CodeAnalysis;
    using System.Linq;

    public static class Utils
    {
        [SuppressMessage("ReSharper", "UnusedMember.Global")]
        public static double[] GetAbsError(double[] exact, double[] num)
        {
            var err = new double[exact.Length];

            for (var i = 0; i < exact.Length; ++i)
            {
                err[i] = Math.Abs(exact[i] - num[i]);
            }

            return err;
        }

        public static Point[] GetAbsError(Point[] exact, Point[] num)
        {
            var err = new Point[exact.Length];

            for (var i = 0; i < exact.Length; ++i)
            {
                try
                {
                    var p = new Point(Math.Abs(exact[i].S - num[i].S), Math.Abs(exact[i].VS - num[i].VS));
                    err[i] = p;
                }
                catch (Exception e)
                {
                    Console.WriteLine(e);
                    throw;
                }
            }

            return err;
        }

        public static IEnumerable<double> GetError(double[] exact, double[] num)
        {
            var err = new double[exact.Length];

            for (var i = 0; i < exact.Length; ++i)
            {
                err[i] = exact[i] - num[i];
            }

            return err;
        }

        public static double GetL1(double h, IEnumerable<double> data)
        {
            var r = data.Sum(Math.Abs);
            return r * h;
        }

        public static double GetLInf(IEnumerable<double> data)
        {
            return data.Concat(new[] {0.0}).Max();
        }

        public static void Print(double[] arr, string header)
        {
            Console.WriteLine(header);
            for (var i = 0; i < arr.Length; i++)
            {
                Console.Write(arr[i] + " ");
            }

            Console.WriteLine();
        }
    }
}