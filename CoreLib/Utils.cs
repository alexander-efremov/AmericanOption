namespace CoreLib.Utils
{
    using System;
    using System.Collections.Generic;
    using System.Linq;

    public static class Utils
    {
        public static void Print(double[] arr, string header)
        {
            Console.WriteLine(header);
            for (var i = 0; i < arr.Length; i++)
            {
                Console.Write(arr[i] + " ");
            }

            Console.WriteLine();
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

        public static double[] GetAbsError(double[] exact, double[] num)
        {
            var err = new double[exact.Length];

            for (var i = 0; i < exact.Length; ++i)
            {
                err[i] = Math.Abs(exact[i] - num[i]);
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
    }
}