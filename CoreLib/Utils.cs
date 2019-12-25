namespace PerpetualAmericanOptions
{
    using System;

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

        public static double[] GetError(double[] exact, double[] num)
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

        public static double GetL1(double h, double[] data)
        {
            var r = 0.0;
            for (var i = 0; i < data.Length; ++i)
            {
                r += Math.Abs(data[i]);
            }

            return r * h;
        }

        public static double GetLInf(double[] data)
        {
            var mx = 0.0;
            for (var i = 0; i < data.Length; ++i)
            {
                mx = Math.Max(data[i], mx);
            }

            return mx;
        }
    }
}