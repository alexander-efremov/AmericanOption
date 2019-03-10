using System;

namespace PerpetualAmericanOptions
{
    internal static class Utils
    {
        internal static void Print(double[] arr, string header)
        {
            Console.WriteLine(header);
            for (int i = 0; i < arr.Length; i++)
            {
                Console.Write(arr[i] + " ");
            }

            Console.WriteLine();
        }

        internal static void print_result_arrays1(int exp_cnt, double[] l1_10, double[] l1_40, double[] l1_160,
            double[] l1_640, double[] l1_2560)
        {
            if (l1_10 != null)
            {
                for (var i = 0; i < exp_cnt; ++i) Console.WriteLine("{0};", l1_10[i]);
                Console.WriteLine("\n");
            }

            if (l1_40 != null)
            {
                for (var i = 0; i < exp_cnt; ++i) Console.WriteLine("{0};", l1_40[i]);
                Console.WriteLine("\n");
            }

            if (l1_160 != null)
            {
                for (var i = 0; i < exp_cnt; ++i) Console.WriteLine("{0};", l1_160[i]);
                Console.WriteLine("\n");
            }

            if (l1_640 != null)
            {
                for (var i = 0; i < exp_cnt; ++i) Console.WriteLine("{0};", l1_640[i]);
                Console.WriteLine("\n");
            }

            if (l1_2560 != null)
            {
                for (var i = 0; i < exp_cnt; ++i) Console.WriteLine("{0};", l1_2560[i]);
                Console.WriteLine("\n");
            }
        }

        internal static double[] GetError(double[] arr1, double[] arr2, int n)
        {
            var err = new double[n];

            for (var i = 0; i < n; ++i) err[i] = arr1[i] - arr2[i];

            return err;
        }

        internal static double GetL1(double h, int n, double[] data)
        {
            var r = 0.0;
            for (var i = 0; i < n; ++i)
            {
                r += Math.Abs(data[i]);
            }
            return r * h;
        }

        internal static double[] FillArrayDiff(double[] arr1, double[] arr2)
        {
            double[] err = new double[arr1.Length];
            for (var i = 0; i < err.Length; ++i)
            {
                err[i] = arr1[i] - arr2[i];
            }

            return err;
        }

        internal static void
            print_header(int norm_type, bool check_solution, bool use_betas)
        {
            if (check_solution)
            {
                if (norm_type == 1)
                    Console.WriteLine("\nL1 of u");
                else if (norm_type == 2)
                    Console.WriteLine("\nL_inf of u");
                else if (norm_type == 3) Console.WriteLine("\nL_2 of u");
                if (use_betas)
                    Console.WriteLine(" using beta approach.\n");
                else
                    Console.WriteLine(" using integration approach,\n");
            }
            else
            {
                if (norm_type == 1)
                    Console.WriteLine("\nL1 of error");
                else if (norm_type == 2)
                    Console.WriteLine("\nL_inf of error");
                else if (norm_type == 3) Console.WriteLine("\nL_2 of error");
                if (use_betas)
                    Console.WriteLine(" using beta approach.\n");
                else
                    Console.WriteLine(" using integration approach.\n");
            }
        }

        internal static void
            print_result_arrays(int exp_cnt, double[] l1_10, double[] l1_40, double[] l1_160, double[] l1_640,
                double[] l1_2560)
        {
            if (l1_10 != null)
            {
                for (var i = 0; i < exp_cnt; ++i) Console.WriteLine("{0};", l1_10[i]);
                Console.WriteLine("\n");
            }

            if (l1_40 != null)
            {
                for (var i = 0; i < exp_cnt; ++i) Console.WriteLine("{0};", l1_40[i]);
                Console.WriteLine("\n");
            }

            if (l1_160 != null)
            {
                for (var i = 0; i < exp_cnt; ++i) Console.WriteLine("{0};", l1_160[i]);
                Console.WriteLine("\n");
            }

            if (l1_640 != null)
            {
                for (var i = 0; i < exp_cnt; ++i) Console.WriteLine("{0};", l1_640[i]);
                Console.WriteLine("\n");
            }

            if (l1_2560 != null)
            {
                for (var i = 0; i < exp_cnt; ++i) Console.WriteLine("{0};", l1_2560[i]);
                Console.WriteLine("\n");
            }
        }
    }
}