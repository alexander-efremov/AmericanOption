using System;

namespace PerpetualAmericanOptions
{
    public static class Utils
    {
        public static void Print(double[] arr, string header)
        {
            Console.WriteLine(header);
            for (int i = 0; i < arr.Length; i++)
            {
                Console.Write(arr[i] + " ");
            }

            Console.WriteLine();
        }

        public static void print_result_arrays1(int exp_cnt, double[] l1_10, double[] l1_40, double[] l1_160,
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

        public static double[] GetError(double[] exact, double[] num)
        {
            var err = new double[exact.Length];

            for (var i = 0; i < exact.Length; ++i) err[i] = exact[i] - num[i];

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

        public static double[] FillArrayDiff(double[] arr1, double[] arr2)
        {
            double[] err = new double[arr1.Length];
            for (var i = 0; i < err.Length; ++i)
            {
                try
                {
err[i] = arr1[i] - arr2[i];
                }
                catch (Exception e)
                {
                    Console.WriteLine(e);
                    throw;
                }
            }

            return err;
        }

        public static void print_header(int norm_type, bool check_solution, bool use_betas)
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

        public static void print_result_arrays(int exp_cnt, double[] l1_10, double[] l1_40, double[] l1_160, double[] l1_640,
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