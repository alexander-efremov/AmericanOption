using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;

namespace CoreLib
{
    public static class Utils
    {
        public static double GetL1(double h, IEnumerable<double> data) => data.Sum(Math.Abs) * h;

        public static double GetLInf(IEnumerable<double> data) => data.Concat(new[] {0.0}).Max();

        [SuppressMessage("ReSharper", "UnusedMember.Global")]
        public static double[] GetAbsError(double[] exact, double[] num)
        {
            var err = new double[exact.Length];
            for (var i = 0; i < exact.Length; ++i)
                err[i] = Math.Abs(exact[i] - num[i]);

            return err;
        }

        public static Point[] GetAbsError(Point[] exact, Point[] num)
        {
            var err = new Point[exact.Length];
            for (var i = 0; i < exact.Length; ++i)
            {
                try
                {
                    err[i] = new Point(Math.Abs(exact[i].S - num[i].S), Math.Abs(exact[i].VS - num[i].VS));
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
                err[i] = exact[i] - num[i];

            return err;
        }

        public static void Print(double[] arr, string header)
        {
            Console.WriteLine(header);
            for (var i = 0; i < arr.Length; i++)
                Console.Write(arr[i] + " ");

            Console.WriteLine();
        }

        public static void PrintAsVector(double[] arr, string header)
        {
            Console.WriteLine(header);
            for (var i = 0; i < arr.Length; i++)
                Console.Write(arr[i] + " ");

            Console.WriteLine();
        }

        public static void PrintMatrix(int n1, double[] a0, double[] b0, double[] c0, string name, bool printOnConsole = false)
        {
            using var writer = new StreamWriter(name, false);
            for (var i = 0; i < n1; i++)
            {
                int k;
                for (k = 0; k < i; k++)
                {
                    writer.Write($"{0d:E10} ");
                    if (printOnConsole)
                        Console.Write($"{0d:E10} ");
                }

                writer.Write($" {a0[k]:E10} {b0[k]:E10} {c0[k]:E10} ");
                if (printOnConsole)
                    Console.Write($" {a0[k]:E10} {b0[k]:E10} {c0[k]:E10} ");
                
                for (var j = k + 1; j < n1; j++)
                {
                    writer.Write($" {a0[j]:E10} {b0[j]:E10} {c0[j]:E10} ");
                    if (printOnConsole)
                        Console.Write($" {a0[j]:E10} {b0[j]:E10} {c0[j]:E10} ");
                }

                writer.Write('\n');
                if (printOnConsole)
                    Console.WriteLine();
            }
        }
        
        public static void PrintMatrix1(int n1, double[] a0, double[] b0, double[] c0, string name )
        {
            using var writer = new StreamWriter(name, false);
            for (var i = 0; i < n1; i++)
            {
                writer.Write($" {a0[i]:E10} {b0[i]:E10} {c0[i]:E10} ");
                writer.Write('\n'); 
            }
        }
    }
}