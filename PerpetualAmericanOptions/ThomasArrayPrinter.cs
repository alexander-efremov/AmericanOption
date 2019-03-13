using System;

namespace PerpetualAmericanOptions
{
    internal class ThomasArrayPrinter
    {
        internal void PrintThomasArrays(double[] b_t, double[] c_t, double[] d_t)
        {
            Console.WriteLine("B:");
            for (int i = 0; i < b_t.Length; i++)
            {
                Console.Write(b_t[i].ToString("E") + " ");
            }

            Console.WriteLine();
            Console.WriteLine("C:");
            for (int i = 0; i < c_t.Length; i++)
            {
                Console.Write(c_t[i].ToString("E") + " ");
            }

            Console.WriteLine();
            Console.WriteLine("D:");
            for (int i = 0; i < d_t.Length; i++)
            {
                Console.Write(d_t[i].ToString("E") + " ");
            }

            Console.WriteLine();
        }
    }
}