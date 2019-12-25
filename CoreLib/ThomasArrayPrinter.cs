namespace CoreLib
{
    using System;

    public class ThomasArrayPrinter
    {
        // ReSharper disable once UnusedMember.Global
        public void PrintThomasArrays(double[] b_t, double[] c_t, double[] d_t)
        {
            Console.WriteLine("B:");
            for (var i = 0; i < b_t.Length; i++)
            {
                Console.Write(b_t[i].ToString("E") + " ");
            }

            Console.WriteLine();
            Console.WriteLine("C:");
            for (var i = 0; i < c_t.Length; i++)
            {
                Console.Write(c_t[i].ToString("E") + " ");
            }

            Console.WriteLine();
            Console.WriteLine("D:");
            for (var i = 0; i < d_t.Length; i++)
            {
                Console.Write(d_t[i].ToString("E") + " ");
            }

            Console.WriteLine();
        }
    }
}