using System;
using System.Diagnostics.CodeAnalysis;

namespace CoreLib
{
    [SuppressMessage("ReSharper", "UnusedType.Global", Justification = "UsedType")]
    public class ThomasArrayPrinter
    {
        // ReSharper disable once UnusedMember.Global
        public void PrintThomasArrays(double[] b, double[] c, double[] d)
        {
            Console.WriteLine("B:");
            for (var i = 0; i < b.Length; i++)
                Console.Write(b[i].ToString("E") + " ");

            Console.WriteLine("\nC:");
            for (var i = 0; i < c.Length; i++)
                Console.Write(c[i].ToString("E") + " ");

            Console.WriteLine("\nD:");
            for (var i = 0; i < d.Length; i++)
                Console.Write(d[i].ToString("E") + " ");

            Console.WriteLine();
        }
    }
}