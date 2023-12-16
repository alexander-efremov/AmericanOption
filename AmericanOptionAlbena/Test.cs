using System;
using System.IO;
using System.Linq;
using NUnit.Framework;

namespace AmericanOptionAlbena
{
    [TestFixture]
    public class Test
    {
        [Test]
        public void T()
        {
            using var writer = new StreamWriter("diff.dat", false);
            writer.WriteLine("TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA');");
            writer.WriteLine("VARIABLES = S0 t");
            writer.WriteLine($"ZONE T='diff'");
            var enumerable = File.ReadLines("1_s0_T-t_K=100_N1=10001_T=1_h_condensed_False_tau_condensed_False_finite_elem_False.dat").Skip(6);
            var file1 = enumerable.Where(s => !string.IsNullOrEmpty(s)).Select(line => double.Parse(line.Split(' ')[0]));
            var skip = File.ReadLines("1_s0_T-t_K=100_N1=10001_T=1_h_condensed_False_tau_condensed_False_finite_elem_True.dat").Skip(6);
            var file2 = skip.Where(s => !string.IsNullOrEmpty(s)).Select(line => double.Parse(line.Split(' ')[0]));
            var differences = file1.Zip(file2, (a, b) => Math.Abs(a - b));
            Console.WriteLine(string.Join(", ", differences));
        }
    }
}