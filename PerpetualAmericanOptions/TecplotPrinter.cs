using System.IO;

namespace PerpetualAmericanOptions
{
    internal class TecplotPrinter
    {
        private readonly int n_1;
        private readonly double h;
        private readonly double a;
        private readonly double b;
        private readonly double tau;

        internal TecplotPrinter(int n_1, double h, double a, double b, double tau)
        {
            this.n_1 = n_1;
            this.h = h;
            this.a = a;
            this.b = b;
            this.tau = tau;
        }

        internal void PrintXY(string filename,
            double t, double[] data)
        {
            var name = string.Format("{0}_nx={1}_hx={2}_t={3}_tau={4}_a={5}_c={6}.dat", filename, n_1, h, t, tau, a, b);
            using (var writer = new StreamWriter(name, false))
            {
                writer.WriteLine(
                    "TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA'\nVARIABLES = 'x' {0}\nZONE T='{1}'", filename,
                    filename);
                writer.WriteLine("\nI={0} K={1} ZONETYPE=Ordered", n_1, 1);
                writer.WriteLine("\nDATAPACKING=POINT\nDT=(SINGLE SINGLE)");
                for (var i = 0; i < n_1; i++)
                {
                    writer.WriteLine("\n{0:e2} {1:e2}", i * h, data[i]);
                }
            }
        }

        internal void PrintXY(string filename,
            double t,
            double[] exact, double[] numerical)
        {
            var name = string.Format("{0}_nx={1}_hx={2}_t={3}_tau={4}_a={5}_c={6}.dat", filename, n_1, h, t, tau, a, b);
            using (var writer = new StreamWriter(name, false))
            {
                writer.WriteLine("TITLE = 'DEM DATA'\nVARIABLES = 'x' {0}", filename);
                writer.WriteLine("\nZONE T='EXACT'");
                writer.WriteLine("\nI={0} K={1} ZONETYPE=Ordered", n_1, 1);
                writer.WriteLine("\nDATAPACKING=POINT\nDT=(SINGLE SINGLE)");
                for (var i = 0; i < n_1; i++)
                {
                    writer.WriteLine("\n{0:e2}  {1:e2}", i * h, exact[i]);
                }

                writer.WriteLine("\n\nZONE T='NUMER'");
                writer.WriteLine("\nI={0} K={1} ZONETYPE=Ordered", n_1, 1);
                writer.WriteLine("\nDATAPACKING=POINT\nDT=(SINGLE SINGLE)");
                for (var i = 0; i < n_1; i++)
                {
                    writer.WriteLine("\n{0:e2}  {1:e2}", i * h, numerical[i]);
                }
            }
        }
    }
}