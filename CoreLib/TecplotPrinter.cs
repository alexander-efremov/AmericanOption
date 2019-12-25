namespace CoreLib
{
    using System.IO;

    public class TecplotPrinter
    {
        private readonly double a;
        private readonly double b;
        private readonly int n1;
        private readonly double tau;

        public TecplotPrinter(int n_1, double a, double b, double tau)
        {
            this.n1 = n_1;
            this.a = a;
            this.b = b;
            this.tau = tau;
        }

        public void PrintXY(string filename, double t, double h, double[] data, double start = 0d)
        {
            var name = string.Format("{0}_nx={1}_hx={2}_t={3}_tau={4}_a={5}_c={6}.dat", filename, this.n1, h, t, this.tau, this.a, this.b);
            using (var writer = new StreamWriter(name, false))
            {
                writer.WriteLine("TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA'\nVARIABLES = 'S' {0}\nZONE T='{1}'", "V", "SubZone");
                writer.WriteLine("I={0} K={1} ZONETYPE=Ordered", this.n1, 1);
                writer.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE)");
                for (var i = 0; i < this.n1; i++)
                {
                    writer.WriteLine("{0:e8} {1:e8}", start + i * h, data[i]);
                }
            }
        }

        public void PrintXY(string filename, double t, double h, double[] exact, double[] numerical, double S0 = 0)
        {
            var name = string.Format("{0}_nx={1}_hx={2}_t={3}_tau={4}_a={5}_c={6}.dat", filename, this.n1, h, t, this.tau, this.a, this.b);
            using (var writer = new StreamWriter(name, false))
            {
                writer.WriteLine("TITLE = 'DEM DATA'\nVARIABLES = 'x' {0}", "u");
                writer.WriteLine("ZONE T='ONE'");
                writer.WriteLine("I={0} K={1} ZONETYPE=Ordered", this.n1, 1);
                writer.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE)");
                for (var i = 0; i < this.n1; i++)
                {
                    writer.WriteLine("{0:e8}  {1:e8}", S0 + i * h, exact[i]);
                }

                writer.WriteLine("\nZONE T='TWO'");
                writer.WriteLine("I={0} K={1} ZONETYPE=Ordered", this.n1, 1);
                writer.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE)");
                for (var i = 0; i < this.n1; i++)
                {
                    writer.WriteLine("{0:e8}  {1:e8}", S0 + i * h, numerical[i]);
                }
            }
        }
    }
}