namespace CoreLib
{
    using System.Diagnostics.CodeAnalysis;
    using System.Globalization;
    using System.IO;

    public class TecplotPrinter
    {
        private readonly double a;
        private readonly double b;
        private readonly double tau;

        public TecplotPrinter(double a, double b, double tau)
        {
            this.a = a;
            this.b = b;
            this.tau = tau;
        }

        public void PrintXY(string filename, double t, double h, double[] data, double start = 0d)
        {
            var name = string.Format("{0}_t={3}_nx={1}_hx={2}_tau={4}_a={5}_c={6}.dat", filename, data.Length, h, t, this.tau, this.a, this.b);
            using (var writer = new StreamWriter(name, false))
            {
                writer.WriteLine(
                    "TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA'\nVARIABLES = 'S' {0}\nZONE T='{1}'",
                    "V",
                    "SubZone");
                writer.WriteLine("I={0} K={1} ZONETYPE=Ordered", data.Length, 1);
                writer.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE)");
                for (var i = 0; i < data.Length; i++)
                {
                    writer.WriteLine($"{(start + i * h).ToString("e8", CultureInfo.InvariantCulture)} {data[i].ToString("e8", CultureInfo.InvariantCulture)}");
                }
            }
        }

        [SuppressMessage("ReSharper", "UnusedMember.Global")]
        public void PrintXY(string filename, double t, double h, double[] data, int id, string mapName)
        {
            var name = string.Format(
                "{0}_id_{7}_t={3}_nx={1}_hx={2}_tau={4}_a={5}_c={6}.dat",
                filename,
                data.Length,
                h,
                t,
                this.tau,
                this.a,
                this.b,
                id);
            using (var writer = new StreamWriter(name, false))
            {
                writer.WriteLine(
                    "TITLE = '{1}'\nVARIABLES = 'S' {0}\nZONE T='{1}'",
                    "V",
                    mapName);
                writer.WriteLine("I={0} K={1} ZONETYPE=Ordered", data.Length, 1);
                writer.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE)");
                for (var i = 0; i < data.Length; i++)
                {
                    var val = data[i];
                    if (double.IsInfinity(val) || val < 0d)
                    {
                        val = 0d;
                    }

                    writer.WriteLine($"{(i * h).ToString("e8", CultureInfo.InvariantCulture)} {val.ToString("e8", CultureInfo.InvariantCulture)}");
                }
            }
        }

        [SuppressMessage("ReSharper", "UnusedMember.Global")]
        public void PrintXY(string filename, double t, double h, double[] data, int id, string mapName, string comment, int commentIndex)
        {
            var name = string.Format(
                "{0}_id_{7}_t={3}_nx={1}_hx={2}_tau={4}_a={5}_c={6}.dat",
                filename,
                data.Length,
                h,
                t,
                this.tau,
                this.a,
                this.b,
                id);
            using (var writer = new StreamWriter(name, false))
            {
                writer.WriteLine(
                    "TITLE = '{1}'\nVARIABLES = 'S' {0}\nZONE T='{1}'",
                    "V",
                    mapName);
                writer.WriteLine("I={0} K={1} ZONETYPE=Ordered", data.Length, 1);
                writer.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE)");
                for (var i = 0; i < data.Length; i++)
                {
                    var val = data[i];
                    if (double.IsInfinity(val) || val < 0d)
                    {
                        val = 0d;
                    }

                    if (i == commentIndex)
                    {
                        writer.WriteLine($"#{comment}");
                    }

                    writer.WriteLine($"{(i * h).ToString("e8", CultureInfo.InvariantCulture)} {val.ToString("e8", CultureInfo.InvariantCulture)}");
                }
            }
        }

        public void PrintXY(string filename, double t, double h, Point[] data, int id, string mapName)
        {
            var name = string.Format(
                "{0}_id_{7}_t={3}_nx={1}_hx={2}_tau={4}_a={5}_c={6}.dat",
                filename,
                data.Length,
                h,
                t,
                this.tau,
                this.a,
                this.b,
                id);
            using (var writer = new StreamWriter(name, false))
            {
                writer.WriteLine(
                    "TITLE = '{1}'\nVARIABLES = 'S' {0}\nZONE T='{1}'",
                    "V",
                    mapName);
                writer.WriteLine("I={0} K={1} ZONETYPE=Ordered", data.Length, 1);
                writer.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE)");
                for (var i = 0; i < data.Length; i++)
                {
                    var val = data[i].VS;
                    if (double.IsInfinity(val) || val < 0d)
                    {
                        val = 0d;
                    }

                    writer.WriteLine($"{data[i].S.ToString("e8", CultureInfo.InvariantCulture)} {val.ToString("e8", CultureInfo.InvariantCulture)}");
                }
            }
        }

        public void PrintXY(string filename, double t, double h, Point[] data, int id, string mapName, string comment, int commentIndex)
        {
            var name = string.Format(
                "{0}_id_{7}_t={3}_nx={1}_hx={2}_tau={4}_a={5}_c={6}.dat",
                filename,
                data.Length,
                h,
                t,
                this.tau,
                this.a,
                this.b,
                id);
            using (var writer = new StreamWriter(name, false))
            {
                writer.WriteLine(
                    "TITLE = '{1}'\nVARIABLES = 'S' {0}\nZONE T='{1}'",
                    "V",
                    mapName);
                writer.WriteLine("I={0} K={1} ZONETYPE=Ordered", data.Length, 1);
                writer.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE)");
                for (var i = 0; i < data.Length; i++)
                {
                    var val = data[i].VS;
                    if (double.IsInfinity(val) || val < 0d)
                    {
                        val = 0d;
                    }

                    if (i == commentIndex)
                    {
                        writer.WriteLine($"#{comment}");
                    }

                    writer.WriteLine($"{data[i].S.ToString("e8", CultureInfo.InvariantCulture)} {val.ToString("e8", CultureInfo.InvariantCulture)}");
                }
            }
        }

        public void PrintXY(string filename, double t, double h, double[] exact, double[] numeric, double S0)
        {
            var name = string.Format("{0}_nx={1}_hx={2}_t={3}_tau={4}_a={5}_c={6}.dat", filename, exact.Length, h, t, this.tau, this.a, this.b);
            using (var writer = new StreamWriter(name, false))
            {
                writer.WriteLine("TITLE = 'DEM DATA'\nVARIABLES = 'x' {0}", "u");
                writer.WriteLine("ZONE T='ONE'");
                writer.WriteLine("I={0} K={1} ZONETYPE=Ordered", exact.Length, 1);
                writer.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE)");
                for (var i = 0; i < exact.Length; i++)
                {
                    writer.WriteLine($"{(S0 + i * h).ToString("e8", CultureInfo.InvariantCulture)}  {exact[i].ToString("e8", CultureInfo.InvariantCulture)}");
                }

                writer.WriteLine("\nZONE T='TWO'");
                writer.WriteLine("I={0} K={1} ZONETYPE=Ordered", exact.Length, 1);
                writer.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE)");
                for (var i = 0; i < exact.Length; i++)
                {
                    writer.WriteLine($"{(S0 + i * h).ToString("e8", CultureInfo.InvariantCulture)}  {numeric[i].ToString("e8", CultureInfo.InvariantCulture)}");
                }
            }
        }
    }
}