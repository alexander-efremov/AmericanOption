namespace CoreLib
{
    using System.Collections.Generic;
    using System.Diagnostics.CodeAnalysis;
    using System.Globalization;
    using System.IO;

    public class TecplotPrinter
    {
        private readonly double a;
        private readonly double b;
        public TecplotPrinter(double a, double b)
        {
            this.a = a;
            this.b = b;
        }

        public void PrintXY(string filename, double tau, int k, double h, double[] data, double start = 0d)
        {
            var t = tau * k;
            var name = string.Format("{0}_t={3}_nx={1}_hx={2}_a={4}_c={5}.dat", filename, data.Length, h, t, this.a, this.b);
            using (var writer = new StreamWriter(name, false))
            {
                writer.WriteLine(
                    "TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA'\nVARIABLES = S {0}\nZONE T='{1}'",
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
                "{0}_id_{6}_t={3}_nx={1}_hx={2}_a={4}_c={5}.dat",
                filename,
                data.Length,
                h,
                t,
                this.a,
                this.b,
                id);
            using (var writer = new StreamWriter(name, false))
            {
                writer.WriteLine(
                    "TITLE = '{1}'\nVARIABLES = S {0}\nZONE T='{1}'",
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
                "{0}_id_{6}_t={3}_nx={1}_hx={2}_a={4}_c={5}.dat",
                filename,
                data.Length,
                h,
                t,
                this.a,
                this.b,
                id);
            using (var writer = new StreamWriter(name, false))
            {
                writer.WriteLine(
                    "TITLE = '{1}'\nVARIABLES = S {0}\nZONE T='{1}'",
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
                "{0}_id_{6}_t={3}_nx={1}_hx={2}_a={4}_c={5}.dat",
                filename,
                data.Length,
                h,
                t,
                this.a,
                this.b,
                id);
            using (var writer = new StreamWriter(name, false))
            {
                writer.WriteLine(
                    "TITLE = '{1}'\nVARIABLES = S {0}\nZONE T='{1}'",
                    "V",
                    mapName);
                int xx = 0;
                for (var i = 0; i < data.Length; i++)
                {
                    if (data[i].S>= this.b) break;
                    xx++;
                }
                writer.WriteLine("I={0} K={1} ZONETYPE=Ordered", xx, 1);
                writer.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE)");
                for (var i = 0; i < xx; i++)
                {
                    var val = data[i].VS;
                    if (double.IsInfinity(val) || val < 0d)
                    {
                        val = 0d;
                    }

                    var d = data[i].S;
                    writer.WriteLine($"{d.ToString("e8", CultureInfo.InvariantCulture)} {val.ToString("e8", CultureInfo.InvariantCulture)}");
                }
            }
        }

        public void PrintXY(string filename, int k, string t, string h, double hd,double K, Point[] data, string mapName,string mapName2, bool printexct)
        {
            // var name = string.Format(
            //     "{0}_t={3}_nx={1}_hx={2}_a={4}_c={5}.dat",
            //     filename,
            //     data.Length,
            //     h,
            //     t,
            //     this.a,
            //     this.b);
            
            var name = string.Format("{0}_k={1}.dat", filename, k);
            using (var writer = new StreamWriter(name, false))
            {
                writer.WriteLine(
                    "TITLE = '{1}'\nVARIABLES = S {0}\nZONE T='{2}'",
                    "V",
                    mapName2,
                    mapName);
                int xx = 0;
                double s = 0;
                if (printexct)
                {
                    while (s < data[0].S)
                    {
                        s += hd;
                        xx++;
                    }
                }
                
                for (var i = 0; i < data.Length; i++)
                {
                    if (data[i].S >= this.b)
                    {
                        break;
                    }
                    xx++;
                }
                
                writer.WriteLine("I={0} K={1} ZONETYPE=Ordered", xx, 1);
                writer.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE)");
                
                s = 0;
                if (printexct)
                {
                    while (s < data[0].S)
                    {
                        var val = K - s;
                        writer.WriteLine($"{s.ToString("e8", CultureInfo.InvariantCulture)} {val.ToString("e8", CultureInfo.InvariantCulture)}");
                        s += hd;
                    }
                }

                xx = 0;
                for (var i = 0; i < data.Length; i++)
                {
                    if (data[i].S>= this.b) break;
                    xx++;
                }

                for (var i = 0; i < xx; i++)
                {
                    var val = data[i].VS;
                    if (double.IsInfinity(val) || val < 0d)
                    {
                        val = 0d;
                    }

                    var d = data[i].S;
                    writer.WriteLine($"{d.ToString("e8", CultureInfo.InvariantCulture)} {val.ToString("e8", CultureInfo.InvariantCulture)}");
                }
            }
        }
        
        public void PrintXYZ(string filename, double h, double tau, List<SolutionData> data, int id, string mapName)
        {
            var name = string.Format(
                "{0}_t={3}_nx={1}_hx={2}_a={4}_c={5}.dat",
                filename,
                data[0].Solution.Length,
                h,
                this.a,
                this.b,
                id);
            using (var writer = new StreamWriter(name, false))
            {
                // each time layer is a zone
                writer.WriteLine(
                    "TITLE = '{1}'\nVARIABLES = S '{0}' '{1}'\n",
                    "t",
                    "V");
                foreach (SolutionData sol in data)
                {
                    // each time layer is a zone
                    writer.WriteLine(
                        "ZONE T='{0}'",
                        mapName + sol.k);
                    writer.WriteLine("I={0} J={1} K={2} ZONETYPE=Ordered", sol.Solution.Length, 1, 1);
                    writer.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE DOUBLE)");
                    
                    for (var i = 0; i < sol.Solution.Length; i++)
                    {
                        var val = sol.Solution[i].VS;
                        if (double.IsInfinity(val) || val < 0d)
                        {
                            val = 0d;
                        }

                        writer.WriteLine($"{sol.Solution[i].S.ToString("e8", CultureInfo.InvariantCulture)} {(sol.k * tau).ToString("e8", CultureInfo.InvariantCulture)} {val.ToString("e8", CultureInfo.InvariantCulture)}");
                    }
                }
            }
        }

        public void PrintXY(string filename, double t, double h, Point[] data, int id, string mapName, string comment, int commentIndex)
        {
            var name = string.Format(
                "{0}_id_{6}_t={3}_nx={1}_hx={2}_a={4}_c={5}.dat",
                filename,
                data.Length,
                h,
                t,
                this.a,
                this.b,
                id);
            using (var writer = new StreamWriter(name, false))
            {
                writer.WriteLine(
                    "TITLE = '{1}'\nVARIABLES = S {0}\nZONE T='{1}'",
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
            var name = string.Format("{0}_nx={1}_hx={2}_t={3}_a={4}_c={5}.dat", filename, exact.Length, h, t, this.a, this.b);
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