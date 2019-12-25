namespace CoreLib
{
    using System.IO;

    public class TecplotPrinterSpecial : TecplotPrinter
    {
        public TecplotPrinterSpecial(int n_1, double a, double b, double tau) : base(n_1, a, b, tau)
        {
        }

        public void PrintXYSpecial(double tau, double a, double b, string filename, double t, double h1, double h2,
            double[] KS, double[] V, double S0 = 0)
        {
            var name = string.Format("{0}_hx={1}_t={2}_tau={3}_a={4}_c={5}.dat", filename, h1, t, tau, a, b);
            using (var writer = new StreamWriter(name, false))
            {
                writer.WriteLine("TITLE = 'DEM DATA'\nVARIABLES = 'x' {0}", "u");
                writer.WriteLine("ZONE T='ONE'");
                writer.WriteLine("I={0} K={1} ZONETYPE=Ordered", KS.Length, 1);
                writer.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE)");
                for (var i = 0; i < KS.Length; i++) writer.WriteLine("{0:e8}  {1:e8}", a + i * h1, KS[i]);

                writer.WriteLine("\nZONE T='TWO'");
                writer.WriteLine("I={0} K={1} ZONETYPE=Ordered", V.Length, 1);
                writer.WriteLine("DATAPACKING=POINT\nDT=(DOUBLE DOUBLE)");
                for (var i = 0; i < V.Length; i++) writer.WriteLine("{0:e8}  {1:e8}", S0 + i * h2, V[i]);
            }
        }
    }
}