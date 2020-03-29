namespace CoreLib
{
    using System.Collections.Generic;

    public class SolutionData
    {
        public SolutionData(IReadOnlyList<double> S0arr, int k)
        {
            this.k = k;
            this.S0 = S0arr[k];
        }

        public SolutionData(double S0, int k)
        {
            this.k = k;
            this.S0 = S0;
        }

        public double k { get; private set; }

        public double S0 { get; private set; }

        public Point[] Solution { get; set; }

        public int StartPosOfS0 { get;  set; }
    }
}