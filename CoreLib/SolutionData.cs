using System.Collections.Generic;

namespace CoreLib
{
    public class SolutionData
    {
        public SolutionData(IReadOnlyList<double> s0, int k)
        {
            K = k;
            S0 = s0[k];
        }

        public SolutionData(double s0, int k)
        {
            K = k;
            S0 = s0;
        }

        public int K { get; private set; }

        public double S0 { get; private set; }

        public Point[] Solution { get; set; }

        public int StartPosOfS0 { get; set; }
    }
}