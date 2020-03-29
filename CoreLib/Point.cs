namespace CoreLib
{
    using System;
    using System.Diagnostics;

    [DebuggerDisplay("{ToString()}")]
    public struct Point
    {
        public Point(double s, double vs)
        {
            this.S = s;
            this.VS = vs;
        }

        public double S { get; }

        public double VS { get; }

        /// <summary>Returns the fully qualified type name of this instance.</summary>
        /// <returns>The fully qualified type name.</returns>
        public override string ToString()
        {
            return string.Format("({0:e2}, {1:e2})", S, VS);
        }
    }
}