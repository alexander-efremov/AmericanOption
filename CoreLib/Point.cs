namespace CoreLib
{
    public struct Point
    {
        public Point(double s, double vs)
        {
            this.S = s;
            this.VS = vs;
        }

        public double S { get; }

        public double VS { get; }
    }
}