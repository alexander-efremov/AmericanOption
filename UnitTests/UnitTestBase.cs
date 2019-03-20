using NUnit.Framework;

namespace PerpetualAmericanOptions
{
    public abstract class UnitTestBase
    {
        protected string WorkingDirPath;

        protected abstract string SetWorkingDir();

        [SetUp]
        protected void SetUp()
        {
            WorkingDirPath = SetWorkingDir();
        }
    }
}