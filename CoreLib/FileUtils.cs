namespace CoreLib.Utils
{
    using System.IO;
    using System.Text;

    public static class FileUtils
    {
        public static void WriteVectorToFile(string path, double[] arr)
        {
            var sb = new StringBuilder();
            for (int i = 0; i < arr.Length; i++)
            {
                sb.AppendLine(arr[i].ToString("e2"));
            }
            File.WriteAllText(path, sb.ToString());     
        }
    }
}