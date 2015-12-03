namespace Canvas.CommandLineParsing
{
    public class BaseOptions
    {
        public bool ShowHelp { get; }
        public bool ShowVersion { get; }

        public BaseOptions(bool showHelp, bool showVersion)
        {
            ShowHelp = showHelp;
            ShowVersion = showVersion;
        }
    }
}