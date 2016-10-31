using Illumina.SecondaryAnalysis;

namespace Canvas
{
    public class SmallPedigreeWorkflow
    {
        private readonly CanvasRunner _runner;

        public SmallPedigreeWorkflow(CanvasRunner runner)
        {
            _runner = runner;
        }


        public void CallPedigree(SmallPedigreeCallset callset)
        {
            _runner.CallPedigree(callset);
        }
    }
}