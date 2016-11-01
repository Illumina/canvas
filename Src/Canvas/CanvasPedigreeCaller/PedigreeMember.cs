using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using CanvasCommon;

namespace CanvasPedigreeCaller
{
    public class PedigreeMember
    {
        public List<CanvasSegment> Segments = new List<CanvasSegment>();
        public double MeanCoverage { get; set; }
        public double MeanMafCoverage { get; set; }
        public string Name { get; set; }
        public List<string> Parents { get; set; }
        public List<string> Offspring { get; set; }
        public PloidyInfo Ploidy { get; set; }
        public double Variance { get; internal set; }
        public CopyNumberModel CnModel { get; set; }

        public double GetCoverage(int segmentIndex)
        {
            return Segments[segmentIndex].Counts.Average();
        }
        public Tuple<int,int> GetAlleleCounts(int segmentIndex)
        {
            double allele1 = Segments[segmentIndex].VariantAlleleCounts.Select(x=>x.Item1).Average();
            double allele2 = Segments[segmentIndex].VariantAlleleCounts.Select(x => x.Item2).Average();
            return new Tuple<int, int>(Convert.ToInt32(allele1), Convert.ToInt32(allele2));
        }
    }

    // requirements :
    // iterate over pedigree
    // get Root parents
    // get Parents 
    // get Children
    // check isRoot
    // check isLeaf
    // return tree by names 


    public class Pedigree<TItem>: IEnumerable<TItem>
    {

        public TItem NodeData { get; set; }
        public Pedigree<TItem> Partner { get; set; }

        public bool IsRoot;

        public List<TItem> Children { get; set; }

        public Pedigree(TItem value, bool isRoot = false)

        {

            this.NodeData = value;

            this.Partner = null;

            this.Children = null;

            this.IsRoot = isRoot;

        }

        public void AddParent(TItem value)
        {
            if (this.Partner == null)
                this.Partner = new Pedigree<TItem>(value, true);
            else
                throw new ArgumentException("Could not have nore than two parents");
        }


        public void AddChildren(TItem newItem)
        {
            if (this.Children == null)
                this.Children = new List<TItem>();
                this.Children.Add(newItem);
        }

        public List<TItem> GetChildren()
        {
            return this.Children;
        }

        public List<TItem> GetParents()
        {
            return new List<TItem>{this.NodeData, this.Partner.NodeData};
        }

        IEnumerator<TItem> IEnumerable<TItem>.GetEnumerator()
        {
           yield return this.NodeData;

           yield return this.Partner.NodeData;

           foreach (TItem item in this.Children)
               yield return item;
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            throw new NotImplementedException();
        }

    }
}