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
        public string Name { get; set; }
        public List<string> Parents { get; set; }
        public List<string> Offspring { get; set; }

        public PloidyInfo Ploidy { get; set; }

        public PedigreeMember()
        {
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


    public class PedigreeNode<TItem>: IEnumerable<TItem>
    {

        public TItem NodeData { get; set; }
        public string NodeName { get; set; }

        public bool IsRoot;

        public List<PedigreeNode<TItem>> Parents { get; set; }

        public List<PedigreeNode<TItem>> Children { get; set; }

        public List<PedigreeNode<TItem>> Siblings { get; set; }


        public PedigreeNode(TItem value, string name, bool isRoot = false)

        {

            this.NodeData = value;

            this.NodeName = name;

            this.Parents = null;

            this.Children = null;

            this.Siblings = null;

            this.IsRoot = isRoot;

        }

        public void AddParent(TItem newItem, string newName)

        {
            if (this.IsRoot)
                throw new ArgumentException($"Could not add parents to root node {this.NodeName}");
            if (this.Parents.Count > 2)
                throw new ArgumentException($"Node {this.NodeName} could not have more than two parents");

            this.Parents.Add(new PedigreeNode<TItem>(newItem, newName));
        }

        public void AddSiblings(TItem newItem, string newName)
        {
            if (this.IsRoot)
                this.Siblings.Add(new PedigreeNode<TItem>(newItem, newName, true));
            else
                this.Siblings.Add(new PedigreeNode<TItem>(newItem, newName, false));
        }

        public void AddChildren(TItem newItem, string newName)
        {
            this.Children.Add(new PedigreeNode<TItem>(newItem, newName));
        }

        IEnumerator<TItem> IEnumerable<TItem>.GetEnumerator()
        {
            if (this.Siblings != null)
            {
                foreach (TItem item in this.Siblings)
                {
                    yield return item;
                }
            }

        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            throw new NotImplementedException();
        }
    }
}