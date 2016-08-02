using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CanvasCommon
{
    /// <summary>
    /// http://allantech.blogspot.com/2007/02/c-and-keeping-your-lists-sorted.html
    /// </summary>
    /// <typeparam name="T"></typeparam>
    public class SortedList<T> : List<T>
        where T : IComparable<T>
    {

        public new void Add(T Item)
        {
            if (Count == 0)
            {
                //No list items
                base.Add(Item);
                return;
            }
            if (Item.CompareTo(this[Count - 1]) > 0)
            {
                //Bigger than Max
                base.Add(Item);
                return;
            }
            int min = 0;
            int max = Count - 1;
            while ((max - min) > 1)
            {
                //Find half point
                int half = min + ((max - min) / 2);
                //Compare if it's bigger or smaller than the current item.
                int comp = Item.CompareTo(this[half]);
                if (comp == 0)
                {
                    //Item is equal to half point
                    Insert(half, Item);
                    return;
                }
                else if (comp < 0) max = half;   //Item is smaller
                else min = half;   //Item is bigger
            }
            if (Item.CompareTo(this[min]) <= 0) Insert(min, Item);
            else Insert(min + 1, Item);
        }
    }
}
