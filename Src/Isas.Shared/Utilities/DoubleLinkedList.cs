using System;
using System.Collections;

namespace Illumina.SecondaryAnalysis
{
    public class DoubleLinkedList<T> : IEnumerator, IEnumerable
    {
        #region members

        public int Count;
        private DoubleLinkedListNode<T> _currentNode;
        public DoubleLinkedListNode<T> First = null;
        public DoubleLinkedListNode<T> Last = null;

        #endregion

        // constructor

        // IEnumerator and IEnumerable require these methods.
        public IEnumerator GetEnumerator()
        {
            return this;
        }

        // IEnumerator
        public bool MoveNext()
        {
            if (_currentNode == null)
            {
                _currentNode = First;
                return true;
            }

            if (_currentNode.Next != null)
            {
                _currentNode = _currentNode.Next;
                return true;
            }

            _currentNode = null;
            return false;
        }

        // IEnumerable
        public void Reset()
        {
            _currentNode = null;
            Console.WriteLine("RESET");
        }

        // IEnumerable
        public object Current
        {
            get { return _currentNode.Value; }
        }

        /// <summary>
        ///     Adds the object to the beginning of the linked list
        /// </summary>
        public void AddFirst(T value)
        {
            if (First == null)
            {
                First = new DoubleLinkedListNode<T>();
                Last = First;
                //CurrentNode = First;
            }
            else
            {
                DoubleLinkedListNode<T> temp = new DoubleLinkedListNode<T> {Next = First};
                First.Previous = temp;
                First = temp;
            }

            First.Value = value;
            ++Count;
        }

        /// <summary>
        ///     Adds the object to the end of the linked list
        /// </summary>
        public void AddLast(T value)
        {
            if (First == null)
            {
                First = new DoubleLinkedListNode<T>();
                Last = First;
                //CurrentNode = First;
            }
            else
            {
                DoubleLinkedListNode<T> temp = new DoubleLinkedListNode<T> {Previous = Last};
                Last.Next = temp;
                Last = temp;
            }

            Last.Value = value;
            ++Count;
        }

        /// <summary>
        ///     Removes all of the entries in this list
        /// </summary>
        public void Clear()
        {
            Count = 0;
            DoubleLinkedListNode<T> currentNode = First;
            while (currentNode != null)
            {
                DoubleLinkedListNode<T> tempNode = currentNode;
                currentNode = currentNode.Next;
                tempNode = null;
            }
        }

        /// <summary>
        ///     Reverses the linked list
        /// </summary>
        public void Reverse()
        {
            if (Count == 0) return;

            Last = First;
            DoubleLinkedListNode<T> start = First;

            // Loop through until null node (next node of the latest node) is found
            while (start != null)
            {
                // Swap the “Next” and “Previous” node properties
                DoubleLinkedListNode<T> temp = start.Next;
                start.Next = start.Previous;
                start.Previous = temp;

                // Head property needs to point to the latest node
                if (start.Previous == null) First = start;

                // Move on to the next node (since we just swapped “Next” and “Previous”)
                // “Next” is actually the “Previous”
                start = start.Previous;
            }
        }

        /// <summary>
        ///     Removes the first node
        /// </summary>
        /// <returns></returns>
        public void RemoveFirst()
        {
            if (Count == 0) return;

            if (Count == 1)
            {
                First = null;
            }
            else
            {
                DoubleLinkedListNode<T> temp = First.Next;
                temp.Previous = null;
                First = temp;
            }

            --Count;
        }
    }

    public class DoubleLinkedListNode<T>
    {
        public DoubleLinkedListNode<T> Next;
        public DoubleLinkedListNode<T> Previous;
        public T Value;
    }
}